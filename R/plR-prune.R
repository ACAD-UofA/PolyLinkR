#' @title Gene set enrichment adjusting for multiple testing and shared genes
#'
#' @description
#' Adjusts gene set scores to account dependencies between gene sets using a
#' pruning routine involving sequential p-value ranking and removal of shared
#' genes. Performs multiple testing correction on the pruned gene set p-values.
#' Note that this step require that at more than one gene set with shared genes
#' is present, otherwise the function will exit with an error message.
#'
#' @param plR.input
#'   \code{plR} class object; output from \code{polylinkR::plR_permute} or
#'   \code{polylinkR::plR_rescale}. Required.
#'
#' @param n.fdr
#'   \code{integer}; number of times to replicate the pruning procedure.
#'   Used to estimate FDR-corrected p (q) values using a histogram method.
#'   Defaults to \code{300L}. Range \code{[100, Inf)} and must be exactly
#'   divisible by \code{100}.
#'
#' @param est.pi0
#'   \code{logical}; should \code{pi0} be estimated during FDR correction?
#'   Defaults to \code{TRUE}. If \code{FALSE}, \code{pi0} is set to \code{1}.
#'
#' @param tolerance
#'   \code{numeric}; minimum difference between successive \code{pi0}
#'   estimates for estimation to stop. Defaults to \code{1e-3}. Range
#'   \code{(0, 0.01]}.
#'
#' @param verbose
#'   \code{logical}; should progress reporting be enabled? Defaults to
#'   \code{TRUE}.
#'
#' @param n.cores
#'   \code{integer}; number of cores for parallel processing. Defaults to
#'   \code{1} or \code{maximum cores - 1}. Must be in \code{[1, maximum cores]}.
#'
#' @param fut.plan
#'   \code{character}; parallel backend from the \code{future} package.
#'   Defaults to user \code{n.cores} choice or checks cores, choosing
#'   \code{"sequential"} on single-core and \code{"multisession"} on
#'   multi-core systems. Options: \code{"multisession"},
#'   \code{"multicore"}, \code{"cluster"}, \code{"sequential"}.
#'
#' @export
#'
#' @import data.table
#' @import foreach
#' @import doFuture
#' @importFrom future plan
#' @importFrom progressr handlers progressor with_progress
#' @importFrom dqrng dqset.seed dqsample.int generateSeedVectors
#' @importFrom Rfast rowMins
#' @importFrom rlang local_options
#' @importFrom cli boxx col_cyan style_italic
#'
#' @return
#' A \code{plR} S3 object containing:
#' \describe{
#'   \item{\code{obj.info}}{\code{data.table} for each gene (object).}
#'   \item{\code{set.info}}{\code{data.table} for each gene set.}
#'   \item{\code{set.obj}}{\code{data.table} mapping genes to gene sets.}
#' }
#'
#' Pruned set scores (i.e., corrected for between–gene-set correlations
#' caused by shared genes) are recorded in \code{set.info} as
#' \code{setScore.pr}. Corrected scores are shown in \code{setScore.pr.p}.
#' Results from multiple testing correction are shown in \code{setScore.pr.q}.
#'
#' All \code{plR} S3 objects include auxiliary information and summaries as
#' attributes:
#' \describe{
#'   \item{\code{plr.data}}{
#'     Reusable datasets and parameters, including GPD fitting results and
#'     autocovariance estimation.
#'   }
#'   \item{\code{plr.args}}{
#'     Argument settings used in the function.
#'   }
#'   \item{\code{plr.summary}}{
#'     Model fitting results and associated data used in diagnostics and the
#'     \code{plot} method.
#'   }
#'   \item{\code{plr.session}}{
#'     R session information and function run time.
#'   }
#'   \item{\code{plr.track}}{
#'     Internal tracking number indicating functional steps.
#'   }
#'   \item{\code{class}}{
#'     S3 object class (i.e., \code{plr}).
#'   }
#' }
#'
#' Each attribute can be accessed using \code{attr()} or \code{attributes()}.
#' The first four attribute classes aggregate information over successive
#' functions. For example, to access the \code{plr.data} attribute for the
#' \code{plR} object output after running \code{plR_permute}, use
#' \code{attributes(X)$plR.data$permute.data}, where \code{X} is the object
#' name. Similarly, the arguments used in \code{plR_read} are in
#' \code{attributes(X)$plR.args$read.args}.
#'
#' The primary data structure of the \code{plR} object can be accessed using
#' \code{print()} or by simply typing the object's name. Diagnostic plots of
#' the various data transformations and enrichment results are available via
#' the \code{plot()} method.
#'
#' @examples
#' \dontrun{
#' # Assuming `my_plr` is the result of `polylinkR::plR_permute` or
#' # `polylinkR::plR_rescale`
#'
#' # Example 1: Basic usage
#' new_plr <- plR_prune(plR.input = my_plr)
#'
#' # Example 2: Increase iterations for null p-value distribution
#' new_plr <- plR_prune(
#'    plR.input = my_plr,
#'    n.fdr = 1000L
#' )
#'
#' # Example 3: Do not estimate pi0 (sets pi0 = 1)
#' new_plr <- plR_prune(
#'    plR.input = my_plr,
#'    est.pi0 = FALSE
#' )
#' }
plR_prune <- function(plR.input, n.fdr = 300L, est.pi0 = TRUE, tolerance = 1e-3,
                      verbose = TRUE, n.cores = "auto", fut.plan = "auto") {

   ##=========================================================================##
   ##PART 1: clean data and run checks
   ##=========================================================================##

   # track function run time
   startT <- Sys.time()

   # control message reporting
   rlang::local_options(verbose = verbose, .frame = environment())

   # perform checks and unpack required plR objects
   .plR_check(f = "prune", ENV = environment())
   .arg_check(f = "prune", ENV = environment())
   .plR_unpack(plr = plR.input, ENV = environment())
   rm(plR.input); gc(verbose = FALSE)

   # set random seed
   seed <- plr.seed$seed
   set.seed(seed) # reinstate seed from plR_permute

   # ensure contiguous IDs for genes and sets
   .file_set(OI = obj.info, SI = set.info, SO = set.obj, pos.info = pos.info,
             ENV = environment())

   # set up for parallel back end and reporting
   .par_params(n.cores = n.cores, fut.plan = fut.plan, verbose = verbose,
               ENV = environment())
   progressr::handlers(prog.hand) # set up parallel progress reporting
   oplan <- future::plan() # obtain default future plan
   on.exit(future::plan(oplan), add = TRUE) # reset future options upon function completion
   if (fut.plan == "sequential") {
      future::plan(strategy = fut.plan)
   } else {
      future::plan(strategy = fut.plan, workers = n.cores)
      rlang::local_options(future.globals.maxSize = 1 * 1e9, # increase max foreach output
                           future.messages = FALSE, # stop printing package loading details during parallelised loops
                           doFuture.rng.onMisuse = "ignore", # ignore false positive RNG calls
                           .frame = environment())
   }

   # checks completed, announce function
   hdr <- "Running plR_prune: FDR correction accounting for shared genes"
   pdg <- (80 - max(nchar(hdr))) / 2
   .vrb(cli::boxx(hdr, padding = c(0, floor(pdg), 0, ceiling(pdg)),
                  align = "center", col = "cyan", border_col = "cyan"))

   # report collated warnings and messages
   .vrb("\n")
   if (!is.null(WMSS)) warning(WMSS, immediate. = T, call. = F)
   if (!is.null(MSS)) .vrb(paste0("\n", MSS, "\n"))
   if (!is.null(pWMSS)) warning(pWMSS, immediate. = T, call. = F)
   if (!is.null(pMSS)) .vrb(paste0("\n", pMSS, "\n"))
   .vrb("\n")

   ##=========================================================================##
   ##PART 2: generate FDR sets and empirical prepare objects for pruning step
   ##=========================================================================##

   .vrb(cli::style_italic("Preparing objects for pruning step:\n"))

   # generate random scores for q-value adjustment
   .vrb(paste("Generating", n.fdr, "permuted gene sets for pruning step\n"))
   # create gene set by gene identity matrix (genes limited to those in sets)
   dqi <- dqrng::generateSeedVectors(n.fdr)
   os0 <- obj.info$objStat.std
   n.max <- max(set.info$setN)
   soi <- set.obj[, .(sID = setID, oID = as.integer(factor(objID)))] # make gene IDs in gene sets contiguous
   data.table::setkey(soi, sID, oID)
   fopt <- list(packages = "data.table", globals = c("soi", "prog"))

   fdr.perm <- foreach::foreach(i = dqi) %do% {
      dqrng::dqset.seed(i)
      ps <- dqrng::dqsample.int(n = n.genes, size = n.set.genes) # permute data
      ss <- soi[, sum(os0[ps[oID]]), keyby = sID]$V1 # unscaled gene set score
      list(ss, ps)
   }

   FDR <- lapply(fdr.perm, "[[", 1) # extract standardised scores
   fdr.perm <- lapply(fdr.perm, "[[", 2) # extract permuted gene IDs

   .vrb("Calculating p values for pre-pruned gene sets\n")
   # create FDR summary object
   FDR.pr <- list(FDR.rep = rep(1:n.fdr, each = n.sets),
                  setID = rep(1:n.fdr, n.sets),
                  setN = rep(set.info$setN, n.fdr))
   data.table::setDT(FDR.pr)
   fpc <- (n.genes - 1:n.max) / (n.genes - 1) # finite population correction
   sN <- set.info$setN

   # if rescaling used, rescale FDR gene set scores
   rescaled <- plr.args$rescale.args$rescale
   if (is.null(rescaled)) rescaled <- FALSE
   if (rescaled) { # include covariance in standardised scores
      .vrb(paste("Rescaling pruning gene sets\n"))
      ac0 <- ac[, .(A = objID.A, B = objID.B, CV = 2 * rho)]
      data.table::setkey(ac0, A, B)
      c.i <- rep(0, n.sets)

      FDR.pr[[4]] <- foreach::foreach(f.i = fdr.perm, F.i = FDR, .combine = c) %do% {
         ss.i <- soi[, .(sID, oID = f.i[oID])]
         cov.i <- .get_cov(SSi = ss.i, ac0 = ac0)
         sN0 <- sN[cov.i$sID]
         c.i[cov.i$sID] <- cov.i$CV
         F.i / sqrt((sN + c.i) * fpc[sN])
      }
   } else { # standardised scores
      FDR.pr[[4]] <- foreach::foreach(F.i = FDR, .combine = c) %do% {
         F.i / sqrt(sN * fpc[sN])
      }
      ac0 <- NULL
   }

   # calculate p values
   n.perm <- plr.args$permute.args$n.perm
   n.min <- min(plr.args$read.args$min.set.n)
   gpd.cutoff <- plr.args$permute.args$gpd.cutoff
   alt <- plr.args$permute.args$alt

   if (rescaled) {
      ecdf0 <- ecdf.rs
      gpd0 <- gpd.rs
      rm(gpd.rs); gc(verbose = FALSE)
   } else {
      ecdf0 <- ecdf.std
      gpd0 <- gpd.std
      rm(gpd.std); gc(verbose = FALSE)
   }

   if (alt == "lower") {
      .get_p <- .get_p_lt
   } else {
      .get_p <- .get_p_ut
   }

   # compute p values
   FDR.pr[[5]] <- NA_real_
   for(j in unique(sN)) {
      fR <- which(FDR.pr$setN == j)
      pI <- .get_p(ss = FDR.pr[[4]][fR], n = j, g0 = gpd0, e0 = ecdf0)
      FDR.pr[[5]][fR] <- pI
   }
   data.table::setDT(FDR.pr)
   nms <- c(names(FDR.pr)[1:3],
            paste0("setScore.", ifelse(rescaled, "rs", "std"), c("", ".p")))
   data.table::setnames(FDR.pr, nms)

   .vrb("Determining shared genes for each pair of gene sets\n")
   progressr::with_progress({
      prog <- progressr::progressor(steps = n.sets)
      suppressPackageStartupMessages(
         PI0 <- foreach::foreach(i = 1:n.sets, .options.future = fopt,
                                 .combine = rbind) %dofuture% {
            Xi <- soi[.(i)]$oID
            out <- soi[.(Xi), on = "oID", nomatch = 0]
            prog()
            cbind(A0 = i, out[, .(sID, oID)])
         }
      )
   })
   data.table::setnames(PI0, c("A", "B", "X"))
   data.table::setkey(PI0, A, B)

   ##=========================================================================##
   ##PART 3: Prune scores and use permuted datasets to estimate FDR
   ##=========================================================================##

   .vrb(cli::style_italic("\nRunning pruning step\n"))
   .vrb("Pruning observed gene set scores\n")

   # limit autocovariance matrix to genes appearing in gene sets
   objName <- paste0("setScore.", ifelse(rescaled, "rs", "std"), ".p")
   Oi <- data.table::data.table(I = set.genes, J = 1:n.set.genes)
   data.table::setkey(Oi, I)
   if (rescaled) {
      acX <- ac0[.(Oi), nomatch = 0][.(Oi), on = .(B = I), nomatch = 0]
      acX <- acX[, .(A = J, B = i.J, CV)]
      data.table::setkey(acX, A, B)
   } else {
      acX <- NULL
   }

   SSO.pr <- .prune_sets(OBS = set.info$setScore.std * sqrt(sN * fpc[sN]), # get raw unscaled scores
                         EXP = ecdf0, GPD = gpd0, acX = acX, nX = sN,
                         osX = obj.info[set.genes]$objStat.std, # limit scores to genes in gene sets
                         pX = set.info[, get(objName)], mss = n.min, fpc = fpc,
                         PI = data.table::copy(PI0), gP = .get_p)

   SSO.pr <- cbind(1:n.sets, SSO.pr) # add set ID
   SSO.pr <- data.table::as.data.table(SSO.pr)
   k.nms <- c("setID", "Rank", "Tot.obj.rem", "Tot.set.rem", "setN.pr",
              "setScore.pr", "setScore.pr.p")
   data.table::setnames(SSO.pr, k.nms)
   SSO.pr[Rank == 0, (k.nms[-1]) := list(NA, NA, NA, NA, NA, NA)] # set 0 to NA

   # FDR set scores
   .vrb("Pruning FDR gene set scores\n")
   fopt <- list(packages = c("data.table", "Rfast", 'fExtremes'),
                globals = c("os0", "ac0", "sN", "ecdf0", "gpd0", "n.min",
                            "PI0", "fit.meth", "fpc", "prog", "n.set.genes",
                            ".prune_sets", ".get_cov", ".get_p"))

   progressr::with_progress({
      prog <- progressr::progressor(steps = n.fdr)
      suppressPackageStartupMessages(
         ps <- foreach::foreach(f.i = fdr.perm, ss.i = FDR,
                                p.i = split(FDR.pr$setScore.rs.p,
                                            rep(1:n.fdr, each = n.sets)),
                                .options.future = fopt) %dofuture% {
            os.i <- os0[f.i] # extract permuted gene scores, retain permuted order
            if (!is.null(ac0)) {
              Oi <- data.table::data.table(I = f.i, J = 1:n.set.genes)
              data.table::setkey(Oi, I)
              acX <- ac0[.(Oi), nomatch = 0][.(Oi), on = .(B = I), nomatch = 0]
              acX <- acX[, .(A = J, B = i.J, CV)]
              data.table::setkey(acX, A, B)
            } else {
              acX <- NULL
            }

            ps.i <- .prune_sets(OBS = ss.i, EXP = ecdf0, GPD = gpd0, acX = acX,
                                pX = p.i, nX = sN, osX = os.i, PI = PI0,
                                mss = n.min, fpc = fpc, gP = .get_p)
            prog()
            ps.i
         }
      )
   })

   # add observed pruned results
   FDR.pr <- cbind(FDR.pr, data.table::as.data.table(do.call(rbind, ps)))
   p.nms <- c("Rank", "n.obj.rem", "n.set.rem",  "setN.pr",
              "setScore.pr", "setScore.pr.p")
   data.table::setnames(FDR.pr, c(nms, p.nms))
   FDR.pr[Rank == 0, (p.nms) := list(NA, NA, NA, NA, NA, NA)]

   # rename gpd0
   assign(x = paste0("gpd.", ifelse(rescaled, "rs", "std")), value = gpd0) # reassign gpd object name
   rm(ps, gpd0); gc(verbose = FALSE)

   ##=========================================================================##
   ##PART 4: Compute adjusted p-values
   ##=========================================================================##

   .vrb("Computing FDR corrected p values [q values]\n")

   #---------------------------------------------------------------------------#
   ## estimate pi0 using histogram method
   #---------------------------------------------------------------------------#

   n.fdr.sets <- nrow(SSO.pr[!is.na(Rank)]) # number of pruned p-values evaluated
   if (est.pi0) { # estimate pi0
      .est_pi0(n.fdr.sets = n.fdr.sets, SSO.pr = SSO.pr, FDR.pr = FDR.pr,
               n.fdr = n.fdr, tolerance = tolerance, env = environment())
   } else {
      pi0 <- 1
   }

   #---------------------------------------------------------------------------#
   ## calculate q values
   #---------------------------------------------------------------------------#

   # bins to evaluate expected number of true nulls
   emp.p.bins <- c(0, sort(unique(SSO.pr$setScore.pr.p)))
   q.dt <- .fdr_bin_counts(O = SSO.pr[!is.na(Rank)],
                           E = FDR.pr[!is.na(Rank)],
                           p.bins = emp.p.bins, n.fdr = n.fdr)
   q.dt[, pi0 := pi0]
   q.dt[, VP := cumsum(EXP)]
   q.dt[, RP := cumsum(OBS)]
   q.dt[, FDR := Rfast::rowMins(cbind(1, pi0 * VP / RP), value = TRUE)]
   q.dt[FDR == 0, FDR := 1 / n.fdr.sets]

   # compute q values from FDR values
   FDR.vect <- q.dt$FDR
   M <- min(FDR.vect)
   I <- max(which(FDR.vect == M))
   K <- 1
   q.vect <- c(rep(M, I), rep(NA, length(FDR.vect) - I))
   while(I < length(q.vect)) { # reset FDRs for larger p values with smaller FDR
      K <- I + 1
      M <- min(FDR.vect[K:length(FDR.vect)])
      I <- max(which(FDR.vect == M))
      q.vect[K:I] <- M
   }

   q.dt[, q := q.vect]
   SSO.pr[, FDR.bin := cut(setScore.pr.p, emp.p.bins)]
   SSO.pr <- merge(SSO.pr,
                   q.dt[, .(FDR.bin, setScore.pr.q = q)],
                   by = "FDR.bin")
   SSO.pr[, FDR.bin := NULL]

   # recover remaining genes in each observed gene set
   soc <- set.obj[, .(A = setID, X = objID)]
   data.table::setkey(soc, A)
   set.obj.pr <- foreach::foreach(i = 1:max(SSO.pr$Rank, na.rm = TRUE),
                                  .combine = rbind) %do% {
      sp.i <- SSO.pr[Rank == i]
      out <- soc[.(sp.i), on = .(A = setID)][, .(setID = A, objID = X)]
      soc <- soc[!.(out), on = .(X = objID)]
      cbind(Rank = i, out)
   }

   # merge set.info with pruned statistics
   set.info <- merge(set.info, SSO.pr, by = "setID", all.x = TRUE)

   ##=========================================================================##
   ##PART 5: Create output
   ##=========================================================================##

   # repack data
   plr.data <- list()
   read.data <- list(n.genes = n.genes, n.sets = n.sets, n.chr = n.chr,
                     n.cov = n.cov, n.set.genes = n.set.genes,
                     file.paths = file.paths, coord = coord,
                     get.objStat = get.objStat, cov.names = cov.names,
                     cov.info = cov.info, no.share = no.share,
                     pos.info = pos.info)
   permute.data <- list(ecdf.std = ecdf.std, gpd.std = gpd.std,
                        no.deconf = no.deconf)
   plr.data$read.data <- read.data
   plr.data$permute.data <- permute.data
   if (rescaled) {
      rescale.data <- list(ac = ac, ecdf.rs = ecdf.rs, gpd.rs = gpd.rs,
                           user.ac = user.ac)
      plr.data$rescale.data <- rescale.data
   }

   # retain arguments and summary statistics
   args0 <- names(formals(get("plR_prune", envir = asNamespace("polylinkR"))))
   prune.args <-  mget(setdiff(args0, "plR.input"), envir = environment()) # collect arguments
   plr.args$prune.args <- prune.args

   prune.summary <- list(set.obj.pr = set.obj.pr, q.dt = q.dt, FDR.pr = FDR.pr)
   plr.summary$prune.summary <- prune.summary

   # capture R session information and record run time
   r0 <- .get_time(start.time = startT)
   plr.session$prune.session <- list(session = sessionInfo(), run.time = r0[[2]])

   # set s3 class and create attributes
   .file_reset(OI = obj.info, SI = set.info, SO = set.obj, pos.info = pos.info) # recreate original file formats
   OUT <- list(set.info = set.info, obj.info = obj.info, set.obj = set.obj)
   plr.out <- .new_plR(BASE = OUT, plR.data = plr.data, plR.args = plr.args,
                       plR.summary = plr.summary, plR.seed = plr.seed,
                       plR.session = plr.session)

   .vrb("\n")
   ftr <- cli::col_cyan(paste0("Finished plR_prune -- run time: ", r0[[1]]))
   pdg <- (80 - nchar(ftr)) / 2
   .vrb(cli::boxx(ftr, padding = c(0, floor(pdg), 0, ceiling(pdg)),
                  border_style = "double", col = "cyan", border_col = "cyan"))
   .vrb("\n")
   return(invisible(plr.out))
}

