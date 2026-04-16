#' @title Gene set enrichment on scores rescaled for genetic autocorrelation
#'
#' @description
#' Rescales gene set scores to account for dependencies within gene sets based
#' on autocorrelated genes and re-evaluates gene set enrichment on the rescaled
#' scores. Note that the \code{obj.info} file must contain appropriately
#' labelled genomic coordinates for autocorrelation to be estimated and gene set
#' rescaling to be performed, otherwise the function will exit with an error
#' message.
#'
#' @param plr_input
#'   \code{plR} class object; output from \code{polylinkR::permute_polylinkr_data}.
#'   Required.
#'
#' @param rescale
#'   \code{logical}; should gene set score (i.e.,\code{setScore.std}) rescaling
#'   be performed? Defaults to \code{TRUE}. If \code{FALSE}, only inter-gene
#'   autocorrelation is estimated without revising enrichment testing for
#'   rescaled (i.e., decorrelated) gene set scores. This is useful for exploring
#'   how parameter settings impact gene set autocorrelation. The rescaling step
#'   can be performed later by passing the output to \code{rescale_polylinkr_data}, or by
#'   providing the \code{ac} object to the \code{ac} argument.
#'
#' @param fast
#'   \code{logical}; should the gene set score rescaling be performed using the
#'   fast mode? Defaults to \code{TRUE}, whereby the analytical scaling factor
#'   \eqn{1 / \sqrt{1 + (m - 1)\hat{\rho}}} is calculated from the covariance
#'   matrix \code{ac} and used to adjust the empirical null and GPD (scale and
#'   threshold) parameters. If \code{FALSE}, the full empirical rescaling is
#'   performed instead, by recalculating the score for every permuted gene set.
#'   This is much more computationally intensive and is intended for users
#'   interested in comparing analytical and empirical results.
#'
#' @param ac
#'   \code{data.frame} or \code{data.table}; user-provided object containing
#'   inter-gene autocorrelations. Include all unique pairs of genes with
#'   non-zero autocorrelation. Genes must be identified by \code{objID}s in
#'   columns \code{objID.A} and \code{objID.B}. The corresponding
#'   autocorrelation value must be in \code{gamma}. Defaults to \code{NULL}.
#'
#' @param cgm_range
#'   \code{numeric}; maximum inter-gene lag used to evaluate autocovariance.
#'   Defaults to \code{2e6} bp or \code{2} cM, depending on genetic distance
#'   measure. Must be within interval \code{[1e5, 5e7]} bp or \code{[0.1, 50]}
#'   cM.
#'
#' @param cgm_bin
#'   \code{numeric}; mimimum number of gene pairs required for bins of empirical
#'   covariances. Default value = \code{30}. Must be in range \code{[10, 1e3]}.
#'   Initially an exponential grid of bin sizes will be generated, favouring
#'   smaller bins at short distances and larger bins for more distant gene pairs.
#'   Bins with too few genes will be successively merged with the larger
#'   bins until the minimum number of genes are reached. This condition is
#'   evaluated across all chromosomes, ensuring no bin is smaller than the
#'   minimum value (other than the final bins).
#'
#' @param cgm_weight_max
#'   \code{numeric}; maximum probability weight for a single lag window.
#'   Defaults to \code{0.1}. Set to \code{1} if no upper bound is desired.
#'   Note that the lower bound is reset to \code{2 / no. fitted lags} if the
#'   chosen value falls below this limit.
#'
#' @param empirical_bayes
#'   \code{numeric}; A character string indicating the empirical Bayes rescaling
#'   framework employed for chromosome-level covariance beta coefficients. Users
#'   may choose to fit either the \code{full} or \code{reduced} model, with the
#'   the former model also allowing for the two coefficient to be correlated.
#'   The default setting (\code{auto}) runs the \code{full} model if more than
#'   \code{15} chromosomes are present, otherwise the \code{reduced} model is
#'   fit. In cases when the \code{full} model is assessed, the  \code{reduced}
#'   model is also fitted and a likelihood-based test is used to decide whether
#'   a correlated random-effects structure is supported.
#'
#' @param min_rho
#'   \code{numeric}; minimum estimated correlation between two gene sets;
#'   values below this are set to \code{0}. Defaults to \code{1e-5}. Range
#'   \code{(0, 0.01]}.
#'
#' @param verbose
#'   \code{logical}; should progress reporting be enabled? Defaults to
#'   \code{TRUE}.
#'
#' @param n_cores
#'   \code{integer}; number of cores for parallel processing. Defaults to
#'   \code{1} or \code{maximum cores - 1}. Must be in \code{[1, maximum cores]}.
#'
#' @param future_plan
#'   \code{character}; parallel backend from the \code{future} package.
#'   Defaults to user \code{n_cores} choice or checks cores, choosing
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
#' @importFrom dqrng dqset.seed dqsample.int
#' @importFrom Rfast rowmeans colsums colCumSums
#' @importFrom robustbase huberM
#' @importFrom rlang local_options
#' @importFrom cli boxx col_cyan col_red style_italic
#'
#' @return
#' A \code{plR} S3 object containing:
#' \describe{
#'   \item{\code{obj.info}}{\code{data.table} for each gene (object).}
#'   \item{\code{set.info}}{\code{data.table} for each gene set.}
#'   \item{\code{set.obj}}{\code{data.table} mapping genes to gene sets.}
#' }
#'
#' Rescaled set scores (i.e., corrected for correlations between genes within
#' sets) are recorded in \code{set.info} as \code{setScore.rs}. Results from
#' the revised enrichment tests are shown in \code{setScore.rs.p}.
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
#'   \item{\code{plr.session}}{
#'   R session information and function run time.
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
#' \code{attributes(X)$plr_data$permute.data}, where \code{X} is the object
#' name. Similarly, the arguments used in \code{plR_read} are in
#' \code{attributes(X)$plr_args$read.args}.
#'
#' The primary data structure of the \code{plR} object can be accessed using
#' \code{print()} or by simply typing the object's name.
#'
#' @examples
#' \dontrun{
#' # Assuming `my_plr` is the result of `polylinkR::permute_polylinkr_data`
#'
#' # Example 1: Basic usage
#' new_plr <- rescale_polylinkr_data(plr_input = my_plr)
#'
#' # Example 2: Estimate autocorrelation only (no rescaling)
#' new_plr <- rescale_polylinkr_data(
#'    plr_input = my_plr,
#'    rescale   = FALSE
#' )
#'
#' # Example 3: Custom variogram model and cores
#' new_plr <- rescale_polylinkr_data(
#'   plr_input = my_plr,
#'   n_cores   = 4
#' )
#'
#' # Example 4: Using a user-provided autocovariance object
#' # Assuming `my_ac` is a valid data.table or data.frame
#' new_plr <- rescale_polylinkr_data(
#'    plr_input = my_plr,
#'    ac        = my_ac
#' )
#'
#' # Or using the `new_plr` object generated by Example 2
#' new_plr <- rescale_polylinkr_data(plr_input = my_rescaled_plr)
#' }
rescale_polylinkr_data <- function(plr_input, rescale = TRUE, fast = TRUE, ac = NULL,
                                    cgm_bin = 30, cgm_range = "auto", cgm_weight_max = 0.05,
                                    empirical_bayes = "auto", min_rho = 1e-5, verbose = TRUE,
                                    n_cores = "auto", future_plan = "auto") {

   ##==========================================================================##
   ## PART 1: Clean data and run checks
   ##==========================================================================##

   # Map snake_case parameters to legacy names for internal use
   plR.input <- plr_input
   cgm.bin <- cgm_bin
   cgm.range <- cgm_range
   cgm.wt.max <- cgm_weight_max
   emp.bayes <- empirical_bayes
   min.rho <- min_rho
   n.cores <- n_cores
   fut.plan <- future_plan

   # track function run time
   startT <- Sys.time()

   # control message reporting
   rlang::local_options(verbose = verbose, .frame = environment())

   # perform checks and unpack required plR objects
   .check_plr_object(f = "rescale", ENV = environment())
   .check_arguments(f = "rescale", ENV = environment())
   .unpack_plr_object(plr = plR.input, ENV = environment())
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

   #checks passed, announce function
   if (rescale) {
      if (is.null(ac)) {
         hdr <- c("estimate intergene autocovariance &",
                  "test enrichment on decorrelated gene set scores")
      } else {
         hdr <- "test enrichment on decorrelated gene set scores"
      }
   } else {
      hdr <- "estimate intergene autocovariance"
   }

   hdr[1] <- paste("Running rescale_polylinkr_data:", hdr[1])
   pdg <- (80 - max(nchar(hdr))) / 2
   .verbose_msg(cli::boxx(hdr, padding = c(0, floor(pdg), 0, ceiling(pdg)),
                  align = "center", col = "cyan", border_col = "cyan"))

   # report collated warnings and messages
   .verbose_msg("\n")
   if (!is.null(WMSS)) warning(WMSS, immediate. = T, call. = F)
   if (!is.null(MSS)) .verbose_msg(paste0("\n", MSS, "\n"))
   if (!is.null(pWMSS)) warning(pWMSS, immediate. = T, call. = F)
   if (!is.null(pMSS)) .verbose_msg(paste0("\n", pMSS, "\n"))
   .verbose_msg("\n")

   ##==========================================================================##
   ## PART 2: Calculate genomic autocovariance
   ##==========================================================================##

   if (is.null(ac)) { # infer genetic autocovariance
      .verbose_msg(cli::style_italic("Estimating genetic autocovariance:\n"))

      # generate exponential lag bin sizes
      if (cgm.range == "auto") { # determine maximum correlogram range if not provided
         maxL <- ifelse(coord == "cM", 3L, 3e6L)
      }

      .verbose_msg(paste0("Creating correlogram lag bins [max ", maxL, coord,
                  " lag and >= ", cgm.bin, " gene pairs per bin]\n"))

      minL <- if (coord == "cM") { # define smallest possible lag bin (capturing 0-lags)
         minL <- min(obj.info[, diff(midpos), by = chr][V1 > 0]$V1)
         minL <- 10 ^ floor(log10(minL))
      } else {
         minL <- 1
      }

      # create input data
      oi0 <- obj.info[, .(objID, chr, objStat, midpos, startpos, endpos)]
      if (!plr.args$read.args$rem.genes) { # exclude genes with identical positions
         oi0 <- data.table::copy(obj.info)
         oi0[, objDup := paste(chr, startpos.base, endpos.base, sep = ".")] # use base coordinates
         oi0 <- oi0[!duplicated(objDup)]
         oi0[, objDup := NULL]
         objDup <- setdiff(1:n.genes, oi0$objID)
      }

      oi0 <- oi0[order(startpos, endpos)]
      oi0 <- split(oi0, oi0$chr)
      cBins <- .get_cgm_bins(maxL = maxL, minL = minL,
                             cgm.bin = cgm.bin, OI = oi0)
      nBins <- length(cBins) - 1

      # calculate empirical covariance and point estimates for each lag bin
      cv0 <- foreach::foreach(x = oi0) %do% {
         xD <- dist(x$midpos)
         wD <- which(xD <= maxL) # identify values within maximum lag range
         n <- x[, .N]

         # Convert linear indices back to pair coordinates (i, j)
         j <- floor(n + 0.5 - sqrt((n - 0.5) ^ 2 - 2 * (wD - 1))) # column index
         i <- wD - (j - 1) * (n - j / 2) + j # row index

         # extract and bin covariances
         xIJ <- cbind(x[j, .(osJ = objStat, sJ = startpos, eJ = endpos)],
                      x[i, .(osI = objStat, sI = startpos, eI = endpos)])
         xIJ[, `:=` (cov.emp = osI * osJ, gI = sI - ifelse(eJ > eI, eI, eJ),
                     gU = sqrt((eI - sI) * (eJ - sJ)))]
         data.table::set(xIJ, j = c("sI", "sJ", "eI", "eJ", "osI", "osJ"),
                         value = NULL)
         xIJ[, cov.ovlp := ifelse(gI <= 0, abs(gI) + 1, 0) / (gU + 1)] # calculate overlaps
         xIJ[, pos.lag := xD[wD]] # add lags
         xIJ[, bin := cut(pos.lag, cBins, labels = 1:nBins, right = FALSE)] # bin lags

         out <- xIJ[, .(.N, pos.lag = mean(pos.lag), cov.ovlp = mean(cov.ovlp),
                        cov.mean = mean(cov.emp), cov.med = median(cov.emp)),
                    keyby = bin] # binned means and medians
         hm <- xIJ[, robustbase::huberM(cov.emp, se = TRUE), keyby = bin] # huber M estimator
         out[hm, `:=` (cov.hubM.mu = mu, cov.hubM.se = SE)]
         if (n.chr > 1) { # keep covariances for global estimation
            td0 <- split(xIJ$cov.emp, f = xIJ$bin)
         } else {
            td0 <- NULL
         }
         list(cbind(chr = x[1]$chr, out) , td0)
      }

      # extract objects
      cN <- do.call(rbind, lapply(cv0, "[[", 1)) # binned estimates

      .verbose_msg("Assigning emprical covariances to lag bins\n")
      if (n.chr > 1) { # add complete genome
         cv0 <- lapply(cv0, "[[", 2)
         gv0 <- lapply(1:nBins, function(b){
            cvB <- unlist(lapply(cv0, "[[", b)) # extract bins
            c(median(cvB), unlist(robustbase::huberM(cvB, se = TRUE)[c(1, 4)]))
         })
         gv0 <- do.call(rbind, gv0) # aggregate results

         gN <- cN[, .(pos.lag = sum(pos.lag * N, na.rm = TRUE),
                      cov.mean = sum(cov.mean * N, na.rm = TRUE),
                      cov.ovlp = sum(cov.ovlp * N, na.rm = TRUE),
                      N = sum(N, na.rm = TRUE)), keyby = bin]
         gN[, `:=` (pos.lag = pos.lag / N, cov.mean = cov.mean / N,
                    cov.ovlp = cov.ovlp / N)]
         gN[, `:=` (cov.med = gv0[, 1], cov.hubM.mu = gv0[, 2],
                    cov.hubM.se = gv0[, 3])]

         cN <- rbind(cN, cbind(chr = "All", gN)) # include genome wide aggregated data
      }
      rm(cv0); gc(verbose = FALSE)

      # add bin boundaries
      cN[, bin.min := cBins[bin]]
      cN[, bin.max := cBins[-1][bin]]
      cN[, bin.size := bin.max - bin.min]

      .verbose_msg(paste("Fitting autocovariance function\n"))

      # determine lag bin weights for covariance function fitting
      cgm.wt.min <- 2 / (nBins - 1) # check sufficient minimum weight (ignore lag = 0 bin)
      if (cgm.wt.max < cgm.wt.min) cgm.wt.max <- cgm.wt.min

      cN[, lag.wts := (N / bin.size), by = chr] # gene density
      cN[, lag.wts := .cap_probs(matrix(lag.wts, nrow = 1),
                                        maxP = cgm.wt.max), by = chr] # avoid too much emphasis on smallest category
      cN[, bin := match(bin, bin) - 1] # set 0-lag bin index to 0

      # fit model to empirical covariances
      if (n.chr > 1) {
         cN0 <- cN[chr == "All"]
      } else {
         cN0 <- cN
      }

      # paramter bounds
      pInit <- c(1, 1, maxL * 1/2)
      pUpper <- c(1e3, 1e3, maxL * 1e3)
      pLower <- c(1e-6, 1e-6,  maxL * 1e-6)

      gFit <- with(cN0,
                   optim(par = pInit, fn = .cgm_fun_fit, h = pos.lag,
                         w = lag.wts, cov_emp = cov.hubM.mu,
                         se = cov.hubM.se, ovlp = cov.ovlp, method = "L-BFGS-B",
                         lower = pLower, upper = pUpper))
      gPar <- gFit$par # extract global parameters

      if (n.chr > 1) { # chromosomal-leval parameter fitting required
         # Estimate chromosome parameters
         fopt <- list(packages = c("foreach", "data.table"),
                      globals = c("cN", "gPar", "pInit", "pLower", "pUpper",
                                  ".cov_fun", ".cgm_fun_fit"))

         cFit <- foreach::foreach(vb.i = split(cN, cN$chr)[1:n.chr],
                                  .options.future = fopt) %dofuture% {
            with(vb.i,
                 optim(par = pInit[1:2], fn = .cgm_fun_fit, fixed = gPar[3],
                       h = pos.lag, w = lag.wts, cov_emp = cov.hubM.mu,
                       se = cov.hubM.se, ovlp = cov.ovlp,  method = "L-BFGS-B",
                       lower = pLower[1:2], upper = pUpper[1:2]))
         }
         cPar <- sapply(cFit, function(x) x$par)

         # combine with global fitting results
         gFit <- c(list(gFit), cFit)
         names(gFit) <- c("All", 1:n.chr)

         # regularise covariance function parameter estimates for each chromosome
         .verbose_msg("Apply empirical Bayes to revise chromosome parameters\n")
         ebPar0 <- .empirical_bayes_estimate(cPar = cPar, gPar = gPar, cN0 = cN,
                                  corr = FALSE)
         if (emp.bayes == "auto") {
            emp.bayes <- ifelse(n.chr >= 15, "full", "reduced")
         }

         if (emp.bayes == "full") {
            .verbose_msg("Comparing full and marginal models\n")
            ebPar <- .empirical_bayes_estimate(cPar = cPar, gPar = gPar, cN0 = cN,
                                    corr = TRUE)
            LR <- 2 * (ebPar0$logLik - ebPar$logLik)
            pval <- 0.5 * pchisq(LR, df = 1, lower.tail = FALSE)
            LR.test <- list(REML = LR, p = pval)
            if (pval < 0.05) {
               cPar.reg <- do.call(cbind, ebPar$beta.chr)
            } else {
               cPar.reg <- do.call(cbind, ebPar0$beta.chr)
            }
            ebPar0 <- list(ebPar, ebPar0)
            names(ebPar0) <- c("full", "reduced")
         } else {
            cPar.reg <- do.call(cbind, ebPar0$beta.chr)
            LR.test <- NULL
         }

         cPar.reg <- pmax(cPar.reg, 0) # set negative beta coefficients to 0

         # create summary only
         cgm.summary <- data.frame(1:n.chr, rep(gPar[3], n.chr),
                                   t(cPar), cPar.reg)
         data.table::setDT(cgm.summary)
         data.table::setnames(cgm.summary, c("chr", "scale", "beta.LD",
                                             "beta.ovlp", "beta.LD.reg",
                                             "beta.ovlp.reg"))
         attr(cgm.summary, "emp.Bayes") <- ebPar0
         attr(cgm.summary, "LR.test") <- LR.test
      } else {  # create summary only
         cgm.summary <- data.table::data.table(chr = obj.info$chr[1], t(gPar))
         data.table::setnames(cgm.summary, c("chr", "beta.LD", "beta.ovlp",
                                             "scale.LD", "scale.ovlp"))
      }

      ckFit <- sapply(gFit, function(gf) gf$convergence == 0) # check convergence
      if (!all(ckFit)) {
         gErr <- c("Global", paste("chr", unique(obj.info$chr.orig)))[!ckFit]
         warning("Correlogram function fitting did not converge for",
                 paste(gErr, collpase = ", "), "\nRecommend checking cgm.fit ",
                 "attribute and rerunning plR_rescale with altered parameters",
                 immediate. = TRUE)
      }

      # determine empirical covariances
      .verbose_msg("Generating intergene covariances\n")

      ac <- foreach::foreach(x = oi0, par = split(cgm.summary, 1:n.chr)) %do% {
         xD <- dist(x$midpos)
         wD <- which(xD <= maxL) # identify values within maximum lag range
         n <- x[, .N]

         # Convert linear indices back to pair coordinates (i, j)
         j <- floor(n + 0.5 - sqrt((n - 0.5) ^ 2 - 2 * (wD - 1))) # column index
         i <- wD - (j - 1) * (n - j / 2) + j # row index

         xIJ <- cbind(x[j, .(objID.A = objID, sJ = startpos, eJ = endpos)],
                      x[i, .(objID.B = objID, sI = startpos, eI = endpos)])
         xIJ[, `:=` (gI = sI - ifelse(eJ > eI, eI, eJ),
                     gU = sqrt((eI - sI) * (eJ - sJ)))]
         data.table::set(xIJ, j = c("sI", "sJ", "eI", "eJ"), value = NULL) # remove columns
         xIJ[, cov.ovlp := ifelse(gI <= 0, abs(gI) + 1, 0) / (gU + 1)] # define overlap
         xIJ <- cbind(xIJ[, .(objID.A, objID.B, cov.ovlp)], pos.lag = xD[wD]) # add covariances

         # calculate pairwise covariances (= correlations)
         xIJ[, rho := .covariance_function(h = pos.lag, ovlp = cov.ovlp, scale = par$scale,
                               beta1 = par$beta.LD.reg,
                               beta2 = par$beta.ovlp.reg)]
         cbind(chr = x$chr[1],  xIJ[rho >= min.rho])
      }
      ac <- do.call(rbind, ac)

      user.ac <- FALSE
   } else { # user has provided valid ac object
      .verbose_msg(paste("Valid ac object detected:",
                 "Skipping autocovariance calculation step\n"))
      ac <- ac[objID.A != objID.B][rho >= min.rho] # remove unused values
      if (exists("plr.summary$rescale.summary")) { # ac previously estimated using plR, retain args and summaries
         list2env(plr.args$rescale.args, envir = environment())
         list2env(plr.summary$rescale.summary, envir = environment())
         plr.summary$rescale.summary <- NULL
         plr.args$rescale.args <- NULL
      } else { # ac not previously estimated using plR
         gFit <- NULL
         cgm.summary <- NULL
      }
      user.ac <- TRUE
   }

   ##==========================================================================##
   ## PART 3: Decorrelate set scores
   ##==========================================================================##

   if (rescale) {
      .verbose_msg(cli::style_italic("\nRescaling null set score distributions:\n"))

      # get common variables for analytical and permutation-based corrections
      n.perm <- plr.args$permute.args$n.perm
      gpd.cutoff <- plr.args$permute.args$gpd.cutoff
      n.max <- max(set.info$setN)
      alt <- plr.args$permute.args$alt

      # create quantiles to sample: use logistic quantiles to ensure higher precision at margins
      q.bnd <- 10 ^ -(floor(log10(n.perm)) - 1) # constrain quantile bounds to avoid poor estimation
      eQnt <- plogis(seq(qlogis(q.bnd), qlogis(1 - q.bnd), length.out = 500))
      eQnt <- sort(unique(c(1 - gpd.cutoff, eQnt))) # add GPD estimation anchor

      rho.hat <- sum(ac$rho) * 2  / (n.genes * (n.genes - 1)) # mean autocorrelation
      rI <- sqrt(1 + (1:n.max - 1) * rho.hat) # analytical rescaling factor

      if (fast) {
         .verbose_msg(paste("fast = TRUE: Rescaling empirical CDFs and GPD scale",
                    "coefficient using analytical scaling factor\n"))

         if (n.sets == 1) {
            qX <- get(x = "x", envir = environment(ecdf.std)) / rI[n.max]
            ecdf0 <- stats::approxfun(x = qX, y = eQnt, method = "linear",
                                      rule = 2, ties = "ordered")
            gpd0 <- data.table::copy(gpd.std)
            gpd0[, scale := scale / rI[n.max]]
            rSumm <- list(setN = n.max, rs.analyt = 1 / rI[n.max])
            data.table::setDT(rSumm)
         } else {
            # rescale ECDF
            ecdf0 <- foreach::foreach(eI = ecdf.std[2:n.max], m = rI[-1]) %do% {
               qX <- get(x = "x", envir = environment(eI)) / m
               stats::approxfun(x = qX, y = eQnt, method = "linear", rule = 2,
                                ties = "ordered")
            }
            ecdf0 <- c(NA, ecdf0)

            # revise GDP tail distribution
            gpd0 <- data.table::copy(gpd.std)
            gpd0[, scale := c(NA, scale[-1] / rI[-1])]
            rSumm <- list(setN = 1:n.max, rs.analyt = 1 / rI)
            data.table::setDT(rSumm)
         }

         gSumm <- NULL
         eSumm <- NULL
      } else {
         n.boot <- plr.args$permute.args$n.boot

         .verbose_msg(paste0("fast = FALSE: Set-wise rescaling of empirical CDFs and",
                     " refitting GPDs for ", n.boot, " x ",
                     n.perm / ifelse(n.perm / 1e6 >= 1, 1e6, 1e3),
                     ifelse(n.perm / 1e6 >= 1, "M", "k"), " gene sets\n"))

         if (n.sets == 1) {
            th0 <- n.max
         } else { # determine the sub-sampled gene set sizes (interpolate the remainder)
            wM <- which.max(cumsum(1:200) * 10 >= n.max)
            th0 <- cumsum(unlist(lapply(1:wM, rep, 10)))
            th0 <- c(th0[2:sum(th0 < n.max)], n.max) # exclude unused values and cap at max set size
         }

         fpc <- sqrt((n.genes - th0) / (n.genes - 1)) # finite population correction
         n.th <- length(th0) # no. set sizes estimated
         n.tail <- ceiling(gpd.cutoff * n.perm) # gpd tail length

         # determine optimal block size for permutations
         minB <- max(n.tail, min(n.perm / n.cores, 5e4)) / 1e4 # minimum block size
         possB <- which(n.perm %% (1:5 * 1e4) == 0) # possible block sizes
         n.block <- possB[which.min(abs(minB - possB))] * 1e4 # final block size
         x.th <- 1 - n.tail / n.block # tail threshold in inner permutation

         eps <- sqrt(eQnt * (1 - eQnt) / n.perm) # interval around quantiles to estimate measurement error
         os0 <- obj.info$objStat.std # gene scores
         if (alt == "lower") { # take extremes from lower tail
            os0 <- -os0
         }
         inI <- list(g0 = matrix(-1e6, nrow = n.th, ncol = n.tail),
                     e0 = replicate(n.th, list()),
                     r0 = 0) # inner loop initial combine object
         inO <- list(g0 = NULL,
                     e0 = list(qX = replicate(n.th, NULL),
                               qE = replicate(n.th, NULL)),
                     r0 = 0) # outher loop initial combine object
         dqi <- replicate(n.boot,
                          dqrng::generateSeedVectors(n.perm / n.block),
                          simplify = FALSE) # reproducible seeds

         rs0 <- matrix(0, nrow = n.max, ncol = n.block) # create rescaling matrix
         SS0 <- data.table::data.table(sID = rep(1:n.block, each = n.max),
                                       sN = rep(1:n.max, n.block), key = "sN") # all possible pairs to this point
         ac0 <- ac[, .(A = objID.A, B = objID.B, CV = 2 * rho)] # ensure 2 * covariance
         data.table::setkey(ac0, A, B)

          fopt <- list(packages = c("data.table", "foreach", "dqrng", "Rfast"),
                      globals = c("n.genes", "n.max", "n.tail", "n.sets", "n.th",
                                  "n.block", "x.th", "fpc", "os0", "th0", "inI",
                                  "dqI", "rs0", "ac0", "SS0", ".cmb", ".get_cov0"),
                      seed = seed)

         progressr::with_progress({
            prog <- progressr::progressor(steps = n.boot)
            std.null <- foreach::foreach(I = dqi, .combine = .cmb, .init = inO) %do% { # iterate over each bootstrap replicate
               suppressPackageStartupMessages(
                  out <- foreach::foreach(x = I, .combine = .cmb, .init = inI,
                                          .options.future = fopt) %dofuture% {
                     dqrng::dqset.seed(x)
                     csX <- replicate(n.block,
                                      dqrng::dqsample.int(n.genes, n.max)) # faster than dqsample::dqsample
                      rsX <- .get_covariance_null(csj = csX, rsX = rs0, set_scores = SS0, autocov_data = ac0) # get covariances
                     csX <- matrix(os0[csX], nrow = n.max, ncol = n.block) # get gene scores

                     if (n.sets == 1) { # only take final entry
                        rsX <- Rfast::colsums(rsX)
                        rs.mean <- sqrt(n.max / (n.max + rsX))
                        csX <- Rfast::colsums(csX)
                        csX <- matrix(csX / (sqrt(n.max + rsX) * fpc), nrow = 1) # standardise set scores
                     } else { # cumulative sum
                        rsX <- Rfast::colCumSums(rsX)
                        rsD <- sqrt(1:n.max + rsX)
                        rs.mean <- sqrt(1:n.max) / Rfast::rowmeans(rsD)
                        csX <- Rfast::colCumSums(csX)
                        csX <- csX[th0, ] / (rsD[th0, ] * fpc) # standardise set scores
                     }

                     # extract tail values
                     if (n.block == n.tail) {
                        g0 <- csX
                        csX <- split(csX, 1:n.th)
                     } else {
                        csX <- split(csX, 1:n.th)
                        g0 <- t(sapply(csX, FUN = function(csI) {
                           x.th <- quantile(csI, x.th)
                           x.id <- which(csI >= x.th)
                           if (length(x.id) > n.tail) {
                              x.out <- which(csI == x.th)[1:(length(x.id) - n.tail)] # identify values to remove
                              x.id <- x.id[!(x.id %in% x.out)]
                           }
                           csI[x.id]
                        }))
                     }

                      e0 <- lapply(csX, FUN = function(csI) { # collect raw values
                         as.list(csI)
                      })

                     list(g0 = g0, e0 = e0, r0 = rs.mean)
                  }
               )
               # estimate ecdf and gpd fit for each null
               out$g0 <- t(apply(out$g0, 1, .fit_gpd))
               e0 <- lapply(out$e0, .get_quant, eQnt = eQnt, eps = eps)
               out$e0 <- list(qX = lapply(e0, "[[", 1), qE = lapply(e0, "[[", 2))
               prog()
               out
            }
         })

         # extract results
         g0 <- std.null$g0
         e0 <- std.null$e0$qX
         eE <- std.null$e0$qE
         r0 <- std.null$r0 * n.block / (n.boot * n.perm)
         rm(std.null); gc(verbose = FALSE)

         # collate rescaling statistics
         if (n.sets == 1) {
            rSumm <- list(setN = n.max, rs.perm = r0, rs.analyt = 1 / rI[n.max])
         } else {
            rSumm <- list(setN = 1:n.max, rs.perm = r0, rs.analyt = 1 / rI)
         }

         # check for inviable fits in gpd estimation
         gpd.miss <- sum(!complete.cases(g0)) / nrow(g0)
         if (gpd.miss > 0.5) {
            warning("GPD fitting failed for >50% of nulls -- halting function:",
                    "\ntry again with different gpd fitting options",
                    "\n[or set fast = TRUE for analytical rescaling]\n",
                    immediate. = T, call. = F)
            permute <- FALSE
            gpd0 <- NULL
            gSumm <- NULL
            ecdf0 <- NULL
            eSumm <- NULL
         } else { # smooth estimates and interpolate missing set sizes
            .verbose_msg(paste("Smoothing empirical CDF and GPD parameter estimates",
                       "and interpolating missing values\n"))

            K <- .est_ss_cov(x = th0, n.genes = n.genes) # covariance matrix
            sm.k <- min(max(10, round(n.th / 4)), 40) # deterministic basis for gam smoother
            sm0 <- mgcv::smoothCon(mgcv::s(x, k = sm.k, bs = "ps"),
                                   data = data.frame(x = log(th0)))[[1]] # get inner smoothing function

            # empirical CDFs
            eE <- lapply(eE, FUN = function(x) { # get quantile measure error
               eQnt * rev(eQnt) * (x / (2 * eps)) ^ 2 / n.perm
            })

            if (alt == "lower") { # ensure correct tail orientation
               e0 <- lapply(e0, function(x) -rev(x))
               eE <- lapply(eE, function(x) rev(x))
            }

            # calculate means
            e0mean <- sapply(e0, Rfast::colmeans)
            eEmean <- sapply(eE, Rfast::colmeans)

            # fit smoothing functions
            fopt <- list(globals = c("th0", "K", "n.th", "n.max", "sm0",
                                     ".par_smooth", "Predict.matrix.wls.smooth",
                                     "smooth.construct.wls.smooth.spec"),
                         packages = "mgcv", seed = seed)

            ecdf0 <- foreach::foreach(ei = split(e0mean, eQnt),
                                      ej = split(eEmean, eQnt),
                                      .options.future = fopt,
                                      .combine = cbind) %dofuture% {
               eSmFit <- .parallel_smooth(y = ei, x = th0, K = K, n.th = n.th,
                                     sm0 = sm0, tau = ej)
               predict(eSmFit, newdata = data.frame(x = log(2:n.max), oneW = 1))
            }
            exc.z <- ecdf0[,  which(eQnt == 1 - gpd.cutoff)] # get smoothed GDP exceedence threshold
            cut.z <- ecdf0[,  length(eQnt)] # get smoothed GDP exceedence threshold
            ecdf0 <- c(NA, lapply(split(ecdf0, 2:n.max), function (qX) { # create ecdf list
               stats::approxfun(x = qX, y = eQnt, method = "linear", rule = 2,
                                ties = "ordered")
            }))

            # compile summary statistics
            eSumm <- cbind(e0mean, sapply(e0, Rfast::colVars, std = TRUE))
            eSumm <- data.table::data.table(Stat = rep(c("Q.mean", "Q.SD"),
                                                       each = n.th),
                                            setN = rep(th0, 2), t(eSumm))
            data.table::setnames(eSumm, c("Stat", "setN", paste0("q", eQnt)))
            rm(e0, eE, e0mean, eEmean); gc(verbose = FALSE) # clean up

            # GPD fits
            g0 <- data.frame(rep(1:n.boot, each = n.th), rep(th0, n.boot), g0)
            data.table::setDT(g0)
            data.table::setnames(g0, c("Boot", "setN", "scale", "shape",
                                       "scale.var", "shape.var"))
            gSumm <- data.table::data.table(setN = th0)
            data.table::setkey(gSumm, setN)
            gpd0 <- data.table::data.table(setN = 1:n.max) # collect smoothed values
            for (p in c("scale", "shape")) {
               gM <- g0[, mean(get(p)), by = setN]
               gV <- g0[, mean(get(paste0(p, ".var"))), by = setN]
               gSmFit <- .parallel_smooth(y = gM$V1, x = th0, K = K, n.th = n.th,
                                     sm0 = sm0, tau = gV$V1)
               yFit <- predict(gSmFit, newdata = data.frame(x = log(2:n.max),
                                                            oneW = 1))
               gpd0[, (p) := c(NA, yFit)] # add smoothed gpd parameter
               # collate summary statistics
               gSumm[g0[, sum(!is.na(scale)), by = setN], bootN := V1] # non-missing bootstraps
               gSumm[gM, (paste0(p, ".mean")) := V1] # parameter mean
               gSumm[g0[, sd(get(p)), by = setN], (paste0(p, ".sd")) := V1] # parameter std dev
            }
            gpd0[, exc.z := c(NA, exc.z)] # add GDP exceedence threshold
            gpd0[, cut.z := c(NA, cut.z)] # add 0.999 quantile threshold

            # include associated cutoff p-value gpd minimum non-0 p values and quantiles
            .verbose_msg(paste("Identifying and setting minimum possible p values\n"))
            fopt <- list(packages = "data.table",
                         globals = c("gpd.cutoff", "q.bnd", ".get_gpd_mins"),
                         seed = seed)

            gpdM <- foreach::foreach(gpdI = split(gpd0[-1], 2:n.max),
                                     .options.future = fopt,
                                     .combine = rbind) %dofuture% {
            .get_gpd_minimums(gpdI, gpd.cutoff = gpd.cutoff, q.bnd = q.bnd)
                                     }
            gpd0[, adj.p := c(NA, gpdM[, "adj.p"])]
            gpd0[, min.z := c(NA, gpdM[, "min.z"])]
            gpd0[, min.p := c(NA, gpdM[, "min.p"])]
            rm(g0); gc(verbose = FALSE) # clean up

            # add gpd and ECDF summaries to rescaling summary
            rSumm$gpd.scale.rs <- gpd0$scale / gpd.std$scale
            rSumm$gpd.shape.rs <- gpd0$shape / gpd.std$shape

            if (n.sets == 1) {
               qS <- get(x = "x", envir = environment(ecdf.std))
               qR <- get(x = "x", envir = environment(ecdf0))
               rSumm$quant.rs.med <- median(qR / qS)
               rSumm$quant.rs.05 <- quantile(qR / qS, 0.05)
               rSumm$quant.rs.95 <- quantile(qR / qS, 0.95)
            } else {
               qRS <- foreach::foreach(eS = ecdf.std[-1], eR = ecdf0[-1],
                                       .combine = cbind) %do% {
                  qS <- get(x = "x", envir = environment(eS))
                  qR <- get(x = "x", envir = environment(eR))
                  qR / qS
               }
               rSumm$quant.rs.med <- c(NA, Rfast::colMedians(qRS))
               rSumm$qunat.rs.05 <- c(NA, apply(qRS, 2, quantile, 0.05))
               rSumm$quant.rs.95 <- c(NA, apply(qRS, 2, quantile, 0.95))
            }
            data.table::setDT(rSumm)
         }
         .verbose_msg("\n")
      }

      # set pre-scaled ECDF and GDP to null
      gpd.std <- NULL
      ecdf.std <- NULL
   } else {
      ecdf0 <- NULL
      eSumm <- NULL
      gpd0 <- NULL
      gSumm <- NULL
      rSumm <- NULL
   }

   ##==========================================================================##
   ## PART 4: Calculate test statistics
   ##==========================================================================##

   if (rescale) {
      .verbose_msg(cli::style_italic("\nCalculating test statistics:\n"))

      # compute rescaling factor for observed set scores
      SSI <- set.obj[, .(oID = objID), keyby = setID][, .(sID = setID, oID)]
      ac0 <- ac[, .(A = objID.A, B = objID.B, CV = 2 * rho)]

      obs.cov <- rep(0, n.sets)
       oc0 <- .get_covariance(SSi = SSI, autocov_data = ac0)
      obs.cov[oc0$sID] <- oc0$CV

      # calculate observed rescaled set scores
      set.info[, setScore.rs := setScore.std * sqrt(setN / (setN + obs.cov))]

      # calculate p values
      .verbose_msg("Calculating p values for each gene set\n")
      if (alt == "lower") {
         .get_p <- .get_p_lt
      } else {
         .get_p <- .get_p_ut
      }
      set.info[, setScore.rs.p := .get_p(ss = setScore.rs, n = setN[1],
                                         g0 = gpd0, e0 = ecdf0), by = setN]
   }

   ##==========================================================================##
   ##PART 5: Create output
   ##==========================================================================##

   # repack data and add new data
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
   rescale.data <- list(ac = ac, ecdf.rs = ecdf0, gpd.rs = gpd0,
                        user.ac = user.ac)
   plr.data$rescale.data <- rescale.data

   # retain arguments and summary statistics
   args0 <- names(formals(get("plR_rescale", envir = asNamespace("polylinkR"))))
   rescale.args <-  mget(setdiff(args0, "plR.input"), envir = environment()) # collect arguments
   plr.args$rescale.args <- rescale.args

   if ("rescale.summary" %in% names(plr.summary)) {
      plr.summary$rescale.summary$gpd.rs.summary <- gSumm
      plr.summary$rescale.summary$ecdf.rs.summary <- eSumm
      plr.summary$rescale.summary$cgm.fit <- gFit
      plr.summary$rescale.summary$cgm.summary <- cgm.summary
   } else {
      rescale.summary <- list(cgm.fit = gFit, cgm.summary = cgm.summary,
                              gpd.rs.summary = gSumm, ecdf.rs.summary = eSumm,
                              rescale.stats = rSumm)
      plr.summary$rescale.summary <- rescale.summary
   }

   # capture R session information and record run time
   r0 <- .get_time(start.time = startT)
   if ("rescale.session" %in% names(plr.session)) {
      rs.sess <- plr.session$rescale.session
      rescale.session <- list()
      rescale.session$session <- list(rs.sess$session, sessionInfo())
      rescale.session$run.time <- c(rs.sess$run.time, r0[[2]])
   } else {
      rescale.session <- list(session = sessionInfo(), run.time = r0[[2]])
   }
   plr.session$rescale.session <- rescale.session

   # set s3 class and create attributes
   .file_reset(OI = obj.info, SI = set.info, SO = set.obj, pos.info = pos.info) # recreate original file formats
   OUT <- list(set.info = set.info, obj.info = obj.info, set.obj = set.obj)
   plr.out <- .new_plr(BASE = OUT, plr_data = plr.data, plr_args = plr.args,
                       plr_summary = plr.summary, plr_seed = plr.seed,
                       plr_session = plr.session)

   .verbose_msg("\n")
   ftr <- cli::col_cyan(paste0("Finished rescale_polylinkr_data -- run time: ", r0[[1]]))
   pdg <- (80 - nchar(ftr)) / 2
   .verbose_msg(cli::boxx(ftr, padding = c(0, floor(pdg), 0, ceiling(pdg)),
                  border_style = "double", col = "cyan", border_col = "cyan"))

   return(invisible(plr.out))
}


#' @title Gene set enrichment on scores rescaled for genetic autocorrelation (deprecated)
#' @description
#' This function is deprecated. Please use \code{rescale_polylinkr_data()} instead.
#' @param plR.input Deprecated. Use \code{plr_input}.
#' @param rescale Deprecated. Use \code{rescale}.
#' @param fast Deprecated. Use \code{fast}.
#' @param ac Deprecated. Use \code{ac}.
#' @param cgm.bin Deprecated. Use \code{cgm_bin}.
#' @param cgm.range Deprecated. Use \code{cgm_range}.
#' @param cgm.wt.max Deprecated. Use \code{cgm_weight_max}.
#' @param emp.bayes Deprecated. Use \code{empirical_bayes}.
#' @param min.rho Deprecated. Use \code{min_rho}.
#' @param verbose Deprecated. Use \code{verbose}.
#' @param n.cores Deprecated. Use \code{n_cores}.
#' @param fut.plan Deprecated. Use \code{future_plan}.
#' @export
plR_rescale <- function(plR.input, rescale = TRUE, fast = TRUE, ac = NULL,
                         cgm.bin = 30, cgm.range = "auto", cgm.wt.max = 0.05,
                         emp.bayes = "auto", min.rho = 1e-5, verbose = TRUE,
                         n.cores = "auto", fut.plan = "auto") {
   .Deprecated("rescale_polylinkr_data", package = "polylinkR",
               msg = "plR_rescale() is deprecated. Use rescale_polylinkr_data() instead.")
   rescale_polylinkr_data(
      plr_input = plR.input,
      rescale = rescale,
      fast = fast,
      ac = ac,
      cgm_bin = cgm.bin,
      cgm_range = cgm.range,
      cgm_weight_max = cgm.wt.max,
      empirical_bayes = emp.bayes,
      min_rho = min.rho,
      verbose = verbose,
      n_cores = n.cores,
      future_plan = fut.plan
   )
}
