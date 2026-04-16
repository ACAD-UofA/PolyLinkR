#' @title Gene set enrichment with control for confounding covariates
#'
#' @description
#' Performs gene set enrichment testing using a gene-wise permutation procedure,
#' while accounting for confounding variables. Note that the deconfounding step
#' requires that appropriately labeled covariate columns are provided in the
#' \code{obj.info} file (otherwise standard gene score residuals are passed to
#' permutations). Applies a Generalized Pareto Distribution (GPD) to the tail of
#' the empirical null for accurate estimation of small p-values.
#'
#' @param plr_input
#'   \code{plR} class object; output from \code{polylinkR::read_polylinkr_data}.
#'   Required.
#'
#' @param permute
#'   \code{logical}; should the permutation step be performed? Defaults to
#'   \code{TRUE}. If \code{FALSE}, only gene score deconfounding is performed
#'   without estimating set scores. If \code{TRUE}, users must provide
#'   confounder covariates in \code{obj.info}. This is useful for exploring
#'   how parameter settings impact deconfounded gene scores; the permutation
#'   step can be performed later by passing the output back into
#'   \code{permute_polylinkr_data}.
#'
#' @param n_permutations
#'   \code{integer}; number of permutations. Defaults to \code{1e5}. Must be
#'   in the range \code{[5e4L, Inf)} and be exactly divisible by \code{1e4}.
#'
#' @param n_bootstraps
#'   \code{integer}; number of bootstrap replicates for null inference.
#'   Defaults to \code{30L}. Must be in the range \code{[5, Inf)}.
#'
#' @param alternative
#'   \code{character}; direction of the hypothesis test. Defaults to
#'   \code{"upper"} for enrichment in the upper tail (large set scores).
#'   Alternatively, \code{"lower"} tests for enrichment in the lower tail
#'   (small values). When \code{"lower"} is chosen, data are internally
#'   negated in functions performing p-value estimation.
#'
#' @param md_method
#'   \code{character}; determines whether raw covariate data or ranks are
#'   used in Mahalanobis distance calculations. Defaults to \code{"robust"},
#'   where Mahalanobis distances are converted to ranks and Spearman's metric
#'   is used to calculate the covariance matrix. The alternate option
#'   \code{"raw"} uses the original covariates and Pearson's covariance for
#'   scaling.
#'
#' @param kernel_boundary
#'   \code{numeric}; flanking region around each gene where weights of
#'   overlapping genes inside the region are set to 0. Weights of partially
#'   overlapping genes are downscaled by the proportion overlapping the
#'   excluded region. Defaults to 0.1 Mbp or 0.1 cM (depending on genetic
#'   coordinates). Range \code{(0, Inf]}; \code{0} denotes standard gene
#'   boundaries and \code{Inf} excludes all genes on the same chromosome.
#'   Ignored if appropriate covariate columns are not detected in
#'   \code{obj.info}.
#'
#' @param kernel_function
#'   \code{character}; kernel function used to generate probability weights
#'   from distances between the focal gene and other genes in confounder
#'   space. Default is \code{"normal"} (Gaussian kernel). Alternate options
#'   include \code{"exponential"} and \code{"inverse"}. Ignored if covariate
#'   columns are not detected in \code{obj.info}.
#'
#' @param kernel_scale
#'   \code{numeric}; scalar used in the kernel function to convert
#'   Mahalanobis distances to probabilities. Defaults to \code{2} for
#'   Gaussian decay, \code{log(10)} for exponential decay, or \code{2} for
#'   inverse decay. Range \code{(0, Inf]}. Ignored if covariate columns are
#'   not detected in \code{obj.info}.
#'
#' @param kernel_weight_max
#'   \code{numeric}; maximum probability weight for a single gene. Defaults
#'   to \code{0.05}. Must be in the range \code{(1 / (n.genes - 1), 1]}.
#'   Set to \code{1} if no upper bound is desired. Ignored if covariate
#'   columns are not detected in \code{obj.info}.
#'
#' @param gpd_cutoff
#'   \code{numeric}; threshold tail probability at which to apply GPD tail
#'   fitting. Defaults to \code{500 / n_permutations}. Must be in the range
#'   \code{[max(c(1e-04, 500 / n_permutations)), 0.05]}. The lower bound constraint
#'   ensures that a minimum of \code{500} exceedances are available for GPD
#'   estimation, while also ensuring compatibility with the empirical CDF, where
#'   the lowest evaluated quantile is \code{1e-4}.
#'
#' @param seed
#'   \code{integer}; random seed for reproducibility. Preserved across
#'   subsequent polylinkR functions (\code{permute_polylinkr_data}, \code{rescale_polylinkr_data}).
#'   Defaults to \code{NULL}, in which case a seed is generated
#'   automatically. Must be within
#'   \code{[-.Machine$integer.max, .Machine$integer.max]}.
#'
#' @param verbose
#'   \code{logical}; should progress messages be printed to the console?
#'   Defaults to \code{TRUE}.
#'
#' @param n_cores
#'   \code{integer}; number of cores for parallel processing. Defaults to
#'   \code{1} or \code{maximum cores - 1}. Must be in the range
#'   \code{[1, maximum cores]}.
#'
#' @param future_plan
#'   \code{character}; parallel backend from the \code{future} package.
#'   Defaults to user \code{n_cores} choice or checks available cores,
#'   choosing \code{"sequential"} on single-core and \code{"multisession"}
#'   on multi-core systems. Options: \code{"multisession"},
#'   \code{"multicore"}, \code{"cluster"}, \code{"sequential"}.
#'
#' @export
#'
#' @import data.table
#' @import foreach
#' @import doFuture
#' @importFrom future plan
#' @importFrom progressr handlers progressor with_progress
#' @importFrom dqrng generateSeedVectors dqset.seed dqsample.int
#' @importFrom Rfast colsums colCumSums rowsums rowmeans
#' @importFrom distances distances distance_columns
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
#' Deconfounded gene scores are recorded in \code{obj.info} as
#' \code{objStat.std}. Set scores are recorded in \code{set.info} as
#' \code{setScore.std}. Results from enrichment tests are in
#' \code{setScore.std.p}.
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
#' \code{attributes(X)$plr_data$permute.data}, where \code{X} is the object
#' name. Similarly, the arguments used in \code{plR_read} are in
#' \code{attributes(X)$plr_args$read.args}.
#'
#' The primary data structure of the \code{plR} object can be accessed using
#' \code{print()} or by simply typing the object's name.
#'
#' @examples
#' \dontrun{
#' # Assuming `my_plr` is the result of `polylinkR::read_polylinkr_data`
#'
#' # Example 1: Basic usage
#' new_plr <- permute_polylinkr_data(plr_input = my_plr)
#'
#' # Example 2: Less permutations, more bootstraps
#' new_plr <- permute_polylinkr_data(
#'    plr_input       = my_plr,
#'    n_permutations  = 1e5,
#'    n_bootstraps    = 100
#')
#'
#' # Example 3: Modified covariate handling, single processor
#' new_plr <- permute_polylinkr_data(
#'   plr_input        = my_plr,
#'   kernel_weight_max = 0.2,
#'   md_method        = "raw",
#'   n_cores          = 1
#' )
#'
#' # Example 4: Modified GPD estimation, user-specified seed
#' new_plr <- permute_polylinkr_data(
#'   plr_input  = my_plr,
#'   gpd_cutoff = 0.01,
#'   seed       = 1000
#' )
#'
#' # Example 5: Only deconfound scores (no enrichment analysis)
#' new_plr <- permute_polylinkr_data(
#'    plr_input = my_plr,
#'    permute   = FALSE
#' )
#'
#' # Example 6: Use deconfounded scores from Example 5 to run enrichment
#' new_plr <- permute_polylinkr_data(plr_input = new_plr)
#' }
permute_polylinkr_data <- function(plr_input, permute = TRUE, n_permutations = 5e5L,
                                    n_bootstraps = 30L, alternative = "upper",
                                    md_method = "robust", kernel_weight_max = 0.05,
                                    kernel_boundary = "auto", kernel_function = "normal",
                                    kernel_scale = "auto", gpd_cutoff = 5e3L / n_permutations,
                                    seed = NULL, verbose = TRUE, n_cores = "auto",
                                    future_plan = "auto") {

   ##=========================================================================##
   ## PART 1: Clean data and run checks
   ##=========================================================================##

   # Map snake_case parameters to legacy names for internal use
   plR.input <- plr_input
   n.perm <- n_permutations
   n.boot <- n_bootstraps
   alt <- alternative
   md.meth <- md_method
   kern.wt.max <- kernel_weight_max
   kern.bound <- kernel_boundary
   kern.func <- kernel_function
   kern.scale <- kernel_scale
   gpd.cutoff <- gpd_cutoff
   n.cores <- n_cores
   fut.plan <- future_plan

   # track function run time
   startT <- Sys.time()

   # control message reporting
   rlang::local_options(verbose = verbose, .frame = environment())

   # perform checks and unpack required plR objects
   .check_plr_object(f = "permute", ENV = environment())
   .check_arguments(f = "permute", ENV = environment())
   .unpack_plr_object(plr = plR.input, ENV = environment())
   rm(plR.input); gc(verbose = FALSE)

   # set random seed
   seed <- .create_seed(user.seed = seed)
   set.seed(seed)

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

   # announce function
   if (permute) {
      if (perm.path == "full") {
         hdr <- c("deconfound gene scores &",
                  "perform gene set enrichment tests")
      } else { # perm.path = "skip"
         hdr <- "perform gene set enrichment tests"
      }
   } else {
      hdr <- "deconfound gene scores"
   }

   hdr[1] <- paste("Running permute_polylinkr_data:", hdr[1])
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

   ##=========================================================================##
   ## PART 2: Calculate standardised gene scores
   ##=========================================================================##

   if (perm.path == "full") { # perform local regression
      .verbose_msg(cli::style_italic(paste("Estimating prognostic gene scores",
                                   "(confounder effect) using local quadratic",
                                   "regression:\n")))
      .verbose_msg("Calculating pairwise Mahalanobis Distances for all genes")

      cv.val <- as.matrix(obj.info[, .SD, .SDcols = cov.names]) # covariate matrix

      if (md.meth == "robust") {
         .verbose_msg(" in robust covariate space (coverting to normal scores)\n")
         cv.val <- Rfast::colRanks(cv.val, method = "average") # generate ranks
         cv.val <- qnorm((cv.val - 0.5) / n.genes) # convert ranks to normal scores
      } else {
         .verbose_msg(" in covariate space\n")
      }

      # check for singularity among covariates
      if(ncol(cv.val) > 1) {
         qrX <- qr(cv.val)
         if (qrX$rank < length(cov.names)) {
            col.keep <- qrX$pivot[1:qrX$rank]
            cv.val <- cv.val[, col.keep] # remove redundant covariates
            .verbose_msg(paste0("Covariate matrix is singular\nRedundant covariates [",
                        paste0("Cov", setdiff(1:ncol(cov.mat), col.keep),
                               collapse = ", "),
                        "] are ommited from all regression models\n"))
         } else {
            col.keep <- 1:length(cov.names) # no redundancy among covariates
         }
      }

      # generate Mahalanobis distances
      M0 <- distances::distances(data = cv.val, normalize = "mahalanobize")

      .verbose_msg(paste("Estimating regression coefficients for each gene using",
                 kern.func, "kernel weights \n"))
      if (kern.scale == "auto") {
         wK <- which(c("normal", "exponential", "inverse") == kern.func)
         kern.scale <- c(2, log(10), 2)[wK]
      }
      .verbose_msg("\n")

      # determine boundary for excluding genes
      if (kern.bound == "auto") {
         kern.bound <- 1e5
      } else {
         if (!is.infinite(kern.bound)) { # check gene boundaries and rescale
            max.chr <- max(obj.info[, max(endpos) - min(startpos), by = chr]$V1)
            kern.bound <- ifelse(max.chr < kern.bound, Inf, kern.bound)
         }
      }

      # break job into chunks
      n.ch <- ceiling(n.genes / 200)
      os0 <- obj.info$objStat # extract gene scores
      ov0 <- obj.info[, .(X = objID, C = chr, A = startpos.base,
                          B = endpos.base)] # for estimating overlaps (use base boundaries)
      data.table::setkey(ov0, C, A, B)
      if (md.meth == "robust") { # reset covariate values
         cv.val <- as.matrix(obj.info[, .SD, .SDcols = cov.names])
      }
      cv.var <- t(chol(stats::var(cv.val))) # for calculating distances between obs and exp covariates

      fopt <- list(packages = c("data.table", "Rfast"),
                   globals = c("M0", "kern.bound",  "kern.scale", "kern.func",
                               "kern.wt.max", "n.cov", "n.genes", "os0", "ov0",
                               "cv.val", "cv.var", "prog", ".fit_lqr",
                               ".cap_probs", ".md2p"), seed = seed)

      progressr::with_progress({
         prog <- progressr::progressor(steps = n.ch)
         suppressPackageStartupMessages(
            ss <- foreach::foreach(I = split(1:n.genes, sort(1:n.genes %% n.ch)),
                                   .options.future = fopt) %dofuture% {
               wt.i <- t(distances::distance_columns(M0, column_indices = I)) # convert to standard matrix

               # covert distances to kernel weights
               wt.i <- .md2p(wt = wt.i, kern.scale = kern.scale,
                             kern.func = kern.func)

               # rescale distance for overlapping genes
               if (is.infinite(kern.bound)) { # set all genes in same chromosome to 0
                  ovX <- ov0[X %in% I]
                  for (cI in unique(ovX$C)) {
                     wt.i[ovX[C == cI]$X - min(I) + 1, ov0[C == cI]$X] <- 0
                  }
               } else { # set fully overlapping genes to 0, scale others by overlap
                  ovX <- ov0[X %in% I, .(X, C, A = pmax(0, A - kern.bound),
                                         B = B + kern.bound)]
                  oL <- data.table::foverlaps(ov0, ovX, type = "any",
                                              mult = "all", nomatch = NULL) # identify overlaps
                  oL[, xN := pmin(B, i.B) - pmax(A, i.A)] # overlap
                  oL[, x0 := ifelse(xN > 0, 1 - (xN / (i.B - i.A)), 0)] # proportion of non-overlapping sequence
                  I.J <- cbind(oL$X - min(I) + 1, oL$i.X) # coordinates to update
                  wt.i[I.J] <- wt.i[I.J] * oL$x0
               }

               # calculate probability weights
               if (kern.wt.max < 1) { # apply probability cap
                  wt.i <- .cap_probs(mm0 = wt.i, maxP = kern.wt.max)
               } else {
                  wt.i <- wt.i / Rfast::rowsums(wt.i)
               }

               # 0-order (Nadaraya–Watson) estimate
               nw <- wt.i %*% os0

               # estimate local regression coefficients
               lqr <- foreach::foreach(i = I) %do% {
                 d.i <- sweep(cv.val, 2, cv.val[i, ], "-", check.margin = FALSE) # centre covariates
                  .fit_local_quad_reg(d1 = d.i, wt0 = wt.i[i - min(I) + 1, ], object_scores = os0,
                           n.cv = n.cov)
               }

               # calculate fitting statistics
               os.Ne <- 1 / Rfast::rowsums(wt.i ^ 2) # effective sample size
               os.cv <- cv.val[I, ] - wt.i %*% cv.val
               os.MD <- sqrt(Rfast::colsums(forwardsolve(cv.var, t(os.cv)) ^ 2)) # Mahalanobis distance to expected value
               prog()
               list(do.call(cbind, lqr), os.Ne, os.MD, c(nw))
            }
         )
      })

      # extract and label coefficient matrix
      lqr.coeff <- t(do.call(cbind, sapply(ss, "[[", 1)))
      oo <- outer(col.keep, col.keep, FUN = "paste", sep = ".")
      colnames(lqr.coeff) <- c("b.0", paste0("b.", col.keep),
                               paste0("b.", oo[upper.tri(oo, diag = TRUE)]))

      # report local regression outcomes in cases where quadratic model failed
      lqr.fit <- Rfast::rowTrue(!is.na(lqr.coeff))
      if (any(lqr.fit < choose(n.cov, 3) + 1)) {
         tfit <- table(lqr.fit)
         .verbose_msg("Local quadratic regression failed for some genes:\n")
         if (any(lqr.fit == n.cov + 1)) {
            lf <- tfit[names(tfit) == n.par[2]]
            .verbose_msg(paste(lf, "gene", ifelse(lf == 1, "", "s"), "had linear fit"))
         }
         .verbose_msg(ifelse(length(tfit) == 3, " and ", " "))
         if (any(lqr.fit == 1)) {
            nf <- tfit[names(tfit) == n.par[3]]
            .verbose_msg(paste(nf, "gene", ifelse(lf == 1, "", "s"), "had 0-order fit"))
         }
         .verbose_msg("\n[check lqr.fit object in attributes]\n")
      }

      # extract NW estimator
      nw.est <- unlist(lapply(ss, "[[", 4))

      # extract model fit metrics
      obj.Ne <- unlist(lapply(ss, "[[", 2))
      obj.MD <- unlist(lapply(ss, "[[", 3))
      rm(ss); gc(verbose = FALSE)

      # update gene scores
      obj.info[, objStat.res := objStat - lqr.coeff[, 1]] # calculate residual values (observed - prognostic score)
      obj.info[, objStat.std := scale(objStat.res)] # ensure mean = 1 sd = 0

      .verbose_msg("Generating deconfounded gene scores and model fitting statistics\n")
      # calculate correlations between gene scores and covariates
      cv.val <- cbind(obj.info[, .(objStat, objStat.std)], cv.val)
      cP <- stats::cor(cv.val, method = "pearson")
      cS <- stats::cor(cv.val, method = "spearman")
      all.cls <- c("objStat", "objStat.std", cov.names)
      os.cov.cor <- data.table::data.table(t(utils::combn(all.cls, 2)),
                                           c(cP[lower.tri(cP)]),
                                           c(cS[lower.tri(cS)]))
      data.table::setnames(os.cov.cor, c("Var.1", "Var.2", "Rho.p", "Rho.s"))

      # compile fit statistics
      set.Ne <- set.obj[, mean(obj.Ne[objID]), keyby = setID]$V1
      set.MD <- set.obj[, mean(obj.MD[objID]), keyby = setID]$V1
      sO <- data.table::data.table(Class = "Gene", ID = obj.info$objID.orig,
                                   N.eff = obj.Ne, Mahal.dist = obj.MD)
      sS <- data.table::data.table(Class = "Gene set", ID = set.info$setID.orig,
                                   N.eff = set.Ne, Mahal.dist = set.MD)
      kern.fit.stats <- rbind(sO, sS)
      no.deconf <- FALSE
      rm(obj.Ne, obj.MD, set.Ne, set.MD, sO, sS); gc(verbose = FALSE)
   } else { # no regression required
      if (perm.path == "partial") { # no prognostic scores present
         .verbose_msg(paste("No covariate columns identified:",
                    "Skipping gene score deconfounding step\n"))
         # update gene scores
         obj.info[, objStat.res := objStat - mean(objStat)] # calculate residual values
         obj.info[, objStat.std := objStat.res / sd(objStat)] # ensure mean = 0 sd = 1
         no.deconf <- TRUE
         lqr.coeff <- NULL
         kern.wt.max <- NULL
         kern.scale <- NULL
         M0 <- NULL
         os.cov.cor <- NULL
         kern.fit.stats <- NULL
      }
   }

   ##=========================================================================##
   ## PART 3: Genereate null gene set score quantiles and gpd parameters
   ##=========================================================================##

   if (permute) {
      .verbose_msg(cli::style_italic("\nGenerating null set score distributions:\n"))
      .verbose_msg(paste0("Estimating empirical CDFs and fitting GPDs using ",
                  n.boot, " x ", n.perm / ifelse(n.perm / 1e6 >= 1, 1e6, 1e3),
                  ifelse(n.perm / 1e6 >= 1, "M", "k"), " permuted gene sets\n"))

      n.max <- max(set.info$setN)

      if (n.sets == 1) {
         th0 <- n.max
      } else { # determine the sub-sampled gene set sizes (interpolate the remainder)
         wM <- which.max(cumsum(1:200) * 10 >= n.max)
         th0 <- cumsum(unlist(lapply(1:wM, rep, 10)))
         th0 <- c(th0[2:sum(th0 < n.max)], n.max) # exclude unused values and cap at max set size
      }

      fpc <- sqrt(th0 * (n.genes - th0) / (n.genes - 1)) # combined variance and finite population correction
      n.th <- length(th0) # no. set sizes estimated
      n.tail <- ceiling(gpd.cutoff * n.perm) # gpd tail length

      # create quantiles to sample: use logistic quantiles to ensure higher precision at margins
      q.bnd <- 10 ^ -(floor(log10(n.perm)) - 1) # constrain quantile bounds to avoid poor estimation
      eQnt <- plogis(seq(qlogis(q.bnd), qlogis(1 - q.bnd), length.out = 500))
      eQnt <- sort(unique(c(1 - gpd.cutoff, eQnt))) # add GPD estimation anchor
      eps <- sqrt(eQnt * (1 - eQnt) / n.perm) # interval around quantiles to estimate measurement error

      # determine optimal block size for permutations
      minB <- max(n.tail, min(n.perm / n.cores, 5e4)) / 1e4 # minimum block size
      possB <- which(n.perm %% (1:5 * 1e4) == 0) # possible block sizes
      n.block <- possB[which.min(abs(minB - possB))] * 1e4 # final block size
      x.th <- 1 - n.tail / n.block # tail threshold in inner permutation

      os0 <- obj.info$objStat.std # gene scores
      if (alt == "lower") { # take extremes from lower tail
         os0 <- -os0
      }
      inI <- list(g0 = matrix(-1e6, nrow = n.th, ncol = n.tail),
                  e0 = replicate(n.th, list())) # inner loop initial combine object
      inO <- list(g0 = NULL,
                  e0 = list(qX = replicate(n.th, NULL),
                            qE = replicate(n.th, NULL))) # outher loop initial combine object
      dqi <- replicate(n.boot,
                       dqrng::generateSeedVectors(n.perm / n.block),
                       simplify = FALSE) # reproducible seeds

       fopt <- list(packages = c("data.table", "foreach", "dqrng", "Rfast"),
                   globals = c("n.genes", "n.max", "n.tail", "n.sets", "n.th",
                               "n.block", "x.th", "fpc", "os0", "th0",  "inI",
                               "dqI", ".cmb"),
                   seed = seed)

      progressr::with_progress({
         prog <- progressr::progressor(steps = n.boot)
         std.null <- foreach::foreach(I = dqi, .combine = .cmb, .init = inO) %do% { # iterate over each bootstrap replicate
            suppressPackageStartupMessages(
               out <- foreach::foreach(x = I, .combine = .cmb, .init = inI,
                                       .options.future = fopt) %dofuture% {
                  dqrng::dqset.seed(x)
                  csX <- replicate(n.block,
                                   os0[dqrng::dqsample.int(n.genes, n.max)]) # faster than dqsample::dqsample
                  if (n.sets == 1) { # only take final entry
                     csX <- matrix(Rfast::colsums(csX) / fpc, nrow = 1) # standardise set scores
                  } else { # cumulative sum
                     csX <- Rfast::colCumSums(csX)[th0, ] / fpc # standardise set scores
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

                  list(g0 = g0, e0 = e0)
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
      rm(std.null); gc(verbose = FALSE)

      # check for inviable fits in gpd estimation
      gpd.miss <- sum(!complete.cases(g0)) / nrow(g0)
      if (gpd.miss > 0.5) {
         warning("GPD fitting failed for >50% of nulls -- halting function:\n",
                 "try again with different gpd fitting options\n",
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
            e0 <- lapply(e0, function(x) -x[, length(eQnt):1])
            eE <- lapply(eE, function(x) x[, length(eQnt):1])
         }

         # calculate means
         e0mean <- sapply(e0, Rfast::colmeans)
         eEmean <- sapply(eE, Rfast::colmeans)

         # fit smoothing functions
         fopt <- list(packages = "mgcv",
                      globals = c("th0", "K", "n.th", "n.max", "sm0",
                                  ".par_smooth", "Predict.matrix.wls.smooth",
                                  "smooth.construct.wls.smooth.spec"),
                      seed = seed)
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

         # calculate associated cutoff p-value gpd minimum non-0 p values and quantiles
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

         rm(g0, gpdM); gc(verbose = FALSE) # clean up
      }
      .verbose_msg("\n")
   } else {
      dqi <- NULL
      ecdf0 <- NULL
      eSumm <- NULL
      gpd0 <- NULL
      gSumm <- NULL
   }

   ##=========================================================================##
   ## PART 4: Calculate test statistics
   ##=========================================================================##

   if (permute) {
      .verbose_msg(cli::style_italic("Calculating test statistics:\n"))
      .verbose_msg(paste0("Calculating observed and null gene set scores\n"))

      # calculate observed standardised set scores
      # mean of standardised gene scores = 0 and variance = 1 -> set score variance = set size (assuming independence)
      ss.obs <- set.obj[, sum(os0[objID]) / sqrt(.N), keyby = setID]
      data.table::setkey(set.info, setID)
      set.info[.(ss.obs), setScore.std := V1]

      # calculate p values
      .verbose_msg("Calculating p values for each gene set\n")
      if (alt == "lower") {
         .get_p <- .get_p_lt
      } else {
         .get_p <- .get_p_ut
      }
      set.info[, setScore.std.p := .get_p(ss = setScore.std, n = setN[1],
                                          g0 = gpd0, e0 = ecdf0), by = setN]
   }

   ##=========================================================================##
   ## PART 5: Create output
   ##=========================================================================##

   # repack old data and add new data
   plr.data <- list()
   read.data <- list(n.genes = n.genes, n.sets = n.sets, n.chr = n.chr,
                     n.cov = n.cov, n.set.genes = n.set.genes,
                     file.paths = file.paths, coord = coord,
                     get.objStat = get.objStat, cov.names = cov.names,
                     cov.info = cov.info, no.share = no.share,
                     pos.info = pos.info)
   plr.data$read.data <- read.data
   permute.data <- list(ecdf.std = ecdf0, gpd.std = gpd0, no.deconf = no.deconf)
   plr.data$permute.data <- permute.data

   # retain arguments and summary statistics
   plr.seed <- list(seed = seed, dqrng.seed = dqi)
   args0 <- names(formals(get("plR_permute", envir = asNamespace("polylinkR"))))
   permute.args <-  mget(setdiff(args0, "plR.input"), envir = environment()) # collect arguments
   plr.args$permute.args <- permute.args

   if ("permute.summary" %in% names(plr.summary)) { # update ecdf and gpd estimation summaries
      plr.summary$permute.summary$gpd.std.summary <- gSumm
      plr.summary$permute.summary$ecdf.std.summary <- eSumm
   } else {
      permute.summary <- list(kern.fit.stats = kern.fit.stats,
                              os.cov.cor = os.cov.cor, mahal.dist = M0,
                              lqr.coeff = lqr.coeff, gpd.std.summary = gSumm,
                              ecdf.std.summary = eSumm)
      plr.summary$permute.summary <- permute.summary
   }

   # capture R session information and record run time
   r0 <- .get_time(start.time = startT)
   if ("permute.session" %in% names(plr.session)) {
      perm.sess <- plr.session$permute.session
      permute.session <- list()
      permute.session$session <- list(perm.sess$session, sessionInfo())
      permute.session$run.time <- c(perm.sess$run.time, r0[[2]])
   } else {
      permute.session <- list(session = sessionInfo(), run.time = r0[[2]])
   }
   plr.session$permute.session <- permute.session

   # set s3 class and create attributes
   .file_reset(OI = obj.info, SI = set.info, SO = set.obj, pos.info = pos.info) # recreate original file formats
   OUT <- list(set.info = set.info, obj.info = obj.info, set.obj = set.obj)
   plr.out <- .new_plr(BASE = OUT, plr_data = plr.data, plr_args = plr.args,
                       plr_summary = plr.summary, plr_seed = plr.seed,
                       plr_session = plr.session)

   .verbose_msg("\n")
   ftr <- cli::col_cyan(paste0("Finished permute_polylinkr_data -- run time: ", r0[[1]]))
   pdg <- (80 - nchar(ftr)) / 2
   .verbose_msg(cli::boxx(ftr, padding = c(0, floor(pdg), 0, ceiling(pdg)),
                  border_style = "double", col = "cyan", border_col = "cyan"))
   .verbose_msg("\n")
   return(invisible(plr.out))
}


#' @title Gene set enrichment with control for confounding covariates (deprecated)
#' @description
#' This function is deprecated. Please use \code{permute_polylinkr_data()} instead.
#' @param plR.input Deprecated. Use \code{plr_input}.
#' @param permute Deprecated. Use \code{permute}.
#' @param n.perm Deprecated. Use \code{n_permutations}.
#' @param n.boot Deprecated. Use \code{n_bootstraps}.
#' @param alt Deprecated. Use \code{alternative}.
#' @param md.meth Deprecated. Use \code{md_method}.
#' @param kern.wt.max Deprecated. Use \code{kernel_weight_max}.
#' @param kern.bound Deprecated. Use \code{kernel_boundary}.
#' @param kern.func Deprecated. Use \code{kernel_function}.
#' @param kern.scale Deprecated. Use \code{kernel_scale}.
#' @param gpd.cutoff Deprecated. Use \code{gpd_cutoff}.
#' @param seed Deprecated. Use \code{seed}.
#' @param verbose Deprecated. Use \code{verbose}.
#' @param n.cores Deprecated. Use \code{n_cores}.
#' @param fut.plan Deprecated. Use \code{future_plan}.
#' @export
plR_permute <- function(plR.input, permute = TRUE, n.perm = 5e5L, n.boot = 30L,
                         alt = "upper", md.meth = "robust", kern.wt.max = 0.05,
                         kern.bound = "auto", kern.func = "normal",
                         kern.scale = "auto", gpd.cutoff = 5e3L / n.perm,
                         seed = NULL, verbose = TRUE, n.cores = "auto",
                         fut.plan = "auto") {
   .Deprecated("permute_polylinkr_data", package = "polylinkR",
               msg = "plR_permute() is deprecated. Use permute_polylinkr_data() instead.")
   permute_polylinkr_data(
      plr_input = plR.input,
      permute = permute,
      n_permutations = n.perm,
      n_bootstraps = n.boot,
      alternative = alt,
      md_method = md.meth,
      kernel_weight_max = kern.wt.max,
      kernel_boundary = kern.bound,
      kernel_function = kern.func,
      kernel_scale = kern.scale,
      gpd_cutoff = gpd.cutoff,
      seed = seed,
      verbose = verbose,
      n_cores = n.cores,
      future_plan = fut.plan
   )
}
