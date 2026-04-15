#' @title Read and validate required files for polylinkR
#'
#' @description
#' Reads and validates the input files for a \code{polylinkR} analysis. Handles
#' file path resolution, column name canonicalisation, data validation, and
#' optional filtering and data generation tasks like oordinate conversion or
#' gene score calculation.
#'
#' @param input.path
#'   \code{character}; path to the directory with input files.
#'   Defaults to \code{NULL}. Not required if using separate file paths.
#'   Automatically searches for required files with the following labels:
#'   \itemize{
#'     \item \code{setinfo}: a data frame with set-level information.
#'     \item \code{objinfo}: a data frame with object-level information.
#'     \item \code{setobj}: a data frame mapping objects to sets.
#'   }
#'   It also searches for optional files:
#'   \itemize{
#'     \item \code{recrate}: an optional file for genetic coordinate
#'       conversion.
#'     \item \code{varinfo}: an optional file for gene score generation.
#'   }
#'   Allowable file names include capitalisation of each word in the label
#'   and use of internal separators (e.g., \code{set.info}, \code{set_info},
#'   \code{SetInfo}, \code{Set.Info}, \code{Set_Info}). Also see the section
#'   on required and optional input files for the columns needed in each
#'   input file.
#'
#' @param obj.info.path
#'   \code{character}; path to the \code{obj.info} file. Defaults to
#'   \code{NULL}. Ignored if \code{input.path} is provided; otherwise all
#'   required file paths must be specified.
#'
#' @param set.info.path
#'   \code{character}; path to the \code{set.info} file. Defaults to
#'   \code{NULL}. Ignored if \code{input.path} is provided; otherwise all
#'   required file paths must be specified.
#'
#' @param set.obj.path
#'   \code{character}; path to the \code{set.obj} file. Defaults to
#'   \code{NULL}. Ignored if \code{input.path} is provided; otherwise all
#'   required file paths must be specified.
#'
#' @param var.info.path
#'   \code{character}; path to the \code{var.info} file. Defaults to
#'   \code{NULL}. Optional file used for gene score generation. Ignored if
#'   \code{input.path} is provided.
#'
#' @param rec.rate.path
#'   \code{character}; path to the \code{rec.rate} file. Defaults to
#'   \code{NULL}. Optional file used for genetic coordinate conversion.
#'   Ignored if \code{input.path} is provided.
#'
#' @param min.set.n
#'   \code{integer}; minimum size of gene sets to be retained. Defaults to
#'   \code{2}. Must be in the range \code{[2L, max.set.n)}.
#'
#' @param max.set.n
#'   \code{integer}; maximum size of gene sets to be retained. Defaults to
#'   \code{Inf}. Must be in the range \code{(min.set.n, Inf)}.
#'
#' @param group
#'   \code{character}; label used to identify input files within a directory.
#'   Defaults to \code{NULL}.
#'
#' @param map.fun
#'   \code{character}; mapping function to convert physical to genetic
#'   distances. Options are \code{"Haldane"}, \code{"Kosambi"} (default),
#'   \code{"Carter-Falconer"}, and \code{"Morgan"}.
#'
#' @param obj.buffer
#'   \code{numeric}; interval around genes (in base pairs) to include when
#'   assigning values from the \code{var.info} file. Defaults to \code{1e4} if
#'   \code{var.info} is provided and user does not set a value, otherwise it is
#'   set to 0 if score assignment is not performed. User values must be in the
#'   range \code{[0, 1e5L]}. Note that if the user provides their own gene
#'   scores in the \code{obj.info } input file (i.e. not computed from the
#'   \code{var.info} file), then the start and end positions must include any
#'   buffer used to bin scores, otherwise polylinkR deconfounding and
#'   autocorrelation inference will not be performed appropriately.
#'
#' @param obj.stat.fun
#'   \code{character}; function used to correct maximum gene scores based on
#'   the number of overlapping summary statistics (SNPs or windows). Default is
#'  `non.param`, a robust non-parametric method that uses binned data to
#'  calculate median and median absolute deviation (MAD) to normalise scores.
#'  Alternatively, `lm.logN` applies a linear regression to the log-transformed
#'  SNP / bin counts (assumes a roughly linear relationship is appropriate). In
#'  both cases, expected scores are estimated and gene scores calculated as the
#'  residual value. Ignored if no var.info file is provided.
#'
#' @param bin.size
#'   \code{integer}; gene set size interval for non-parametric correction.
#'   Defaults to \code{250L}. Must be in the range \code{[50L, 1e3L]};
#'   ignored if the parametric function is used.
#'
#' @param obj.in
#'   \code{character} or \code{numeric} vector; \code{objID}s of genes to
#'   explicitly retain. Defaults to \code{NULL}.
#'
#' @param obj.out
#'   \code{character} or \code{numeric} vector; \code{objID}s of genes to
#'   explicitly remove. Defaults to \code{NULL}.
#'
#' @param set.in
#'   \code{character} or \code{numeric} vector; \code{setID}s of gene sets
#'   to explicitly retain. Defaults to \code{NULL}.
#'
#' @param set.out
#'   \code{character} or \code{numeric} vector; \code{setID}s of gene sets
#'   to explicitly remove. Defaults to \code{NULL}.
#'
#' @param set.merge
#'   \code{numeric}; minimum proportion of shared genes for merging gene sets.
#'   Defaults to \code{0.95}. Must be in the range \code{(0, 1]}.
#'
#' @param rem.genes
#'   \code{logical}; should genes with identical genomic positions be
#'   removed? Defaults to \code{FALSE}.
#'
#' @param verbose
#'   \code{logical}; should progress messages be printed to the console?
#'   Defaults to \code{TRUE}.
#'
#' @export
#'
#' @import data.table
#' @import foreach
#' @importFrom rlang local_options
#' @importFrom cli boxx col_cyan style_italic
#'
#' @section Required input file structure:
#' The following three files are compulsory for a \code{polylinkR} analysis
#' and extend the format introduced in Polysel. All files must be comma-
#' separated (\code{.csv}) or tab-separated (\code{.tsv}), with a header
#' (see examples at \url{https://github.com/CMPG/polysel/tree/master/data}).
#'
#' \describe{
#'   \item{\code{set.info}}{
#'     A \code{data.table} (and \code{data.frame}) with gene set information.
#'     \itemize{
#'       \item \strong{setID}: \code{character} or \code{factor} vector of
#'         unique gene set identifiers. Required.
#'     }
#'   }
#'   \item{\code{obj.info}}{
#'     A \code{data.table} (and \code{data.frame}) with gene (objects)
#'     information.
#'     \itemize{
#'       \item \strong{objID}: \code{character} or \code{factor} vector of
#'         unique gene identifiers. Required.
#'       \item \strong{objStat}: optional \code{numeric} vector of
#'         pre-computed gene scores. Required if \code{var.info} is absent;
#'         otherwise computed from \code{var.info}.
#'       \item \strong{CovX}: \code{numeric} vector of covariate scores,
#'         where \code{X} is a positive integer denoting covariate number.
#'         Required for deconfounding gene scores (\code{objStat}) in
#'         \code{plR_permute}.
#'       \item \strong{chr}: \code{character} or \code{numeric} chromosome /
#'         contig labels. Required for gene score deconfounding and gene set
#'         score decorrelation in and \code{plR_permute} and \code{plR_rescale},
#'         respectively.
#'       \item \strong{startpos}: \code{numeric} gene start position in base
#'         pairs. If the user is providing their own gene scores, then this
#'         position represents the original position minus any buffer used to
#'         assign summary statistics to genes. Required for gene score
#'         deconfounding and gene set score decorrelation in and
#'         \code{plR_permute} and \code{plR_rescale}, respectively.
#'       \item \strong{endpos}: \code{numeric} gene end position in base pairs.
#'         If the user is providing their own gene scores, then this position
#'         represents the original position plus any buffer used to assign
#'         summary statistics to genes. Required for gene score deconfounding
#'         and gene set score decorrelation in and
#'         \code{plR_permute} and \code{plR_rescale}, respectively.
#'       \item \strong{startpos.base}: \code{numeric} gene start position in
#'         base pairs. Default lower gene boundary (ignoring buffer used in gene
#'         assignment). Created internally if absent or gene scores are
#'         estimated (i.e. var.info provided). Required for gene score
#'         deconfounding and set score decorrelation in \code{plR_permute}
#'         and \code{plR_rescale}, respectively.
#'       \item \strong{endpos.base}: \code{numeric} gene end position base
#'         pairs. Default upper gene boundary (ignoring buffer used in gene
#'         assignment). Created internally if absent or gene scores are
#'         estimated (i.e. var.info provided). Required for gene score
#'         deconfounding and set score decorrelation in \code{plR_permute}
#'         and \code{plR_rescale}, respectively.
#'     }
#'   }
#'   \item{\code{set.obj}}{
#'     A \code{data.table} (and \code{data.frame}) mapping genes to gene sets.
#'     \itemize{
#'       \item \strong{setID}: \code{character} or \code{factor} vector of
#'         unique gene set identifiers. Required.
#'       \item \strong{objID}: \code{character} or \code{factor} vector of
#'         unique gene identifiers. Required.
#'     }
#'   }
#' }
#'
#' @section Optional input file structure:
#' These optional comma-separated (\code{.csv}) or tab-separated (\code{.tsv})
#' files provide the following \code{plR_read} functionality:
#'
#' \describe{
#'   \item{\code{var.info}}{
#'     Contains information used to compute gene scores (\code{objStat}).
#'     \itemize{
#'       \item \strong{chr}: \code{character} or \code{numeric} chromosome /
#'         contig labels. Required.
#'       \item \strong{pos}: \code{numeric} position where statistic was
#'         evaluated (e.g., SNP or central point in window). Required.
#'       \item \strong{value}: \code{numeric} statistical value from SNP /
#'         window score. Required.
#'     }
#'   }
#'   \item{\code{rec.rate}}{
#'     Used to transform genetic coordinates from physical to genetic
#'     distances. Requires \code{startpos} and \code{endpos} in \code{obj.info}
#'     to be base pair coordinates. Uses HapMap format (see labelled examples
#'     at \url{https://zenodo.org/records/11437540}).
#'     \itemize{
#'       \item \strong{chr}: \code{character} or \code{numeric} chromosome /
#'         contig labels. Required.
#'       \item \strong{pos}: \code{numeric} base pair position of upstream
#'         marker for the recombination interval. Required.
#'       \item \strong{rate}: \code{numeric} recombination rate in cM per bp
#'         in the downstream interval. Required.
#'       \item \strong{map}: \code{numeric} genetic distance in centiMorgans
#'         (cM). Optional; will be calculated from recombination rates if
#'         absent.
#'     }
#'   }
#' }
#'
#' @return
#' A \code{plR} S3 object containing three complementary datasets:
#' \describe{
#'   \item{\code{obj.info}}{
#'     \code{data.table} containing information for each gene (object).
#'   }
#'   \item{\code{set.info}}{
#'     \code{data.table} containing information for each gene set.
#'   }
#'   \item{\code{set.obj}}{
#'     \code{data.table} containing the mapping of genes to gene sets.
#'   }
#' }
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
#' # Example 1: Read all files from a single folder
#' my_plr <- plR_read(input.path = "path/to/files")
#'
#' # Example 2: Read all files with "POP1" in label from a single folder
#' my_plr <- plR_read(
#'    input.path = "path/to/files",
#'    group      = "POP1"
#' )
#'
#' # Example 3: Relax gene set merging criteria
#' my_plr <- plR_read(
#'    input.path = "path/to/files",
#'    set.merge  = 0.50
#' )
#'
#' # Example 4: Specify separate file paths and remove user-specified sets
#' my_plr <- plR_read(
#'   set.info.path = "path/to/set.info",
#'   set.obj.path  = "path/to/set.obj",
#'   obj.info.path = "path/to/obj.info",
#'   set.out       = c("set1", "set2")
#' )
#'
#' # Example 5: Generate gene scores from var.info using regression
#' # and include a 50 kb buffer around genes
#' my_plr <- plR_read(
#'   set.info.path = "path/to/set.info",
#'   set.obj.path  = "path/to/set.obj",
#'   obj.info.path = "path/to/obj.info",
#'   var.info.path = "path/to/var.info",
#'   obj.stat.fun   = "lm.logN", obj.buffer = 5e4
#' )
#'
#' # Example 6: Convert distances to cM using rec.rate file
#' my_plr <- plR_read(
#'   set.info.path = "path/to/set.info",
#'   set.obj.path  = "path/to/set.obj",
#'   obj.info.path = "path/to/obj.info",
#'   rec.rate.path = "path/to/rec.rate"
#' )
#' }
plR_read <- function(input.path = NULL, obj.info.path = NULL,
                     set.info.path = NULL, set.obj.path = NULL,
                     var.info.path = NULL, rec.rate.path = NULL, min.set.n = 2L,
                     max.set.n = Inf, group = NULL, map.fun = "kosambi",
                     obj.buffer = "auto", obj.stat.fun = "non.param",
                     bin.size = 250L, obj.in = NULL, obj.out = NULL,
                     set.in = NULL, set.out = NULL, set.merge = 0.95,
                     rem.genes = FALSE, verbose = TRUE) {

   # track function run time
   startT <- Sys.time()

   # control message reporting
   rlang::local_options(verbose = verbose, .frame = environment())

   # check arguments
   .arg_check(f = "read", ENV = environment())

   # announce function
   hdr <- "Running plR_read -- Loading polylinkR input files"
   pdg <- (80 - nchar(hdr)) / 2
   .vrb(cli::boxx(hdr, padding = c(0, floor(pdg), 0, ceiling(pdg)),
                  border_style = "double", col = "cyan", border_col = "cyan"))
   .vrb("\n\n")

   #---------------------------------------------------------------------------#
   ## Read and check plR input files
   #---------------------------------------------------------------------------#

   .path_check(input.path = input.path, oi.path = obj.info.path,
               si.path = set.info.path, so.path = set.obj.path,
               vi.path = var.info.path, rr.path = rec.rate.path,
               group = group, ENV = environment())

   .vrb(cli::style_italic("Reading files:\n"))
   for (fni in names(file.paths)) {
      fpi <- file.paths[[fni]]
      .vrb(paste0("--> Loading ", fni, " (",  fpi, ")\n"))
      assign(x = fni, value = data.table::fread(fpi), envir = environment())
   }

   #---------------------------------------------------------------------------#
   ## Check presence of canonical columns
   #---------------------------------------------------------------------------#

   .vrb(cli::style_italic("\nPerforming file checks\n"))

   # check appropriate columns are present and no missing data
   CHK0 <- c("obj.info", "set.info", "set.obj")
   CHK1 <- list(c("obj_ID"), "set_ID", c("set_ID", "obj_ID"))
   CHK2 <- list(names(obj.info), names(set.info), names(set.obj))
   canon <- list(c("objID"), "setID", c("setID", "objID"))

   bp.to.cM <- exists("rec.rate") # determine if rec.rate file has been loaded
   get.objStat <- exists("var.info")

   if (get.objStat | bp.to.cM) {
      CHK1[[1]] <- c(CHK1[[1]], "chr", "startpos", "endpos") # include genomic coordinates
      canon[[1]] <- c(canon[[1]], "chr", "startpos", "endpos")
   }

   if (get.objStat) {
      CHK0 <- c(CHK0, "var.info")
      CHK1 <- c(CHK1, list(c("chr", "pos", "value")))
      CHK2 <- c(CHK2, list(names(var.info)))
      canon <- c(canon, list(c("chr", "pos", "value")))
   } else { # required if var.info not provided (user has computed their own gene scores)
      CHK1[[1]] <- c(CHK1[[1]], "objStat")
      canon[[1]] <- c(canon[[1]], "objStat")
   }

   if (bp.to.cM) {
      CHK0 <- c(CHK0, "rec.rate")
      CHK1 <- c(CHK1, list(c("chr", "pos", "rate")))
      CHK2 <- c(CHK2, list(names(rec.rate)))
      canon <- c(canon, list(c("chr", "pos", "rate")))
   }

   c0 <- foreach::foreach(i = seq_along(CHK0)) %do% { # loop over different input files
      s0 <- strsplit(CHK1[[i]], "_")
      col.names <- rep(NA, length(s0))
      for(j in seq_along(s0)) { # check required columns
         x <- s0[[j]]
         if (length(x) > 1) { # check multiple options for labels with two components
            for (k in c("", ".", "_")) { # check for alternative column labels
               if (x[2] != "ID") {
                  for (l in 1:3) {
                     if (l %in% 1:2) { # allow upper and lower case for second label component
                        if (l == 1) {
                           x2 <- tolower(substr(x[2], 1, 1))
                        } else {
                           x2 <- toupper(substr(x[2], 1, 1))
                        }
                        x0 <- paste0(x[1], k, x2, substr(x[2], 2, 10))
                     } else {
                        if (x[1] %in% c("start", "end")) { # allow for start and end as viable options
                           x0 <- x[1]
                        } else {
                           x0 <- NA
                        }
                     }
                     if (x0 %in% CHK2[[i]]) col.names[j] <- x0
                  }
               } else {
                  x0 <- paste(x, collapse = k)
                  if (x0 %in% CHK2[[i]]) col.names[j] <- x0
               }
            }
         } else { # label has only single component
            if (x %in% CHK2[[i]]) col.names[j] <- x
         }
      }
      col.names
   }

   col.checks <- sapply(c0, function(x) any(is.na(x)))
   if (any(col.checks)) {
      e.mss <- foreach::foreach(i = which(col.checks), .combine = c) %do% {
         paste0(CHK0[[i]], " [",
                paste(canon[[i]][is.na(c0[[i]])], collapse = ", "), "]")
      }
      stop("Missing required data column in ",
           paste0(e.mss, collapse = " and "),
           "\nSee ?polylinkR::plR_read for required columns", call. = FALSE)
   }

   # record user defined columns and replace non-canonical labels
   for (i in seq_along(CHK0)) {
      if (!all(c0[[i]] %in% canon[[i]])) { # check for noncanonical labels
         .vrb(paste0("Replacing non-canonical column labels for ",
                     CHK0[[i]], "\n"))
         data.table::setnames(x = get(CHK0[[i]]), old = c0[[i]],
                              new = canon[[i]])
      }
   }

   # ensure class consistency across shared data objects
   if (class(set.info$setID) != class(set.obj$setID)) {
      set.info[, setID := as.character(setID)]
      set.obj[, setID := as.character(setID)]
   }

   if (class(obj.info$objID) != class(set.obj$objID)) {
      obj.info[, objID := as.character(objID)]
      set.obj[, objID := as.character(objID)]
   }

   # determine presence of auxiliary columns
   cov.info <- any(grepl(c("COV|Cov|cov"), names(obj.info)))
   pos.info <- all(c("chr", "startpos", "endpos") %in% names(obj.info))

   # apply buffers if gene score assignment chosen otherwise check user option
   if (get.objStat) { # include buffer to gene boundaries
      if (obj.buffer == "auto") {
         obj.buffer <- 1e4 # capture signals 10 kbp on either side of gene margins
      }
      obj.info[, startpos.base := startpos] # retain base lower boundary
      obj.info[, endpos.base := endpos] # retain base upper boundary
      obj.info[, startpos := pmax(startpos - obj.buffer, 0)] # update lower boundary
      obj.info[, endpos := endpos + obj.buffer] # update upper boundary
   } else { # user has supplied own gene scores
      if (pos.info) { # position information present, check for base position info
         bNames <- lapply(c(".", "_"), function(i) {
            paste0(c("startpos", "endpos"), i, "base")
         })
         bCols <- sapply(bNames, function(x) {
            all(x %in% names(obj.info))
         })
         if (any(bCols)) { # required columns exist
            bNI <- bNames[[which(bCols)]]
            ckC <- obj.info[, lapply(.SD, is.numeric), .SDcols = bNI]
            if (!all(ckC)) { # ensure columns are numeric
               stop("User provided ", paste(bNI[ckC], collapse = " and "),
                    " columns in obj.info are not numeric\n", call. = FALSE)
            }
            if (bCols[2]) { # ensure canonical naming conventions
               data.table::setnames(obj.info, old = bNames[[2]], new = bNames[[1]])
            }
         } else { # add if user has not included base position columns
            .vrb(paste0("Note: ", paste(bNames[!bCols], collapse = " and "),
                        " columns not detected in obj.info\n",
                        "[polylinkR analyses will assume standard gene ",
                        " boundaries were used in gene score assignment]"))
            obj.info[, startpos.base := startpos] # retain base lower boundary
            obj.info[, endpos.base := endpos] # update base upper boundary
         }
         obj.buffer <- 0
      }
   }
   obj.info[, midpos := (startpos.base + endpos.base) / 2] # add midpos position

   #---------------------------------------------------------------------------#
   ## Check row contents
   #---------------------------------------------------------------------------#

   # check for shared chromosome labels
   w.mss <- NULL
   oc <- unique(obj.info$chr)
   if (bp.to.cM) {
      rc <- unique(rec.rate$chr)
      worc <- oc %in% rc
      if (!all(worc)) {
         w.mss <- c(w.mss, paste0("Some obj.info chromosomes [",
                                  paste0(oc[!worc], collapse = ", "),
                                  "] are not present in rec.rate\n"))
      }
   }

   if (get.objStat) {
      vc <- unique(var.info$chr)
      wovc <- oc %in% vc
      if (!all(wovc)) {
         w.mss <- c(w.mss, paste0("Some obj.info chromosomes [",
                                  paste0(vc[!wovc], collapse = ", "),
                                  "] are not present in var.info\n"))
      }
   }

   if (!is.null(w.mss)) {
      stop(paste(w.mss, collapse = "\n"),
           "\nAll obj.info chromosomes must be present in auxillary files",
      call. = FALSE)
   }

   if (cov.info) { # ensure that the covariate columns follow conventions
      cov.names <- names(obj.info)[grep("Cov|COV|cov", names(obj.info))]
      cov.int <- as.numeric(substr(cov.names, 4, 4)) # check for appropriate labels
      n.cov <- length(cov.names)
      if (!all(cov.int %in% 1:n.cov)) {
         stop("Covariate column labels in obj.info must be numbered from 1 to",
              " no. of covariates", call. = FALSE)
      } else { # check for numeric data
         cov.num <- obj.info[, lapply(.SD, is.numeric), .SDcols = cov.names]
         wC <- which(!cov.num)
         if (length(wC) > 0) {
            stop("Covariate values in obj.info must be numeric", call. = FALSE)
         }
      }
   } else {
      cov.names <- NULL
      n.cov <- 0
   }

   # check for incomplete rows
   rW <- NULL
   rE <- NULL
   removed.ID <- NULL
   for (i in seq_along(CHK0))  {
      fi <- CHK0[[i]]
      rNA <- !stats::complete.cases(get(fi)[, .SD, .SDcols = canon[[i]]]) # extract required columns
      if (any(rNA)) {
         if (all(rNA)) {
            rE <- c(rE, paste("No complete rows identified in", fi))
         } else {
            rm.rows <- get(fi)[rNA]
            rW <- c(rW, paste(sum(rNA), "incomplete rows detected and removed from",
                              fi))
            # remove incomplete rows
            if (fi == "obj.info") {
               rID <- data.table::data.table(file = fi, setID = NA,
                                             objID = rm.rows$objID,
                                             flag = "incomplete")
            } else {
               if (fi == "set.info") {
                  rID <- data.table::data.table(file = fi,
                                                setID = rm.rows$setID,
                                                objID = NA, flag = "incomplete")
               } else {
                  if (fi == "set.obj") {
                     rID <- data.table::data.table(file = fi,
                                                   setID = rm.rows$setID,
                                                   objID = rm.rows$objID,
                                                   flag = "incomplete")
                  }
               }
            }
            removed.ID <- rbind(removed.ID, rID) # track removed rows
            assign(x = fi, value = get(fi)[!rNA], envir = environment()) # update table
         }
      }
   }

   if (!is.null(rE)) {
      stop(paste0(rE, collapse = "\n"), call. = FALSE)
   }

   # check for shared gene labels
   xov <- set.obj$objID %in% obj.info$objID
   if(!all(xov)) {
      rm.rows <- unique(set.obj[!xov]$objID)
      rW <- c(rW,
              paste(length(rm.rows), "Some set.obj genes [objID] were not ",
                    "found in obj.info: removing orphan genes from set.obj"))
      rID <- data.table::data.table(file = "set.obj", setID = NA,
                                    objID = rm.rows, flag = "orphan")
      removed.ID <- rbind(removed.ID, rID)
      set.obj <- set.obj[xov]
   }

   # check for shared gene set labels
   xsv <- set.obj$setID %in% set.info$setID
   if (!all(xsv)) {
      rm.rows <- unique(set.obj[!xsv]$setID)
      rW <- c(rW,
              paste(length(rm.rows), "Some set.obj gene sets were not found ",
                    "in set.info: removing orphan gene sets from set.obj"))
      rID <- data.table::data.table(file = "set.obj", setID = rm.rows,
                                    objID = NA, flag = "orphan")
      removed.ID <- rbind(removed.ID, rID)
      set.obj <- set.obj[xsv]
   }

   sxv <- set.info$setID %in% set.obj$setID
   if(!all(sxv)) {
      rm.rows <- set.info[!sxv]$setID
      rW <- c(rW,
              paste(length(rm.rows), "Some set.info gene sets were not found ",
                    "in set.obj: removing orphan gene sets from set.info"))
      rID <- data.table::data.table(file = "set.info", setID = rm.rows,
                                    objID = NA, flag = "orphan")
      removed.ID <- rbind(removed.ID, rID)
      set.info <- set.info[sxv]
   }

   if (!is.null(rW)) {
      .vrb(paste0(paste(rW, collapse = "\n"),
                  "\n[info available in removed.ID]\n"))
   }

   # check that all genes and gene sets are unique
   repeat.ID <- NULL
   si.check <- set.info[, .N, by = setID][N > 1]
   Ndup.si <-  nrow(si.check)
   if (Ndup.si > 0) {
      .vrb(paste0("Indentified and removed ", Ndup.si, " duplicated gene set",
                  ifelse(Ndup.si == 1, "", "s"), "\n"))
      rID <- data.table::data.table(file = "set.info",
                                    setID = si.check[, setID],
                                    objID = NA, N = si.check[, N])
      repeat.ID <- rbind(repeat.ID, rID)
      set.info <- set.info[!duplicated(setID)] # remove duplicated gene sets
   }

   oi.check <- obj.info[, .N, by = objID][N > 1]
   Ndup.oi <- nrow(oi.check)
   if (Ndup.oi > 0) {
      .vrb(paste0("Indentified and removed ", Ndup.oi, " duplicated gene",
                  ifelse(Ndup.oi == 1, "", "s"), "\n"))
      rID <- cbind(file = "obj.info", setID = NA, oi.check)
      repeat.ID <- rbind(repeat.ID, rID)
      obj.info <- obj.info[!duplicated(objID)] # remove duplicated genes
   }

   so.check <- set.obj[, .N, by = .(setID, objID)][N > 1]
   Ndup.so <- data.table::uniqueN(so.check$setID)
   if (Ndup.so > 0) {
      .vrb(paste0("Indentified ", Ndup.so, " gene set",
                  ifelse(Ndup.so == 1, "", "s"),
                  " with duplicated gene IDs, removing duplicates\n"))
      rID <- cbind(file.type = "set.obj", so.check)
      repeat.ID <- rbind(repeat.ID, rID)
      set.obj <- set.obj[, .(objID = unique(objID)), by = setID] # remove repeated genes
   }

   if (Ndup.si > 0 | Ndup.oi > 0 | Ndup.so > 0) {
      .vrb("[info available in repeat.ID]\n")
   }

   # check for duplicated genes by position
   if (pos.info) {
      obj.per.set <- obj.info[, .(.N, objID),
                              by = c("chr", "startpos.base", "endpos.base")]
      obj.dups <- obj.per.set[N > 1]
      if (nrow(obj.dups) > 0) {
         .vrb(paste0("Found ", nrow(obj.dups), " genes that have an identical ",
                     "genomic position with other genes\n"))
         if (rem.genes) { # Optional: remove genes with identical start and end positions on same chromosome
            # privilege retention of genes in gene sets (or that appear in most sets)
            obj.dups[, ID := as.factor(paste(chr, startpos.base, endpos.base))]
            obj.dups <- merge(obj.dups[, .(objID, ID)],
                              set.obj[objID %in% obj.dups$objID, .N, by = objID],
                              by = "objID", all.x = TRUE) # enumerate number of times gene occurs in sets
            obj.dups[is.na(N), N := 0]

            obj.dups[, objID.new := objID[which.max(N)], by = ID]
            obj.keep <- unique(obj.dups$objID.new)
            obj.rem <- setdiff(obj.dups$objID, obj.keep)

            obj.info <- obj.info[!(objID %in% obj.rem)]
            set.obj <- set.obj[!(objID %in% obj.rem)]

            rID <- cbind(file = "obj.info",
                         obj.dups[, .(ID.new = objID.new, ID.orig = objID)])
            rID <- data.table::data.table(file = "obj.info", setID = NA,
                                          objID = obj.rem, flag = "coord_dup")
            removed.ID <- rbind(removed.ID, rID)
            .vrb(paste0(length(obj.keep), " unique genes remain after merging\n",
                        "[info available in merged.ID]\n"))
         } else {
            .vrb("rem.genes = FALSE, identical genes will be retained\n")
         }
      }
   }

   #---------------------------------------------------------------------------#
   ## Optional operations: data filtering
   #---------------------------------------------------------------------------#

   #NOTE: removal of sets will not remove orphan genes from gene list
   # since these can still be used for permutations
   # however, gene removal can result in removal of sets
   # if all genes are removed, or when gene sets fall below a specific size
   km.ch <- list(obj.in, obj.out, set.in, set.out)
   w.ch <- !sapply(km.ch, is.null)
   if (any(w.ch)) {
      .vrb(cli::style_italic("\nFiltering files:\n"))
      W.ch <- which(w.ch)
      nor0 <- nrow(obj.info)
      sor0 <- nrow(set.info)
      nsr.d <- nor.d <- 0
      if (1 %in% W.ch) { # genes to keep
         o.in <- obj.info$objID %in% obj.in
         rID <- data.table::data.table(file = "obj.info", setID = NA,
                                       objID = obj.info[!o.in]$objID,
                                       flag = "filter")
         removed.ID <- rbind(removed.ID, rID)
         obj.info <- obj.info[o.in]
         nor <- nor0 - nrow(obj.info)
         .vrb(paste0("obj.in option enacted: ", nor, " gene",
                     ifelse(nor == 1, "", "s"), " removed\n"))
         nor.d <- length(obj.in) - nrow(obj.info)
      }
      if (2 %in% W.ch) { # genes to remove
         o.out <- obj.info$objID %in% obj.out
         rID <- data.table::data.table(file = "obj.info", setID = NA,
                                       objID = obj.info[o.out]$objID,
                                       flag = "filter")
         removed.ID <- rbind(removed.ID, rID)
         obj.info <- obj.info[!o.out]
         nor <- nor0 - nrow(obj.info)
         .vrb(paste0("obj.out option enacted: ", nor, " gene",
                     ifelse(nor == 1, "", "s"), " removed\n"))
         nor.d <- length(obj.out) - nor
      }
      if (nor.d > 0) {
         .vrb(paste0("   [NB: ", nor.d, " gene ID",
                     ifelse(nor.d == 1, " was", "s were"),
                     " not detected in obj.info]\n"))
      }
      if (3 %in% W.ch) {
         s.in <- set.info$setID %in% set.in
         rID <- data.table::data.table(file = "set.info",
                                       setID = set.info[!s.in]$setID,
                                       objID = NA, flag = "filter")
         removed.ID <- rbind(removed.ID, rID)
         set.info <- set.info[s.in]
         nsr <- sor0 - nrow(set.info)
         .vrb(paste0("set.in option enacted: ", nsr, " gene set ID",
                     ifelse(nsr == 1, "", "s"), " removed\n"))
         nsr.d <- length(set.in) - nrow(set.info)
      }
      if (4 %in% W.ch) {
         s.out <- set.info$setID %in% set.out
         rID <- data.table::data.table(file = "set.info",
                                       setID = set.info[s.out]$setID,
                                       objID = NA, flag = "filter")
         removed.ID <- rbind(removed.ID, rID)
         set.info <- set.info[!s.out]
         nsr <- sor0 - nrow(set.info)
         .vrb(paste0("set.out option enacted: ", nsr, " gene set",
                     ifelse(nsr == 1, "", "s"), " removed\n"))
         nsr.d <- length(set.out) - nsr
      }
      if (nsr.d > 0) {
         .vrb(paste0("   [NB: ", nsr.d, " gene set ID",
                     ifelse(nsr.d == 1, " was", "s were"),
                     " not detected in set.info]\n"))
      }


      # ensure sufficient genes remain after filtering
      if (nrow(obj.info) < min.set.n) {
         stop("Fewer than min.set.n genes [<",  min.set.n,
              "] remaining after filtering", call. = FALSE)
      }

      # ensure sufficient gene sets remain after filtering
      if (nrow(set.info) == 0) {
         stop("No gene sets remaining after filtering", call. = FALSE)
      }

      # update set.obj
      set.obj <- set.obj[setID %in% set.info$setID][objID %in% obj.info$objID]
   }

   #---------------------------------------------------------------------------#
   ## Optional operations: gene score generation and coordinate conversion
   #---------------------------------------------------------------------------#

   # get chromosome lengths (also used in final section)
   chr.len <- obj.info[, max(endpos) - min(startpos), by = chr]$V1
   if (bp.to.cM | get.objStat) {
      .vrb(cli::style_italic(paste("\nAuxillary files detected:",
                                   "updating obj.info file\n")))

      # generate gene scores
      if (get.objStat) {
         .vrb(paste0("Generating gene scores using residuals from ",
                     ifelse(obj.stat.fun == "non.param", "non-", ""),
                     "parametric regression: `max value` ~ `",
                     ifelse(obj.stat.fun == "non.param", "", "log "),
                     "|overlap|`\n"))

         # generate maximum scores in each gene interval
         os0 <- .get_obj_stat(OI = obj.info, VI = var.info,
                              FUN = obj.stat.fun, binN = bin.size)
         obj.info[os0[[1]], on = .(objID), objStat := objMax.res]
         obj.stat.param <- os0[[2]]
         rm(os0); gc(verbose = FALSE)

         # check and remove any genes without scores
         nst <- stats::complete.cases(obj.info)
         if (!all(nst)) {
            rw <- paste("Detected and removed", sum(!nst),
                        ifelse(sum(!nst) == 1, "gene", "genes"),
                        "lacking overlapping regions in var.info file")
            rID <- data.table::data.table(file = "obj.info", setID = NA,
                                          objID = obj.info[!nst]$objID,
                                          flag = "no_score")
            removed.ID <- rbind(removed.ID, rID)
            obj.info <- obj.info[nst]
         }
      } else {
         nst <- TRUE
         obj.stat.param <- NULL
      }

      # convert genome coordinates from physical (bp) to genetic (cM)
      if (bp.to.cM) {
         .vrb("Coverting coordinates from physical to genetic system")
         gD <- c("cM", "cm", "CM", "map", "Map", "MAP") %in% colnames(rec.rate)
         if (any(gD)) { # user has provided relevant genetic distance column
            map.fun <- NULL
            .vrb("\n")
         } else {
            .vrb(paste0(" [using the ",
                       toupper(substring(map.fun, 1, 1)), substring(map.fun, 2),
                       " map function]\n"))
         }

         if (max(chr.len) < 1000) {
            warning("Short chromosome lengths detected;",
                    "coordinates may already be genetic\n", call. = FALSE)
         }

         # reset midpos position to cM for autocovariance calculations
         obj.info[, midpos.bp := midpos]
         .bp2cM(OI = obj.info, pos = "midpos", RR = rec.rate, map.fun = map.fun) # update midpos to genetic coordinates

         # check and remove any genes lying outside rec.rate boundaries
         nrr <- stats::complete.cases(obj.info)
         if (!all(nrr)) {
            rw <- c(rw,
                    paste("Detected and removed", sum(!nrr),
                          ifelse(sum(!nrr) == 1, "gene", "genes"),
                          "lacking new coordinates in obj.info"))
            rID <- data.table::data.table(file = "obj.info", setID = NA,
                                          objID = obj.info[!nrr]$objID,
                                          flag = "no_cM")
            removed.ID <- rbind(removed.ID, rID)
            obj.info <- obj.info[nrr]
         }
      } else {
         nrr <- TRUE
      }

      if (!all(nst) | !all(nrr)) { # check revised input files
         # update set.obj
         xov <- set.obj$objID %in% obj.info$objID
         if (!all(xov)) {
            rID <- data.table::data.table(file = "set.obj", setID = NA,
                                          objID = unique(set.obj[!xov]$objID),
                                          flag = "orphan")
            removed.ID <- rbind(removed.ID, rID)
            set.obj <- set.obj[xov]
         }

         # update set.info
         xsv <- set.info$setID %in% set.obj$setID
         if (!all(xsv)) {
            rID <- data.table::data.table(file = "set.info",
                                          setID = set.info[!xsv]$setID,
                                          objID = NA, flag = "orphan")
            removed.ID <- rbind(removed.ID, rID)
            set.info <- set.info[setID %in% set.obj$setID]
         }

         removed.ID[order(file, setID, objID)]
         .vrb(paste(rw, collapse = "\n"))
         .vrb("\n[info available in removed.ID]\n")
      }
   } else {
      obj.stat.param <- NULL
   }

   #---------------------------------------------------------------------------#
   ## Ensure gene sets meet size and gene sharing specifications
   #---------------------------------------------------------------------------#

   .vrb(cli::style_italic("\nEnsuring gene sets are in prescribed range\n"))
   n.sets.orig <- nrow(set.info)
   merged.ID <- NULL
   if (n.sets.orig > 1) { # run merging
      if (set.merge == 1) {
         .vrb("Identifying and merging and identical gene sets: ")
      } else {
         .vrb(paste0("Identifying and merging gene sets with >= ",
                     round(set.merge * 100, 1), "% similarity: "))
      }

      # run merging function
      .merge_sim_sets(SI = data.table::copy(set.info),
                      SO = data.table::copy(set.obj),
                      min.sim = set.merge, ENV = environment())

      if (!is.null(set.info.merged)) {
         mID <- cbind(file = "set.info",
                      set.info.merged[, .(ID.new = setID.new, ID.orig = setID)])
         merged.ID <- rbind(merged.ID, mID)
      }

      if (R == 0) {
         .vrb(paste("none detected\n"))
      } else {
         .vrb(paste0(nrow(set.info.merged), " sets merged into ",
                     data.table::uniqueN(set.info.merged$setID.new),
                     " set", ifelse(n.sets.orig - nrow(set.info) == 1, "", "s"),
                     "\n[info available in merged.ID]\n"))
      }
   } else {
      set.info[, setN := nrow(set.obj)]
   }

   # check set size filtering criteria
   maxN <- max(set.info$setN)
   if (min.set.n > maxN) {
      stop("min.set.n [", min.set.n, "] exceeds no. genes in largest set [",
           maxN, "]\nCheck inputs and argument settings", call. = FALSE)
   } else {
      # check if sets are outside of chosen range
      set.size.rng <- min.set.n > min(set.info$setN) | !is.infinite(max.set.n)
      if (set.size.rng) { # remove gene sets that have too many / too few genes
         .vrb(paste0("Identifying and removing gene sets with < ", min.set.n,
                     ifelse(max.set.n == Inf, "", paste(" or >", max.set.n)),
                     " genes: "))
         s.rng <- set.info[, setN < min.set.n | setN > max.set.n]
         if (all(s.rng)) {
            stop("All gene sets removed\nCheck filtering criteria",
                 call. = FALSE)
         } else {
            if (sum(s.rng) == 0)  { # update objects
               .vrb("no gene sets fall outside prescribed range\n")
            } else {
               rID <- data.table::data.table(file = "set.info",
                                             setID = set.info[s.rng]$setID,
                                             objID = NA, flag = "set_size")
               removed.ID <- rbind(removed.ID, rID)
               set.info <- set.info[!s.rng]
               set.obj <- set.obj[setID %in% set.info$setID] # update set.obj
               .vrb(paste0(sum(s.rng), " set",
                           ifelse(sum(s.rng) == 1, "", "s"),
                           " removed\n[info available in removed.ID]\n"))
            }
         }
      } else {
         .vrb(paste0("All gene sets within specified size range [>= ",
                     min.set.n, ifelse(max.set.n == Inf, "",
                                       paste(" and <=", max.set.n)),
                     " genes]\n"))
      }
   }

   #---------------------------------------------------------------------------#
   ## Prepare input for polylinkR
   #---------------------------------------------------------------------------#

   # determine if there is gene sharing across final configuration of gene sets
   no.share <- all(set.obj[, data.table::uniqueN(setID), by = objID]$V1 == 1)

   # determine genetic coordinates
   coord <- ifelse(bp.to.cM, "cM", "bp")

   # retain auxiliary data, function arguments, and summaries
   n.genes <- obj.info[, .N]
   n.sets <- set.info[, .N]
   n.set.genes <- data.table::uniqueN(set.obj$objID)
   if (pos.info) {
      n.chr <- data.table::uniqueN(obj.info$chr)
   } else {
      n.chr <- NULL
   }
   plr.data <- list()
   read.data <- list(n.genes = n.genes, n.sets = n.sets, n.chr = n.chr,
                     n.cov = n.cov, n.set.genes = n.set.genes,
                     file.paths = file.paths, coord = coord,
                     get.objStat = get.objStat, cov.names = cov.names,
                     cov.info = cov.info, no.share = no.share,
                     pos.info = pos.info)
   plr.data$read.data <- read.data

   plr.args <- list()
   args0 <- names(formals(get("plR_read", envir = asNamespace("polylinkR"))))
   read.args <-  mget(args0, envir = environment()) # collect arguments
   plr.args$read.args <- read.args

   plr.summary <- list()
   read.summary <- list(removed.ID = removed.ID, repeat.ID = repeat.ID,
                        merged.ID = merged.ID, obj.stat.param = obj.stat.param)
   plr.summary$read.summary <- read.summary

   # capture R session information and record run time
   plr.session <- list()
   r0 <- .get_time(start.time = startT)
   plr.session$read.session <- list(session = sessionInfo(), run.time = r0[[2]])

   # retain core files
   data.table::setkey(obj.info, key = NULL)
   data.table::setkey(set.info, key = NULL)
   data.table::setkey(set.obj, key = NULL)
   OUT <- list(set.info = set.info, obj.info = obj.info, set.obj = set.obj)

   # set s3 class and create attributes
   plr.out <- .new_plR(BASE = OUT, plR.data = plr.data, plR.args = plr.args,
                       plR.summary = plr.summary, plR.session = plr.session)

   .vrb(cli::style_italic("\nFile loading summary:\n"))
   .vrb(paste0("Loaded ", nrow(set.info), " unique gene set",
               ifelse(n.sets == 1, "", "s"), " and ", nrow(obj.info),
               " unique genes\n"))

   obj.nx <- nrow(obj.info) - data.table::uniqueN(set.obj$objID)
   if(obj.nx > 0) {
      .vrb(paste0("[NB: ", obj.nx, " genes (",
                  round(obj.nx / nrow(obj.info) * 100, 1),
                  "% of total genes) are not in any gene set ",
                  "but are retained for computing null and autocovariance]\n"))
   }

   # report on data required specific processes
   if (cov.info) {
      c.mss <- paste("Covariate columns detected in obj.info:",
                     "plR_permute will initially run confounder correction and",
                     "use deconfounded gene scores in gene set enrichment")
   } else {
      c.mss <- paste("Covariate columns not detected in obj.info:",
                     "plR_permute will not run confounder correction and",
                     "use raw gene scores in gene set enrichment")
   }

   if (pos.info) {
      p.mss <- paste("Genomic coordinate columns detected in obj.info:",
                     "plR_rescale will compute intergene autocovariance and",
                     "recalculate enrichment for rescaled gene set scores")
   } else {
      p.mss <- paste("Genomic coordinate columns not detected in obj.info:",
                     "plR_rescale requires user-provided autocovariance",
                     "object (ac) to revise gene set scores and enrichment")
   }

   if (no.share) {
      if (n.sets == 1) {
         q.mss <- paste("Only one gene set in set.info:",
                        "pruning using plR_prune not required")
      } else {
         q.mss <- paste("No shared genes detected between gene sets in set.obj:",
                        "No pruning required; standard q values will be",
                        "reported when using plR_prune")
      }
   } else {
      q.mss <- paste("Shared genes detected between gene sets in set.obj:",
                     "plR_prune will estimate q values for gene set scores",
                     "after correcting for shared genes")
   }
   .vrb(paste(c(paste("**", c.mss, "**"), paste("**", p.mss, "**"),
                paste("**", q.mss, "**")), collapse = "\n"))
   .vrb("\n\n")

   ftr <- paste0("Finished loading files -- run time: ", r0[[1]])
   pdg <- (80 - nchar(ftr)) / 2
   .vrb(cli::boxx(ftr, padding = c(0, floor(pdg), 0, ceiling(pdg)),
                  border_style = "double", col = "cyan", border_col = "cyan"))
   .vrb("\n")
   return(plr.out)
}


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
#' @param plR.input
#'   \code{plR} class object; output from \code{polylinkR::plR_read}.
#'   Required.
#'
#' @param permute
#'   \code{logical}; should the permutation step be performed? Defaults to
#'   \code{TRUE}. If \code{FALSE}, only gene score deconfounding is performed
#'   without estimating set scores. If \code{TRUE}, users must provide
#'   confounder covariates in \code{obj.info}. This is useful for exploring
#'   how parameter settings impact deconfounded gene scores; the permutation
#'   step can be performed later by passing the output back into
#'   \code{plR_permute}.
#'
#' @param n.perm
#'   \code{integer}; number of permutations. Defaults to \code{1e5}. Must be
#'   in the range \code{[5e4L, Inf)} and be exactly divisible by \code{1e4}.
#'
#' @param n.boot
#'   \code{integer}; number of bootstrap replicates for null inference.
#'   Defaults to \code{30L}. Must be in the range \code{[5, Inf)}.
#'
#' @param alt
#'   \code{character}; direction of the hypothesis test. Defaults to
#'   \code{"upper"} for enrichment in the upper tail (large set scores).
#'   Alternatively, \code{"lower"} tests for enrichment in the lower tail
#'   (small values). When \code{"lower"} is chosen, data are internally
#'   negated in functions performing p-value estimation.
#'
#' @param md.meth
#'   \code{character}; determines whether raw covariate data or ranks are
#'   used in Mahalanobis distance calculations. Defaults to \code{"robust"},
#'   where Mahalanobis distances are converted to ranks and Spearman's metric
#'   is used to calculate the covariance matrix. The alternate option
#'   \code{"raw"} uses the original covariates and Pearson's covariance for
#'   scaling.
#'
#' @param kern.bound
#'   \code{numeric}; flanking region around each gene where weights of
#'   overlapping genes inside the region are set to 0. Weights of partially
#'   overlapping genes are downscaled by the proportion overlapping the
#'   excluded region. Defaults to 0.1 Mbp or 0.1 cM (depending on genetic
#'   coordinates). Range \code{(0, Inf]}; \code{0} denotes standard gene
#'   boundaries and \code{Inf} excludes all genes on the same chromosome.
#'   Ignored if appropriate covariate columns are not detected in
#'   \code{obj.info}.
#'
#' @param kern.func
#'   \code{character}; kernel function used to generate probability weights
#'   from distances between the focal gene and other genes in confounder
#'   space. Default is \code{"normal"} (Gaussian kernel). Alternate options
#'   include \code{"exponential"} and \code{"inverse"}. Ignored if covariate
#'   columns are not detected in \code{obj.info}.
#'
#' @param kern.scale
#'   \code{numeric}; scalar used in the kernel function to convert
#'   Mahalanobis distances to probabilities. Defaults to \code{2} for
#'   Gaussian decay, \code{log(10)} for exponential decay, or \code{2} for
#'   inverse decay. Range \code{(0, Inf]}. Ignored if covariate columns are
#'   not detected in \code{obj.info}.
#'
#' @param kern.wt.max
#'   \code{numeric}; maximum probability weight for a single gene. Defaults
#'   to \code{0.05}. Must be in the range \code{(1 / (n.genes - 1), 1]}.
#'   Set to \code{1} if no upper bound is desired. Ignored if covariate
#'   columns are not detected in \code{obj.info}.
#'
#' @param gpd.cutoff
#'   \code{numeric}; threshold tail probability at which to apply GPD tail
#'   fitting. Defaults to \code{500 / n.perm}. Must be in the range
#'   \code{[max(c(1e-04, 500 / n.perm)), 0.05]}. The lower bound constraint
#'   ensures that a minimum of \code{500} exceedances are available for GPD
#'   estimation, while also ensuring compatibility with the empirical CDF, where
#'   the lowest evaluated quantile is \code{1e-4}.
#'
#' @param seed
#'   \code{integer}; random seed for reproducibility. Preserved across
#'   subsequent polylinkR functions (\code{plR_permute}, \code{plR_rescale}).
#'   Defaults to \code{NULL}, in which case a seed is generated
#'   automatically. Must be within
#'   \code{[-.Machine$integer.max, .Machine$integer.max]}.
#'
#' @param verbose
#'   \code{logical}; should progress messages be printed to the console?
#'   Defaults to \code{TRUE}.
#'
#' @param n.cores
#'   \code{integer}; number of cores for parallel processing. Defaults to
#'   \code{1} or \code{maximum cores - 1}. Must be in the range
#'   \code{[1, maximum cores]}.
#'
#' @param fut.plan
#'   \code{character}; parallel backend from the \code{future} package.
#'   Defaults to user \code{n.cores} choice or checks available cores,
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
#' # Assuming `my_plr` is the result of `polylinkR::plR_read`
#'
#' # Example 1: Basic usage
#' new_plr <- plR_permute(plR.input = my_plr)
#'
#' # Example 2: Less permutations, more bootstraps
#' new_plr <- plR_permute(
#'    plR.input = my_plr,
#'    n.perm    = 1e5,
#'    n.boot    = 100
#')
#'
#' # Example 3: Modified covariate handling, single processor
#' new_plr <- plR_permute(
#'   plR.input = my_plr,
#'   kern.wt.max  = 0.2,
#'   md.meth   = "raw",
#'   n.cores   = 1
#' )
#'
#' # Example 4: Modified GPD estimation, user-specified seed
#' new_plr <- plR_permute(
#'   plR.input  = my_plr,
#'   gpd.cutoff = 0.01,
#'   seed       = 1000
#' )
#'
#' # Example 5: Only deconfound scores (no enrichment analysis)
#' new_plr <- plR_permute(
#'    plR.input = my_plr,
#'    permute =   FALSE
#' )
#'
#' # Example 6: Use deconfounded scores from Example 5 to run enrichment
#' new_plr <- plR_permute(plR.input = new_plr)
#' }
plR_permute <- function(plR.input, permute = TRUE, n.perm = 5e5L, n.boot = 30L,
                        alt = "upper", md.meth = "robust", kern.wt.max = 0.05,
                        kern.bound = "auto", kern.func = "normal",
                        kern.scale = "auto", gpd.cutoff = 5e3L / n.perm,
                        seed = NULL, verbose = TRUE, n.cores = "auto",
                        fut.plan = "auto") {

   ##=========================================================================##
   ## PART 1: Clean data and run checks
   ##=========================================================================##

   # track function run time
   startT <- Sys.time()

   # control message reporting
   rlang::local_options(verbose = verbose, .frame = environment())

   # perform checks and unpack required plR objects
   .plR_check(f = "permute", ENV = environment())
   .arg_check(f = "permute", ENV = environment())
   .plR_unpack(plr = plR.input, ENV = environment())
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

   hdr[1] <- paste("Running plR_permute:", hdr[1])
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
   ## PART 2: Calculate standardised gene scores
   ##=========================================================================##

   if (perm.path == "full") { # perform local regression
      .vrb(cli::style_italic(paste("Estimating prognostic gene scores",
                                   "(confounder effect) using local quadratic",
                                   "regression:\n")))
      .vrb("Calculating pairwise Mahalanobis Distances for all genes")

      cv.val <- as.matrix(obj.info[, .SD, .SDcols = cov.names]) # covariate matrix

      if (md.meth == "robust") {
         .vrb(" in robust covariate space (coverting to normal scores)\n")
         cv.val <- Rfast::colRanks(cv.val, method = "average") # generate ranks
         cv.val <- qnorm((cv.val - 0.5) / n.genes) # convert ranks to normal scores
      } else {
         .vrb(" in covariate space\n")
      }

      # check for singularity among covariates
      if(ncol(cv.val) > 1) {
         qrX <- qr(cv.val)
         if (qrX$rank < length(cov.names)) {
            col.keep <- qrX$pivot[1:qrX$rank]
            cv.val <- cv.val[, col.keep] # remove redundant covariates
            .vrb(paste0("Covariate matrix is singular\nRedundant covariates [",
                        paste0("Cov", setdiff(1:ncol(cov.mat), col.keep),
                               collapse = ", "),
                        "] are ommited from all regression models\n"))
         } else {
            col.keep <- 1:length(cov.names) # no redundancy among covariates
         }
      }

      # generate Mahalanobis distances
      M0 <- distances::distances(data = cv.val, normalize = "mahalanobize")

      .vrb(paste("Estimating regression coefficients for each gene using",
                 kern.func, "kernel weights \n"))
      if (kern.scale == "auto") {
         wK <- which(c("normal", "exponential", "inverse") == kern.func)
         kern.scale <- c(2, log(10), 2)[wK]
      }
      .vrb("\n")

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
                 .fit_lqr(d1 = d.i, wt0 = wt.i[i - min(I) + 1, ], os0 = os0,
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
         .vrb("Local quadratic regression failed for some genes:\n")
         if (any(lqr.fit == n.cov + 1)) {
            lf <- tfit[names(tfit) == n.par[2]]
            .vrb(paste(lf, "gene", ifelse(lf == 1, "", "s"), "had linear fit"))
         }
         .vrb(ifelse(length(tfit) == 3, " and ", " "))
         if (any(lqr.fit == 1)) {
            nf <- tfit[names(tfit) == n.par[3]]
            .vrb(paste(nf, "gene", ifelse(lf == 1, "", "s"), "had 0-order fit"))
         }
         .vrb("\n[check lqr.fit object in attributes]\n")
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

      .vrb("Generating deconfounded gene scores and model fitting statistics\n")
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
         .vrb(paste("No covariate columns identified:",
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
      .vrb(cli::style_italic("\nGenerating null set score distributions:\n"))
      .vrb(paste0("Estimating empirical CDFs and fitting GPDs using ",
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

      fopt <- list(packages = c("data.table", "foreach", "dqrng", "tdigest",
                                "Rfast"),
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

                  e0 <- lapply(csX, FUN = function(csI) { # extract digests
                     as.list(tdigest::tdigest(csI, compression = 500))
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
         .vrb(paste("Smoothing empirical CDF and GPD parameter estimates",
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
            eSmFit <- .par_smooth(y = ei, x = th0, K = K, n.th = n.th,
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
            gSmFit <- .par_smooth(y = gM$V1, x = th0, K = K, n.th = n.th,
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
         .vrb(paste("Identifying and setting minimum possible p values\n"))
         fopt <- list(packages = "data.table",
                      globals = c("gpd.cutoff", "q.bnd", ".get_gpd_mins"),
                      seed = seed)
         gpdM <- foreach::foreach(gpdI = split(gpd0[-1], 2:n.max),
                                  .options.future = fopt,
                                  .combine = rbind) %dofuture% {
            .get_gpd_mins(gpdI, gpd.cutoff = gpd.cutoff, q.bnd = q.bnd)
         }
         gpd0[, adj.p := c(NA, gpdM[, "adj.p"])]
         gpd0[, min.z := c(NA, gpdM[, "min.z"])]
         gpd0[, min.p := c(NA, gpdM[, "min.p"])]

         rm(g0, gpdM); gc(verbose = FALSE) # clean up
      }
      .vrb("\n")
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
      .vrb(cli::style_italic("Calculating test statistics:\n"))
      .vrb(paste0("Calculating observed and null gene set scores\n"))

      # calculate observed standardised set scores
      # mean of standardised gene scores = 0 and variance = 1 -> set score variance = set size (assuming independence)
      ss.obs <- set.obj[, sum(os0[objID]) / sqrt(.N), keyby = setID]
      data.table::setkey(set.info, setID)
      set.info[.(ss.obs), setScore.std := V1]

      # calculate p values
      .vrb("Calculating p values for each gene set\n")
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

   if (("permute.summary" %in% names(plr.summary))) { # update ecdf and gpd estimation summaries
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
   plr.out <- .new_plR(BASE = OUT, plR.data = plr.data, plR.args = plr.args,
                       plR.summary = plr.summary, plR.seed = plr.seed,
                       plR.session = plr.session)

   .vrb("\n")
   ftr <- cli::col_cyan(paste0("Finished plR_rescale -- run time: ", r0[[1]]))
   pdg <- (80 - nchar(ftr)) / 2
   .vrb(cli::boxx(ftr, padding = c(0, floor(pdg), 0, ceiling(pdg)),
                  border_style = "double", , col = "cyan", border_col = "cyan"))
   .vrb("\n")
   return(invisible(plr.out))
}


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
#' @param plR.input
#'   \code{plR} class object; output from \code{polylinkR::plR_permute}.
#'   Required.
#'
#' @param rescale
#'   \code{logical}; should gene set score (i.e.,\code{setScore.std}) rescaling
#'   be performed? Defaults to \code{TRUE}. If \code{FALSE}, only inter-gene
#'   autocorrelation is estimated without revising enrichment testing for
#'   rescaled (i.e., decorrelated) gene set scores. This is useful for exploring
#'   how parameter settings impact gene set autocorrelation. The rescaling step
#'   can be performed later by passing the output to \code{plR_rescale}, or by
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
#' @param cgm.range
#'   \code{numeric}; maximum inter-gene lag used to evaluate autocovariance.
#'   Defaults to \code{2e6} bp or \code{2} cM, depending on genetic distance
#'   measure. Must be within interval \code{[1e5, 5e7]} bp or \code{[0.1, 50]}
#'   cM.
#'
#' @param cgm.bin
#'   \code{numeric}; mimimum number of gene pairs required for bins of empirical
#'   covariances. Default value = \code{30}. Must be in range \code{[10, 1e3]}.
#'   Initially an exponential grid of bin sizes will be generated, favouring
#'   smaller bins at short distances and larger bins for more distant gene
#'   pairs. Bins with too few genes will be successively merged with the larger
#'   bins until the minimum number of genes are reached. This condition is
#'   evaluated across all chromosomes, ensuring no bin is smaller than the
#'   minimum value (other than the final bins).
#'
#' @param cgm.wt.max
#'   \code{numeric}; maximum probability weight for a single lag window.
#'   Defaults to \code{0.1}. Set to \code{1} if no upper bound is desired.
#'   Note that the lower bound is reset to \code{2 / no. fitted lags} if the
#'   chosen value falls below this limit.
#'
#' @param emp.bayes
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
#' @param min.rho
#'   \code{numeric}; minimum estimated correlation between two gene sets;
#'   values below this are set to \code{0}. Defaults to \code{1e-5}. Range
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
#' # Assuming `my_plr` is the result of `polylinkR::plR_permute`
#'
#' # Example 1: Basic usage
#' new_plr <- plR_rescale(plR.input = my_plr)
#'
#' # Example 2: Estimate autocorrelation only (no rescaling)
#' new_plr <- plR_rescale(
#'    plR.input = my_plr,
#'    rescale = FALSE
#' )
#'
#' # Example 3: Custom variogram model and cores
#' new_plr <- plR_rescale(
#'   plR.input = my_plr,
#'   vg.model  = c("Exp","Sph"),
#'   n.cores   = 4
#' )
#'
#' # Example 4: Using a user-provided autocovariance object
#' # Assuming `my_ac` is a valid data.table or data.frame
#' new_plr <- plR_rescale(
#'    plR.input = my_plr,
#'    ac        = my_ac
#' )
#'
#' # Or using the `new_plr` object generated by Example 2
#' new_plr <- plR_rescale(plR.input = my_rescaled_plr)
#' }
plR_rescale <- function(plR.input, rescale = TRUE, fast = TRUE, ac = NULL,
                        cgm.bin = 30, cgm.range = "auto", cgm.wt.max = 0.05,
                        emp.bayes = "auto", min.rho = 1e-5, verbose = TRUE,
                        n.cores = "auto", fut.plan = "auto") {

   ##==========================================================================##
   ## PART 1: Clean data and run checks
   ##==========================================================================##

   # track function run time
   startT <- Sys.time()

   # control message reporting
   rlang::local_options(verbose = verbose, .frame = environment())

   # perform checks and unpack required plR objects
   .plR_check(f = "rescale", ENV = environment())
   .arg_check(f = "rescale", ENV = environment())
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

   hdr[1] <- paste("Running plR_rescale:", hdr[1])
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

   ##==========================================================================##
   ## PART 2: Calculate genomic autocovariance
   ##==========================================================================##

   if (is.null(ac)) { # infer genetic autocovariance
      .vrb(cli::style_italic("Estimating genetic autocovariance:\n"))

      # generate exponential lag bin sizes
      if (cgm.range == "auto") { # determine maximum correlogram range if not provided
         maxL <- ifelse(coord == "cM", 3L, 3e6L)
      }

      .vrb(paste0("Creating correlogram lag bins [max ", maxL, coord,
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

      .vrb("Assigning emprical covariances to lag bins\n")
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

      .vrb(paste("Fitting autocovariance function\n"))

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
         .vrb("Apply empirical Bayes to revise chromosome parameters\n")
         ebPar0 <- .emp_bayes_est(cPar = cPar, gPar = gPar, cN0 = cN,
                                  corr = FALSE)
         if (emp.bayes == "auto") {
            emp.bayes <- ifelse(n.chr >= 15, "full", "reduced")
         }

         if (emp.bayes == "full") {
            .vrb("Comparing full and marginal models\n")
            ebPar <- .emp_bayes_est(cPar = cPar, gPar = gPar, cN0 = cN,
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
      .vrb("Generating intergene covariances\n")

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
         xIJ[, rho := .cov_fun(h = pos.lag, ovlp = cov.ovlp, scale = par$scale,
                               beta1 = par$beta.LD.reg,
                               beta2 = par$beta.ovlp.reg)]
         cbind(chr = x$chr[1],  xIJ[rho >= min.rho])
      }
      ac <- do.call(rbind, ac)

      user.ac <- FALSE
   } else { # user has provided valid ac object
      .vrb(paste("Valid ac object detected:",
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
      .vrb(cli::style_italic("\nRescaling null set score distributions:\n"))

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
         .vrb(paste("fast = TRUE: Rescaling empirical CDFs and GPD scale",
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

         .vrb(paste0("fast = FALSE: Set-wise rescaling of empirical CDFs and",
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

         fopt <- list(packages = c("data.table", "foreach", "dqrng", "tdigest",
                                   "Rfast"),
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
                     rsX <- .get_cov0(csj = csX, rsX = rs0, SS0 = SS0, ac0 = ac0) # get covariances
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

                     e0 <- lapply(csX, FUN = function(csI) { # extract digests
                        as.list(tdigest::tdigest(csI, compression = 500))
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
            .vrb(paste("Smoothing empirical CDF and GPD parameter estimates",
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
               eSmFit <- .par_smooth(y = ei, x = th0, K = K, n.th = n.th,
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
               gSmFit <- .par_smooth(y = gM$V1, x = th0, K = K, n.th = n.th,
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
            .vrb(paste("Identifying and setting minimum possible p values\n"))
            fopt <- list(packages = "data.table",
                         globals = c("gpd.cutoff", "q.bnd", ".get_gpd_mins"),
                         seed = seed)

            gpdM <- foreach::foreach(gpdI = split(gpd0[-1], 2:n.max),
                                     .options.future = fopt,
                                     .combine = rbind) %dofuture% {
            .get_gpd_mins(gpdI, gpd.cutoff = gpd.cutoff, q.bnd = q.bnd)
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
         .vrb("\n")
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
      .vrb(cli::style_italic("\nCalculating test statistics:\n"))

      # compute rescaling factor for observed set scores
      SSI <- set.obj[, .(oID = objID), keyby = setID][, .(sID = setID, oID)]
      ac0 <- ac[, .(A = objID.A, B = objID.B, CV = 2 * rho)]

      obs.cov <- rep(0, n.sets)
      oc0 <- .get_cov(SSi = SSI, ac0 = ac0)
      obs.cov[oc0$sID] <- oc0$CV

      # calculate observed rescaled set scores
      set.info[, setScore.rs := setScore.std * sqrt(setN / (setN + obs.cov))]

      # calculate p values
      .vrb("Calculating p values for each gene set\n")
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
   plr.out <- .new_plR(BASE = OUT, plR.data = plr.data, plR.args = plr.args,
                       plR.summary = plr.summary, plR.seed = plr.seed,
                       plR.session = plr.session)

   .vrb("\n")
   ftr <- cli::col_cyan(paste0("Finished plR_rescale -- run time: ", r0[[1]]))
   pdg <- (80 - nchar(ftr)) / 2
   .vrb(cli::boxx(ftr, padding = c(0, floor(pdg), 0, ceiling(pdg)),
                  border_style = "double", col = "cyan", border_col = "cyan"))

   return(invisible(plr.out))
}


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

