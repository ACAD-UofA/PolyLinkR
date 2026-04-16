#' @title Read and validate required files for polylinkR
#'
#' @description
#' Reads and validates the input files for a \code{polylinkR} analysis. Handles
#' file path resolution, column name canonicalisation, data validation, and
#' optional filtering and data generation tasks like oordinate conversion or
#' gene score calculation.
#'
#' @param input_path
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
#' @param object_info_path
#'   \code{character}; path to the object info file. Defaults to
#'   \code{NULL}. Ignored if \code{input_path} is provided; otherwise all
#'   required file paths must be specified.
#'
#' @param gene_set_info_path
#'   \code{character}; path to the gene set info file. Defaults to
#'   \code{NULL}. Ignored if \code{input_path} is provided; otherwise all
#'   required file paths must be specified.
#'
#' @param gene_set_mapping_path
#'   \code{character}; path to the gene set mapping file. Defaults to
#'   \code{NULL}. Ignored if \code{input_path} is provided; otherwise all
#'   required file paths must be specified.
#'
#' @param variant_info_path
#'   \code{character}; path to the variant info file. Defaults to
#'   \code{NULL}. Optional file used for gene score generation. Ignored if
#'   \code{input_path} is provided.
#'
#' @param recombination_rate_path
#'   \code{character}; path to the recombination rate file. Defaults to
#'   \code{NULL}. Optional file used for genetic coordinate conversion.
#'   Ignored if \code{input_path} is provided.
#'
#' @param min_set_size
#'   \code{integer}; minimum size of gene sets to be retained. Defaults to
#'   \code{2}. Must be in the range \code{[2L, max_set_size)}.
#'
#' @param max_set_size
#'   \code{integer}; maximum size of gene sets to be retained. Defaults to
#'   \code{Inf}. Must be in the range \code{(min_set_size, Inf)}.
#'
#' @param group_label
#'   \code{character}; label used to identify input files within a directory.
#'   Defaults to \code{NULL}.
#'
#' @param mapping_function
#'   \code{character}; mapping function to convert physical to genetic
#'   distances. Options are \code{"Haldane"}, \code{"Kosambi"} (default),
#'   \code{"Carter-Falconer"}, and \code{"Morgan"}.
#'
#' @param object_buffer
#'   \code{numeric}; interval around genes (in base pairs) to include when
#'   assigning values from the variant info file. Defaults to \code{1e4} if
#'   variant info is provided and user does not set a value, otherwise it is
#'   set to 0 if score assignment is not performed. User values must be in the
#'   range \code{[0, 1e5L]}. Note that if the user provides their own gene
#'   scores in the object info input file (i.e. not computed from the
#'   variant info file), then the start and end positions must include any
#'   buffer used to bin scores, otherwise polylinkR deconfounding and
#'   autocorrelation inference will not be performed appropriately.
#'
#' @param object_statistic_function
#'   \code{character}; function used to correct maximum gene scores based on
#'   the number of overlapping summary statistics (SNPs or windows). Default is
#'  \code{non.param}, a robust non-parametric method that uses binned data to
#'  calculate median and median absolute deviation (MAD) to normalise scores.
#'  Alternatively, \code{lm.logN} applies a linear regression to the log-transformed
#'  SNP / bin counts (assumes a roughly linear relationship is appropriate). In
#'  both cases, expected scores are estimated and gene scores calculated as the
#'  residual value. Ignored if no variant info file is provided.
#'
#' @param bin_size
#'   \code{integer}; gene set size interval for non-parametric correction.
#'   Defaults to \code{250L}. Must be in the range \code{[50L, 1e3L]};
#'   ignored if the parametric function is used.
#'
#' @param objects_to_include
#'   \code{character} or \code{numeric} vector; object IDs of genes to
#'   explicitly retain. Defaults to \code{NULL}.
#'
#' @param objects_to_exclude
#'   \code{character} or \code{numeric} vector; object IDs of genes to
#'   explicitly remove. Defaults to \code{NULL}.
#'
#' @param sets_to_include
#'   \code{character} or \code{numeric} vector; set IDs of gene sets
#'   to explicitly retain. Defaults to \code{NULL}.
#'
#' @param sets_to_exclude
#'   \code{character} or \code{numeric} vector; set IDs of gene sets
#'   to explicitly remove. Defaults to \code{NULL}.
#'
#' @param merge_threshold
#'   \code{numeric}; minimum proportion of shared genes for merging gene sets.
#'   Defaults to \code{0.95}. Must be in the range \code{(0, 1]}.
#'
#' @param remove_duplicate_genes
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
#' \code{print()} or by simply typing the object's name.
#'
#' @examples
#' \dontrun{
#' # Example 1: Read all files from a single folder
#' my_plr <- read_polylinkr_data(input_path = "path/to/files")
#'
#' # Example 2: Read all files with "POP1" in label from a single folder
#' my_plr <- read_polylinkr_data(
#'    input_path = "path/to/files",
#'    group_label = "POP1"
#' )
#'
#' # Example 3: Relax gene set merging criteria
#' my_plr <- read_polylinkr_data(
#'    input_path = "path/to/files",
#'    merge_threshold = 0.50
#' )
#'
#' # Example 4: Specify separate file paths and remove user-specified sets
#' my_plr <- read_polylinkr_data(
#'   gene_set_info_path = "path/to/set.info",
#'   gene_set_mapping_path = "path/to/set.obj",
#'   object_info_path = "path/to/obj.info",
#'   sets_to_exclude = c("set1", "set2")
#' )
#'
#' # Example 5: Generate gene scores from variant info using regression
#' # and include a 50 kb buffer around genes
#' my_plr <- read_polylinkr_data(
#'   gene_set_info_path = "path/to/set.info",
#'   gene_set_mapping_path = "path/to/set.obj",
#'   object_info_path = "path/to/obj.info",
#'   variant_info_path = "path/to/var.info",
#'   object_statistic_function = "lm.logN",
#'   object_buffer = 5e4
#' )
#'
#' # Example 6: Convert distances to cM using recombination rate file
#' my_plr <- read_polylinkr_data(
#'   gene_set_info_path = "path/to/set.info",
#'   gene_set_mapping_path = "path/to/set.obj",
#'   object_info_path = "path/to/obj.info",
#'   recombination_rate_path = "path/to/rec.rate"
#' )
#' }
read_polylinkr_data <- function(input_path = NULL, object_info_path = NULL,
                                 gene_set_info_path = NULL, gene_set_mapping_path = NULL,
                                 variant_info_path = NULL, recombination_rate_path = NULL, 
                                 min_set_size = 2L,
                                 max_set_size = Inf, group_label = NULL, mapping_function = "kosambi",
                                 object_buffer = "auto", object_statistic_function = "non.param",
                                 bin_size = 250L, objects_to_include = NULL, objects_to_exclude = NULL,
                                 sets_to_include = NULL, sets_to_exclude = NULL, merge_threshold = 0.95,
                                 remove_duplicate_genes = FALSE, verbose = TRUE) {
   
   # Backward compatibility: map new snake_case parameters to legacy names
   # This allows the function body to continue using legacy variable names
   input.path <- input_path
   obj.info.path <- object_info_path
   set.info.path <- gene_set_info_path
   set.obj.path <- gene_set_mapping_path
   var.info.path <- variant_info_path
   rec.rate.path <- recombination_rate_path
   min.set.n <- min_set_size
   max.set.n <- max_set_size
   group <- group_label
   map.fun <- mapping_function
   obj.buffer <- object_buffer
   obj.stat.fun <- object_statistic_function
   bin.size <- bin_size
   obj.in <- objects_to_include
   obj.out <- objects_to_exclude
   set.in <- sets_to_include
   set.out <- sets_to_exclude
   set.merge <- merge_threshold
   rem.genes <- remove_duplicate_genes
   # verbose unchanged

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


#' @title Read and validate required files for polylinkR (deprecated)
#' @description
#' This function is deprecated. Please use \code{read_polylinkr_data()} instead.
#' @param ... All arguments passed to \code{read_polylinkr_data()}.
#' @export
#' @keywords internal
plR_read <- function(...) {
   .Deprecated("read_polylinkr_data", package = "polylinkR",
               msg = "plR_read() is deprecated. Use read_polylinkr_data() instead.")
   read_polylinkr_data(...)
}
