#' @title Internal constructor for \code{plR} objects
#'
#' @description
#' Internal S3 method that creates and validates new \code{plR} objects. Checks
#' function history and creates a tracking attribute that is used in
#' validation.
#'
#' @param BASE
#'   \code{list}; core data for the \code{plR} object, i.e. \code{set.info},
#'   \code{obj.info}, and \code{set.obj}.
#'
#' @param plR.data
#'   \code{list}; ancillary data required by downstream \code{polylinkR}
#'   functions.
#'
#' @param plR.args
#'   \code{list}; arguments passed to \code{polylinkR} functions.
#'
#' @param plR.summary
#'   \code{list}; \code{polylinkR} summary data from internal model fitting used
#'   in diagnostics and plotting.
#'
#' @param plR.seed
#'   \code{list}; random seed information used internally.
#'
#' @param plR.session
#'   \code{list}; function run time and session information from
#'   \code{sessionInfo}.
#'
#' @return
#' A newly constructed \code{plR} object.
#'
#' @keywords internal
#' @noRd
.new_plr <- function(BASE = NA, plR.data = NULL, plR.args = NULL,
                     plR.summary = NULL, plR.seed = NULL, plR.session = NULL) {
   if (all(c("set.info", "obj.info", "set.obj") %in% names(BASE))) {
      plr.track <- c(0, 0, 0)
      # diagnose input file structure
      if (!is.null(plR.args$permute_args)) {
         permute <- plR.args$permute_args$permute
         no_deconf <- plR.data$permute_data$no_deconf
         plr.track[1] <- ifelse(permute, ifelse(no_deconf, 3, 2), 1)
      }

      if (!is.null(plR.args$rescale_args)) {
         rescale <- plR.args$rescale_args$rescale
         user_ac <- plR.data$rescale_data$user_ac
         plr.track[2] <- ifelse(rescale, ifelse(user_ac, 3, 2), 1)
      }

      if (!is.null(plR.args$prune_args)) {
         plr.track[3] <- 1
      }

      plr.track <- paste(plr.track, collapse = "")
   } else { # core data missing, not a proper plR object
      plr.track <- "INVALID"
   }
   structure(.Data = BASE, plR.data = plR.data, plR.args = plR.args,
             plR.summary = plR.summary, plR.seed = plR.seed,
             plR.session = plR.session, plr.track = plr.track, class = "plR")
}


#' @title Validate \code{plR} class object
#'
#' @description
#' Internal function performing validation checks on \code{plR} objects used
#' as input to core \code{polylinkR} functions. Checks depend on the target
#' function \code{f}.
#'
#' @param f
#'   \code{character}; name of the \code{polylinkR} function being validated.
#'
#' @param ENV
#'   \code{environment}; environment in which the \code{plR} object is
#'   stored.
#'
#' @import data.table
#'
#' @return
#' No return value. If validation fails, the function exits with a
#' descriptive error message. If validation succeeds, a vectorised plR track is
#' assigned to \code{ENV}.
#'
#' @keywords internal
#' @noRd
.check_plr_object <- function(f, ENV) {
   plr <- deparse(substitute(plr_input, env = ENV))
   if (plr == "") {
      stop("plr_input is empty; please provide valid plr input", call. = FALSE)
   } else {
      pT <- attributes(get("plr_input", envir = ENV))$plr.track
      if (is.null(pT)) {
         stop("plr_input = ", plr, " is not a plr class object", call. = FALSE)
      } else {
         pT.all <- .get_processing_history()
         f0 <- paste0("plr_", f)
         req.track <- unlist(strsplit(pT.all[FUNCTION == f0]$INPUT, "; "))
         if (pT %in% req.track) { # return output to parent environment
            plr.track <- pT
            assign(x = "plr.track",
                   value = as.numeric(unlist(strsplit(pT, split = ""))),
                   envir = ENV)
         } else {
            stop("plr_input = ", plr, " is not valid input for plr_",
                 f, "\nCheck header of print(", plr, ") for valid usage options",
                 call. = FALSE)
         }
      }
   }
}


#' @title Print method for \code{plR} objects
#'
#' @description
#' Prints an overview of a \code{plR} object, including its processing
#' history and available contents.
#'
#' @param x
#'   \code{plR} class object; typically output from a \code{polylinkR}
#'   function.
#'
#' @param ...
#'   Additional arguments passed to \code{print}.
#'
#' @exportS3Method print plR
#' @keywords internal
#'
#' @import data.table
#' @importFrom cli boxx col_cyan col_red
#'
#' @return
#' An overview of the \code{plR} object.
#'
#' @examples
#' \dontrun{
#' # Provide a snapshot of core files (assuming `my_plR` is a valid plR object)
#' print(my_plR)
#'
#' # Simply typing the object name also prints the snapshot
#' my_plR
#' }
print.plR <- function(x, ...) {
   pT <- attributes(x)$plR.track # plR track info
    if (pT == "INVALID") { # check if empty plR object
       cat(cli::col_red("Empty plR object\n"))
    } else {
       # unpack input file information
       list2env(attributes(x)$plR.data$read.data, envir = environment())
       pT.all <- .plR_track()
       # check information
      pI <- sapply(lapply(strsplit(pT.all$INPUT, "; "), '%in%', pT), any)
      pO <- sapply(lapply(strsplit(pT.all$OUTPUT, "; "), '%in%', pT), any)
      path <- c("plR_read", "plR_permute", "plR_rescale", "plR_prune")
      if (any(pO)) {
         if (any(pI)) {
            # create valid polylinkR function input message
            pI0 <- pT.all$FUNCTION[pI]
            if ("plR_prune" %in% pI0 & n.sets == 1) { # check if pruning possible
               pI0 <- setdiff(pI0, "plR_prune")
            }
            mss <- paste("output from", pT.all[pO]$FUNCTION)
            if (length(pI0) > 0) {
               mss <- paste(mss, "and valid input for",
                            paste(pI0, collapse = " or "))
            }

            # create analysis requirements message
            wmss <- "NOTE:"
            if (which(pO) == 1 & !cov.info) {
               w0 <- paste("Standard covariate columns not detected in",
                           "obj.info\n** gene score deconfounding requires",
                           "valid 'wt.mat' object (see ?polylinkR::plr_permute)")
               wmss <- c(wmss, w0)
            }
            if (which(pO) <= 2 %in% pI0 & !pos.info) {
               w0 <- paste("Standard coordinate columns not detected in",
                           "obj.info\n** gene set score decorrelation requires",
                           "valid 'ac' object (see ?polylinkR::plr_rescale)")
               wmss <- c(wmss, w0)
            }
         }
         # report plR object options
         cat(cli::col_cyan(paste("plR class object --", mss, "\n")))
         if (length(wmss) > 1) {
            cat(cli::col_red(paste(wmss, collapse = "\n"), "\n"))
         }
         # report plR object information
         for (name in c("set.info", "obj.info", "set.obj")) {
            pdg <- (min(80, getOption("width")) - nchar(name)) / 2
            cat(cli::boxx(name, padding = c(0, floor(pdg), 0, ceiling(pdg)),
                          col = "cyan", border_col = "cyan"))
            cat("\n")
            print(x[[name]], topn = 3)
         }
      } else {
         mss <- "plR object not produced by canonical polylinkR function\n"
         cat(cli::col_red(mss))
      }
   }
   invisible(x)
}


#' @title Summary method for \code{plR} objects
#'
#' @description
#' Provides a list of significant gene sets from the most recently applied
#' step of the \code{polylinkR} workflow.
#'
#' @param object
#'   \code{plR} class object; typically the output of a \code{polylinkR}
#'   function.
#'
#' @param sig
#'   \code{numeric}; significance threshold for filtering p- or q-values.
#'   Must be in the range \code{(0, 1)}. Defaults to \code{0.05}.
#'
#' @param ...
#'   Additional arguments passed to \code{summary}.
#'
#' @exportS3Method summary plR
#' @keywords internal
#'
#' @import data.table
#' @importFrom cli boxx col_cyan col_red
#'
#' @return
#' A summary of significant gene sets in the \code{plR} object.
#'
#' @examples
#' \dontrun{
#' # Basic usage (assuming `my_plR` is a valid plR object)
#' summary(my_plR)
#'
#' # Using a stricter significance threshold
#' summary(my_plR, sig = 0.01)
#' }
summary.plR <- function(object, sig = 0.05, ...) {
   if (sig <= 0 || sig >= 1) {
      stop("significance value (sig argument) must be between 0 and 1",
           call. = FALSE)
   } else {
       pT <- attributes(object)$plR.track # plR track info
       if (pT == "INVALID") { # check if empty plR object
          cat(cli::col_red("Empty plR object\n"))
       } else {
          acN <- colnames(object$set.info)
         gM <- grep("setScore", acN)
         if (length(gM) == 0) {
            cat(cli::col_red("No tests performed. Nothing to summarise\n"))
         } else {
            cN <- acN[max(gM)] # last p / q value column
            if (cN == "setScore.pr.q") { # additionally sort on p values
               SI <- object$set.info[order(get(cN), setScore.pr.p)]
            } else {
               SI <- object$set.info[order(setScore.pr.p)]
            }
            wSig <- SI[, get(cN) <= sig]
            nS <- sum(wSig, na.rm = TRUE)
            if (nS == 0) {
               mss <- paste0("No gene sets with ", cN, " <= ", sig)
               nn <- nrow(SI)
               if (nn >= 10) {
                  mss <- paste0(mss, "; showing top 10 gene sets:\n\n")
                  nn <- 10
               } else {
                  mss <- paste0(mss, "; showing all gene sets:\n\n")
               }
               cat(cli::col_cyan(mss))
               print(SI[1:nn])
            } else {
               mss <- paste0("Showing ", nS, " gene set",
                             ifelse(nS == 1, "", "s"), " with ", cN,
                             " <= ", sig, ":\n\n")
               cat(cli::col_cyan(mss))
                print(SI[wSig], nrows = nS)
            }
            invisible(SI)
         }
      }
   }
}


# Deprecated aliases for backward compatibility (v0.6.0)
# These will be removed in v1.0.0

#' @keywords internal
#' @noRd
.plR_check <- function(f, ENV) {
   .Deprecated(".check_plr_object", package = "polylinkR")
   .check_plr_object(f, ENV)
}

