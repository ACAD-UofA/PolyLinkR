#' @title Message printing based on verbose option check
#' @description Internal helper function the checks the `verbose` option in the
#'  parent function environment. Prints messages if this option set to `TRUE`.
#' @param m character; The message to print.
#' @importFrom rlang peek_option
#' @keywords internal
#' @noRd
.verbose_msg <- function(m) {
   if (isTRUE(rlang::peek_option("verbose"))) { # access verbose setting from parent function
      cat(m)
   }
}


#' @title Internal file path checker
#' @description An internal function to check for correct file paths for
#'  the required input files. User can either provide path to folder with all
#'  files or separate paths to each input file (applies to both required and
#'  optional files). Error returned if both options are used.
#' @param input.path character; Path to a directory containing all input files.
#' @param oi.path character; Path to the 'obj.info' file.
#' @param si.path character; Path to the 'set.info' file.
#' @param so.path character; Path to the 'set.obj' file.
#' @param rr.path character; Path to the 'rec.rate' file.
#' @param vi.path character; Path to the 'var.info' file.
#' @param group character; Identifying a specific label in file names.
#' @param ENV environment; The environment where results are returned.
#' @import data.table
#' @import foreach
#' @return No return value. Rather, if all required (and optional) paths are
#'  valid, a list of file paths is returned to the specified environment,
#'  otherwise the function stops with an error.
#' @keywords internal
#' @noRd
.check_file_paths <- function(oi.path, si.path, so.path, rr.path, vi.path,
                        input.path, group, ENV) {

   # check single or specific file paths
   emss <- c("input file paths\nEither specify path to each file ",
             "using obj.info.path, set.info.path, and set.obj.path arguments\n",
             "Or specify folder containing all files using file.path argument")
   missing_paths <- sapply(list(oi.path, si.path, so.path), is.null)
   if (is.null(input.path)) {
      if (any(missing_paths)) { # not all required separate paths specified
         stop("Incomplete ", emss, call. = FALSE)
      } else { # required paths specified
         multi.path <- TRUE
      }
   } else { # check
       if (any(!missing_paths)) { # not all required separate paths specified
         stop("Conflicting ", emss, call. = FALSE)
      } else { # required paths specified
         multi.path <- FALSE
      }
   }

   EMSS <- character(0) # collect error messages
   if (multi.path) { # check for individually labelled files
      fn0 <- c(oi.path, si.path, so.path)
      if (!is.null(vi.path)) {
         fn0 <- c(fn0, vi.path)
      }
      if (!is.null(rr.path)) {
         fn0 <- c(fn0, rr.path)
      }

      # check for non-empty files
      fE <- vapply(fn0, file.exists, FUN.VALUE = logical(1L))
      if (!all(fE)) {
         missing <- fn0[!fE]
         EMSS <- c(EMSS, paste(paste(missing, collapse = ", "), "does not exist\n"))
      }
   } else { # search in designated folder for polylingR files
      if (!dir.exists(input.path)) {
         EMSS <- c(EMSS, "input.path directory does not exist")
      } else {
         ll <- list.files(input.path)
         if (!is.null(group)) {
            ll <- ll[grep(group, ll)]
         }
         if (length(ll) == 0) {
            EMSS <- c(EMSS, "No input files detected")
         }
      }

      # check for too many or no copies of each requested input file
      fn0 <- NULL # collect valid file names
      ll <- file.info(list.files(input.path, full.names = TRUE))
      fl <- list.files(input.path)[!ll$isdir] # remove directories
      req.files <- c("ObjInfo", "SetInfo", "SetObj", "VarInfo", "RecRate") # valid file paths and names

      for (rfi in req.files) {
         # check multiple alternatives
         input.names <- foreach::foreach(j = 1:2, .combine = c) %do% {
            a <- ifelse(j == 1, substr(rfi, 1, 3), tolower(substr(rfi, 1, 3)))
            b <- ifelse(j == 1, substr(rfi, 4, 7), tolower(substr(rfi, 4, 7)))
            unlist(foreach::foreach(k = c("", "\\.", "_")) %do% {
               paste(a, b, sep = k)
            })
         }

         # check different file name options
         af0 <- sapply(input.names, FUN = function(x) {
            FL <- fl[grep(x, fl)]
            if (length(FL) > 0 & !is.null(group)) {
               FL <- FL[grep(group, FL)]
            }
            return(FL)
         })
         nT <- sapply(af0, length)

         if (sum(nT) == 1) { # one version of file
            fn0 <- c(fn0, file.path(input.path, af0[nT == 1]))
         } else {
            if (sum(nT) > 1) { # multiple file versions
               mss0 <- paste("Multiple copies of", rfi, "detected\n")
               EMSS <- c(EMSS, mss0)
            } else {
               if (sum(nT) == 0 & rfi %in% req.files[1:3]) { # no versions of file (ignore for rec.rate and var.info files)
                  mss0 <- paste("No copies of", rfi, "detected\n")
                  EMSS <- c(EMSS, mss0)
               }
            }
         }
      }
   }

   if (length(EMSS) > 0L) {
      stop(paste0(EMSS),
           "Check path and ensure that file names are correct", call. = FALSE)
   }

   file.names <- c("obj.info", "set.info", "set.obj", "var.info", "rec.rate")
   names(fn0) <- file.names[1:length(fn0)]
   assign(x = "file.paths", value = fn0, envir = ENV)
}


#' @title Return required datasets to the parent environment
#' @description Internal function to read data from various file formats.
#' @param plr character; plR class object.
#' @param ENV environment; The environment in which the `plR` object is stored.
#' @return No return value. Instead, all required primary and auxiliary
#'  datasets from plR object are unpacked within the specified environment.
#' @keywords internal
#' @noRd
.unpack_plr_object <- function(plr, ENV) {
   list2env(plr, envir = ENV) # unpack core objects
   p.unp <- c("plR.data", "plR.args", "plR.summary", "plR.seed",
              "plR.track", "plR.session")
   for (i in p.unp) {
      plr.i <- attr(plr, i)
      if (i == "plR.data") { # unpack ancillary datasets
         for (plr.j in plr.i) {
            for (j in names(plr.j)) {
               plr.k <- plr.j[[which(names(plr.j) == j)]]
               assign(x = j, value = plr.k, envir = ENV)
            }
         }
      } else {
         assign(x = gsub("plR", "plr", i), value = plr.i, envir = ENV)
      }
   }
}


#' @title Validate arguments for core polylinkR functions
#' @description Internal function that performs checks on the formal arguments
#'  of the parent `polylinkR` function before analysis. Checks are based on the
#'  target function (`f`).
#' @param f character; Name of the `polylinkR` function.
#' @param ENV environment; The environment in which the `plR` object is stored.
#' @import data.table
#' @importFrom gstat vgm
#' @return No return value. Rather, the function either stops with a
#'  descriptive error message if a check fails or, if all checks pass, it
#'  assigns validated objects, messages, and warnings back to the specified
#'  environment (`ENV`).
#' @keywords internal
#' @noRd
.check_arguments <- function(f, ENV) {
   f0 <- ifelse(f == "plot", paste0(f, ".plR"), paste0("plR_", f)) # function names
   f0.args <- names(formals(get(f0, envir = asNamespace("polylinkR")))) # retrieve argument names
   list2env(mget(setdiff(f0.args, "..."), envir = ENV), envir = environment()) # collect function arguments
   obj2env <- list() # collect objects to return to function environment
   MSS <- NULL # collect messages
   WMSS <- NULL # collect warnings

   if (f %in% c("permute", "rescale", "prune")) {
      plr_input_obj <- get("plR.input", envir = ENV)
      attrs <- attributes(plr_input_obj)
      # Handle both old and new attribute names for backward compatibility
      track <- if (!is.null(attrs$plr_track)) attrs$plr_track else attrs$plR.track
      list2env(plr_input_obj, envir = environment()) # retrieve input objects
      plr_data <- if (!is.null(attrs$plr_data)) attrs$plr_data else attrs$plR.data
      read_data <- if (!is.null(plr_data$read_data)) plr_data$read_data else plr_data$read.data
      params <- read_data
      list2env(params, envir = environment())
   }

   if (f == "read") { # check plR_read argumants
      # check lengths of single value options
       single_value_params <- c("min.set.n", "max.set.n", "group", "set.merge", "map.fun",
               "obj.stat.fun", "obj.buffer", "bin.size")
       valid_params <- sapply(single_value_params, function(x) ifelse(is.null(get(x)), TRUE,
                                            length(get(x)) == 1))
       if (!all(valid_params)) {
          stop("Length of ", paste(single_value_params[!valid_params], collapse = " & "), " != 1",
               call. = FALSE)
       }

      e.mss <- paste0("set.merge incorrectly specified [",
                      "must be numeric value between 0 and 1]\n")
      if(!is.numeric(set.merge)) {
         stop(e.mss, call. = FALSE)
      } else {
         if (set.merge > 1 | set.merge <= 0) {
            stop(e.mss, call. = FALSE)
         }
      }

      if (!is.numeric(min.set.n)) {
         stop("min.set.n incorrectly specified ",
              "[must be numeric value between 2 and Inf]", call. = FALSE)
      }

      if (!is.numeric(max.set.n)) {
         stop("max.set.n incorrectly specified [",
              "must be numeric value between 10 and Inf]", call. = FALSE)
      }

      # check set size filtering criteria
      ss.ch <- c(min.set.n < 2, max.set.n < 10, min.set.n > max.set.n)
      mss.ch <- c("min.set.n < 2", "max.set.n < 10", "min.set.n > max.set.n")
      if (any(ss.ch)) {
         W <- mss.ch[which(ss.ch)]
         stop("Set size filters incorrectly specified: ",
              paste(W, collapse = " and "), call. = FALSE)
      }

      if (!is.logical(rem.genes) | length(rem.genes) > 1){
         stop("rem.genes must be logical", call. = FALSE)
      }

      # check gene set and gene filters
      if (!is.null(set.in) & !is.null(set.out)) {
         stop("Use either set.in or set.out option, not both", call. = FALSE)
      }

      if (!is.null(obj.in) & !is.null(obj.out)) {
         stop("Use either obj.in or obj.out option, not both", call. = FALSE)
      }

      # check gene score generation functions
      if (obj.buffer %in% c("auto", "Auto", "AUTO")) {
         obj2env$obj.buffer <- "auto"
      } else {
         emss <- "obj.buffer incorrectly specified: "
         if (is.numeric(obj.buffer)) {
            if (obj.buffer < 0 | obj.buffer > 1e5) {
               stop(emss, "value outside acceptable range [0, 1e5]",
                    call. = FALSE)
            }
         } else {
            stop(emss, "non-numeric value used", call. = FALSE)
         }
      }

      of.names <- c("non.param", "lm.logN")
      if (!(obj.stat.fun %in% of.names)) {
         stop(c("obj.stat.fun is not valid\nOptions: ",
                paste(of.names, collapse = ", ")), call. = FALSE)
      }

      # check gene score generation functions
      emss <- "bin.size incorrectly specified: "
      if (is.numeric(bin.size)) {
         if (bin.size < 50 | bin.size > 1e3) {
            stop(emss, "value outside acceptable range [50, 1e3]", call. = FALSE)
         }
      } else {
         stop(emss, "non-numeric value used", call. = FALSE)
      }

      # check genetic mapping function
      mf.names <- c("haldane", "kosambi")
      if (!(map.fun %in% mf.names)) {
         stop("map.int not appropriately specified\n",
              "Must be either 'kosambi' or 'haldane'", call. = FALSE)
      }
   }

   if (f == "permute") { # check plR_permute arguments
      # check lengths of single value options
       single_value_params <- c("md.meth", "n.perm", "n.boot", "alt", "kern.func", "kern.scale",
               "kern.wt.max", "gpd.cutoff", "n.cores", "fut.plan", "permute")
       valid_params <- sapply(single_value_params, function(x) ifelse(is.null(get(x)), TRUE,
                                            length(get(x)) == 1))
       if (!all(valid_params)) {
          stop("Length of ", paste(single_value_params[!valid_params], collapse = " & "), " != 1",
               call. = FALSE)
       }

      # check for violations in logical settings
      if (!is.logical(permute)) {
         stop("permute parameter must be logical", call. = FALSE)
      }

      # check appropriate testing choices
      if (!(alt %in% c("upper", "lower"))) {
         stop("alt parameter should be either 'upper' or 'lower'\n",
              "i.e. testing for large set scores in upper tail or ",
              "small set scores in lower tail, respectively", call. = FALSE)
      }

      # check for appropriate local regression parameter choices
      if (track == "000") { # first run of plR_permute
         if (cov.info) { # covariates must be present
            if (!(md.meth %in% c("raw", "robust"))) {
               stop("kern.func must be either 'robust' (covariates converted ",
                    "to normal scores) or 'raw' (base covariate values)",
                    call. = FALSE)
            }

            all.func <- c("normal", "exponential", "inverse")
            if (!(kern.func %in% all.func)) { # determine md.meth
               stop("kern.func must be either '",
                    paste(all.func, collapse = "' or '"), "'", call. = FALSE)
            }

            if ((kern.bound %in% c("auto", "Auto", "AUTO"))) {
               obj2env$kern.bound <- "auto"
            } else {
               e.mss <- "kern.bound must numeric value > 0"
               if (is.numeric(kern.bound)) {
                  if (kern.bound <= 0) {
                     stop(e.mss, call. = FALSE)
                  }
               } else {
                  stop(e.mss, call. = FALSE)
               }
            }

            if (kern.scale %in% c("auto", "Auto", "AUTO")) {
               obj2env$kern.scale <- "auto"
            } else {
               e.mss <- "kern.scale must numeric value > 0"
               if (is.numeric(kern.scale)) {
                  if (kern.scale <= 0) {
                     stop(e.mss, call. = FALSE)
                  } else {
                     if (kern.scale >= 100) {
                        wms <- paste("kern.scale >= 100;",
                                     "Large scale can increase bias",
                                     "in corrected gene scores\n")
                     } else {
                        wms <- paste("kern.scale <= 0.01;",
                                     "Small scale can increase variance",
                                     "in corrected gene scores\n")
                     }
                     WMSS <- c(WMSS, wms)
                  }
               } else {
                  stop(e.mss, call. = FALSE)
               }
            }

            uni.p <- 1 / (n.genes - 1)
            e.mss <- paste("kern.wt.max must numeric value >",  uni.p,
                           "[i.e. 1 / (n.genes - 1)] and <= 1")
            if (is.numeric(kern.wt.max)) {
               if (kern.wt.max <= uni.p | kern.wt.max > 1) {
                  stop(e.mss, call. = FALSE)
               }
            } else {
               stop(e.mss, call. = FALSE)
            }

            obj2env$perm.path <- "full" # perform deconfounding
         } else {
            if (permute) { # permutation requested
               MSS <- c(MSS,
                        paste("Note: no covariate information detected;",
                              "confounder correction not possible\n"))
               obj2env$md.meth <- "none"
               obj2env$kern.func <- "none"
               obj2env$perm.path <- "partial" # no deconfounding
            } else { # no permutation requested
               stop("Deconfounding already performed but permute = FALSE. ",
                    "Nothing to do.",  call. = FALSE)
            }
         }
      } else { # deconfounding not performed
         if (permute) { # deconfounding performed, permutation step required
            obj2env$perm.path <- "skip" # go straight to permutation step
         } else {
            stop("Deconfounding already performed but permute = FALSE. ",
                 "Nothing to do.",  call. = FALSE)
         }
      }

      # permuation parameters
      if (permute) {
         e.mss <- paste0("n.perm incorrectly specified\n",
                         "n.perm must be integer value >= 1e5 ",
                         "and be exactly divisible by 1e5")
         if (!is.numeric(n.perm)) {
            stop(e.mss, call. = FALSE)
         } else {
            if (n.perm < 1e5 | n.perm %% 1e5 != 0) {
               stop(e.mss, call. = FALSE)
            }
            if (n.perm >= 1e7) {
               wmss <- paste(n.perm, "permuations requested; operations",
                             "may be slow without sufficient parallelisation")
               WMSS <- c(WMSS, wmss)
            }
         }

         e.mss <- paste0("n.boot incorrectly specified\n",
                         "n.boot must be integer value >= 10")
         if (!is.numeric(n.boot)) {
            stop(e.mss, call. = FALSE)
         } else {
            if (n.boot < 10) {
               stop(e.mss, call. = FALSE)
            }
            if (n.boot >= 100) {
               wmss <- paste(n.boot,
                             "bootstrap replicates requested; operations",
                             "may be slow without sufficient parallelisation")
               WMSS <- c(WMSS, wmss)
            }
         }

         n.mss <- paste("gpd.cutoff incorrectly specified\nMust be value in",
                        "range [1e3 / n.perm, min(5e4 / n.perm, 5e-2)]")
         if (is.numeric(gpd.cutoff)) {
            max.gpd <- min(5e4 / n.perm, 5e-2)
            if (gpd.cutoff < 1e3 / n.perm | gpd.cutoff > max.gpd) {
               stop(n.mss, call. = FALSE)
            }
         } else {
            stop(n.mss, call. = FALSE)
         }
      }
   }

   if (f == "rescale") { # check plR_rescale arguments
      # check lengths of single value options
       single_value_params <- c("cgm.bin", "cgm.range", "cgm.wt.max", "emp.bayes", "min.rho",
               "rescale", "fast", "n.cores", "fut.plan")
       valid_params <- sapply(single_value_params, function(x) ifelse(is.null(get(x)), TRUE,
                                            length(get(x)) == 1))
       if (!all(valid_params)) {
          stop("Length of ", paste(single_value_params[!valid_params], collapse = " & "), " != 1",
               call. = FALSE)
       }

      # check for violations in logical settings
      lp <- c("rescale", "fast")
      clp <- sapply(lp, function(x) is.logical(get(x)))
      if (!all(clp)) {
         stop(paste(lp[!clp], collapse = " & "), " parameter must be logical",
              call. = FALSE)
      }

      track.vec <- strsplit(track, "")[[1]]
      if (track.vec[2] == 1) { # autocovariance already estimated
         if (!rescale) {
            stop("ac object present but rescale = FALSE; ",
                 "nothing to do\n", call. = FALSE)
         }
      } else { # autocovariance not estimated
         if (is.null(ac)) { # check for violations in parameter settings for variogram estimation
            if (!pos.info) {
               stop("No covariate information and no ac object provided: ",
                    "rescaling not possible", call. = FALSE)
            }

            # check minimum correlation coefficient (min.rho)
            mr.mss <- c("min.rho incorrectly entered\n",
                        "Must be numerical value between 0 and 0.001")
            if (!is.numeric(min.rho) | length(min.rho) > 1) {
               stop(mr.mss, call. = FALSE)
            } else {
               if (min.rho < 0 | min.rho > 0.001) {
                  stop(mr.mss, call. = FALSE)
               }
            }

            # check minimum lag bin sizes in correlogram estimation
            if (!is.numeric(cgm.bin) | cgm.bin < 10 | cgm.bin > 1e3) {
               stop("cgm.bin must be numeric value >= 10 and <= 100",
                    call. = FALSE)
            } else {
               obj2env$cgm.bin <- as.integer(round(cgm.bin)) # force integer value
            }

            # check maximum lag used in correlogram estimation
            if (cgm.range %in% c("auto", "Auto", "AUTO")) {
               obj2env$cgm.range <- "auto"
            } else {
               upM <- ifelse(coord == "cM", 5, 5e7)
               lowM <- ifelse(coord == "cM", 0.01, 1e5)
               if (!is.numeric(cgm.range) | cgm.bin < 10 | cgm.bin > 1e3) {
                  stop("cgm.range must be numeric value >= ", lowM, coord,
                       " and <= ", upM, coord, call. = FALSE)
               }
            }

            # check correlogram fitting weights
            if (!is.numeric(cgm.wt.max) | cgm.wt.max <= 0 | cgm.wt.max > 1) {
               stop("cgm.wt.max must numeric value > 0 and <= 1", call. = FALSE)
            }

            # check emprical Bayes choice
            if (emp.bayes %in% c("auto", "Auto", "AUTO")) {
               obj2env$emp.bayes <- "auto"
            } else {
               if (emp.bayes %in% c("full", "reduced")) {
                  if (emp.bayes == "full" & n.chr < 15) {
                     w.mss <- paste0("emp.bayes = 'full' but only ",  n.chr,
                                     " chromosomes present\nReduced model will",
                                     " also be fit and best model determined\n")
                     WMSS <- c(WMSS, w.mss)
                  }
               } else {
                  stop("emp.bayes must be one of 'auto', 'full' or 'reduced'",
                       call. = FALSE)
               }
            }
         } else { # user has provided autocovariance object
            if (rescale) { # check user autocovariance object
               exp.col <- c("objID.A", "objID.B", "distance", "rho")
               emm <- paste0("ac must be a data.table or data.frame ",
                             "with the following columns: '",
                             paste(exp.col, collapse = "', '"), "'",
                             "\nAND include separate rows for all ",
                             "gene pairs with non-zero autocovariance")
               if (!is.data.frame(ac)) {
                  stop(emm, call. = FALSE)
               } else {
                  data.table::setDT(ac) # ensure data.table format
                  if (!all(exp.col %in% colnames(ac))) {
                     stop(emm, call. = FALSE)
                  } else {
                     if (any(ac$rho <= 0)) {
                        stop("corrletions for all genes pairs must be ",
                             "greater than 0", call. = FALSE)
                     }
                  }
                  obj2env$ac <- ac
               }
            } else { # no autocovariance object and rescale option = FALSE
               stop("ac object present but rescale = FALSE; ",
                    "nothing to do\n", call. = FALSE)
            }
         }
      }
   }

   if (f == "prune") { # check plR_prune arguments
      if (no.share) { # check for shared genes
         if (n.sets == 1) {
            emss <- "only single gene set present"
         } else {
            emss <- "no shared genes across gene sets"
         }
         stop(paste("No need for pruning --", emss), call. = FALSE)
      }

      # check lengths of single value options
       single_value_params <- c("n.fdr", "tolerance", "n.cores", "fut.plan")
       valid_params <- sapply(single_value_params, function(x) ifelse(is.null(get(x)), TRUE,
                                            length(get(x)) == 1))
       if (!all(valid_params)) {
          stop("Length of ", paste(single_value_params[!valid_params], collapse = " & "), " != 1",
               call. = FALSE)
       }

      # check logical parameters
      if (!is.logical(est.pi0)) {
         stop("est.pi0 must be logical value", call. = FALSE)
      }

      e.mss <- paste0("n.fdr incorrectly specified\n",
                      "n.fdr must be integer value >= 100",
                      " and be exactly divisible by 100")
      if (!is.numeric(n.fdr) | length(n.fdr) > 1) {
         stop(e.mss, call. = FALSE)
      } else {
         if (n.fdr < 100 | n.fdr %% 100 != 0) {
            stop(e.mss, call. = FALSE)
         } else {
            if (n.fdr > 1000) {
               w.mss <- paste("n.fdr > 1000; pruning may be slow ",
                              "without sufficient parallelisation\n")
               WMSS <- c(WMSS, w.mss)
            }
         }
      }

      # check tolerance on pi0 convergence
      if (est.pi0) {
         if (!is.numeric(tolerance)) {
            stop("tolerance must be numeric\n", call. = FALSE)
         } else {
            if (tolerance <= 0 | tolerance > 0.01) {
               stop("tolerance must be in interval (0, 0.01]\n", call. = FALSE)
            }
         }
      }
   }

if (f == "plot") { # check plotting arguments
       x_obj <- get("x", envir = ENV)
       x_attrs <- attributes(x_obj)
       # Handle both old and new attribute names for backward compatibility
       track <- if (!is.null(x_attrs$plr_track)) x_attrs$plr_track else x_attrs$plR.track
       x_plr_data <- if (!is.null(x_attrs$plr_data)) x_attrs$plr_data else x_attrs$plR.data
       x_read_data <- if (!is.null(x_plr_data$read_data)) x_plr_data$read_data else x_plr_data$read.data
       params <- x_read_data
       list2env(params, envir = environment())

       single_value_params <- c("log.dist", "log.Ne", "log.gamma", "output.path", "plot.name",
               "max.facets", "plot.all")
       valid_params <- sapply(single_value_params, function(x) ifelse(is.null(get(x)), TRUE,
                                            length(get(x)) == 1))
       if (!all(valid_params)) {
          stop("Length of ", paste(single_value_params[!valid_params], collapse = " & "), " != 1",
               call. = FALSE)
       }

      if (!is.null(output.path)) {
         if (!dir.exists(output.path)) {
            stop("output.path does not exist\n", call. = FALSE)
         }
      }

      if (!is.null(plot.name)) {
         if (!is.character(plot.name)) {
            stop("plot.name must be character string\n", call. = FALSE)
         }
      }

      # check common logical parameters
      lp <- c("verbose", "args.to.title", "show.Ne", "log.dist",
              "log.Ne", "log.gamma", "plot.all")
      clp <- sapply(lp, function(x) is.logical(get(x)))
      if (!all(clp)) {
         stop(paste(lp[!clp], collapse = " & "), " must be logical",
              call. = FALSE)
      }

      # check facet specification
      if (is.numeric(max.facets)) {
         if (max.facets < 1 | max.facets > 50) {
            stop("max.facets must be >= 1 and <= 50", call. = FALSE)
         } else {
            obj2env$max.facets <- ceiling(max.facets) # ensure integer
         }
      } else {
         stop("max.facets must be integer value >= 1 and <= 50",
              call. = FALSE)
      }

      # check highlighted sets option
      if (!is.null(show.sets) & track != "000") {
         emss <- paste0("show.sets incorrectly specified\n",
                        "Must be vector where each entry is a set ID")
         hs.check <- show.sets %in% get("x", envir = ENV)$set.info$setID
         if (!all(hs.check)) {
            hsN <- sum(!hs.check)
            if (hsN <= 20) {
               stop("Some gene IDs [",
                    paste(show.sets[!hs.check], collapse = ", "),
                    "] in show.sets are not present in set.info",
                    call. = FALSE)
            } else {
               stop(hsN,
                    "set IDs in show.sets are not present in set.info",
                    call. = FALSE)
            }
         }
      }

      # determine output type
      track.vec <- as.numeric(strsplit(track, "")[[1]])
      if (!plot.all) { # only plot output for most recent function
         ss <- sum(track.vec > 0)
         if (ss > 1) {
            for (i in 1:(ss - 1)) { # ignore output of prior plR_functions
               track.vec[i] <- 0
            }
         }
      }

      x_attrs2 <- attributes(get("x", envir = ENV))
      x_plr_summary <- if (!is.null(x_attrs2$plr_summary)) x_attrs2$plr_summary else x_attrs2$plR.summary
      x_read_summary <- if (!is.null(x_plr_summary$read_summary)) x_plr_summary$read_summary else x_plr_summary$read.summary
      rsum <- x_read_summary
      osp <- !is.null(rsum$obj.std.param) # gene scores were computed in plR_read

      if (track.vec[1] %in% 1:2) { # deconfounding was performed
         # check that covariates are correctly defined
         if (!is.null(cov.labels)) {
            if (length(cov.labels) != n.cov) {
               stop(paste0("length of cov.labels [", length(cov.labels),
                           "] != number of covariate columns in obj.info [" ,
                           n.cov, "]\n"), call. = FALSE)
            }
         }
         if (!is.null(show.sets) & track == "100") {
            show.sets <- NULL
            MSS <- c(MSS,
                     "show.sets requires p values from enrichment tests")
         }
      } else { # no deconfounding plots
         if (sum(track.vec) == 0) { # read data only
            if (osp) {
               if (!is.null(show.sets)) {
                  show.sets <- NULL
                  MSS <- c(MSS,
                           "show.sets requires p values from enrichment tests")
               }
            } else {
               stop("Nothing to plot", call. = FALSE)
            }
         } else { # additional plR operations performed
            if (track == "300" & is.null(show.sets)) { # permuted but no deconfounding
               stop("Nothing to plot", call. = FALSE)
            }
         }
      }

      # define plot outputs
      obj2env$obj.std <- ifelse(track == "000", osp, osp & plot.all)
      obj2env$obj.dist <- !is.null(show.sets) & track.vec[1] > 1

      obj2env$deconf <- track.vec[1] %in% 1:2
      obj2env$permuted <- track.vec[1] %in% 2:3
      obj2env$vario <- track.vec[2] %in% 1:2
      obj2env$rescaled <- track.vec[2] %in% 2:3
      obj2env$pruned <- track.vec[3] > 0

   }

   # Return info to function environment
   if (!is.null(MSS)) MSS <- paste(MSS, collapse = "\n")
   assign("MSS", MSS, envir = ENV)
   if (!is.null(WMSS)) WMSS <- paste(WMSS, collapse = "\n")
   assign("WMSS", WMSS, envir = ENV)

   # return information and warning messages to parent environment
   list2env(obj2env, envir = ENV)
}


#' @title Track polylinkR output types
#' @description Internal function to create a table of valid inputs and
#'  outputs for all core `polylinkR` functions.
#' @import data.table
#' @return A data.table with a summary of `polylinkR` functions and their inputs
#'  / outputs.
#' @keywords internal
#' @noRd
.plR_track <- function() {
   # define possible output tracks for each function
   plr_track_read <- "000"
   plr_track_permute <- paste0(1:3, 0, 0)
   plr_track_rescale <- paste0(outer(2:3, 1:3, paste0), 0)
   plr_track_prune <- paste0(outer(2:3, c(0, 2:3), paste0), 1)
   plr_track_names <- c("read_polylinkr_data", "permute_polylinkr_data", "rescale_polylinkr_data", "prune_polylinkr_data")
   # assign valid inputs
   in0 <- c("Input data",
            paste(c(plr_track_read, plr_track_permute[1]), collapse = "; "),
            paste(c(plr_track_permute[-1], plr_track_rescale[1:2]), collapse = "; "),
            paste(c(plr_track_permute[-1], plr_track_rescale[-(1:2)]), collapse = "; "))
   # assign valid outputs
   out0 <- c(paste(plr_track_read, collapse = "; "),
              paste(plr_track_permute, collapse = "; "),
              paste(plr_track_rescale, collapse = "; "),
              paste(plr_track_prune, collapse = "; "))
   plr.tab <- data.table::data.table(FUNCTION = plr_track_names, INPUT = in0,
                                      OUTPUT = out0)
   return(plr.tab)
}


#' @title Check and optionally create random seed
#' @description Internal function that either validates user provided random
#'  seed or generates a new random integer seed.
#' @param user.seed numeric; An optional value to be used as a random seed.
#' @return A single integer value to be used as a random seed.
#' @noRd
.create_seed <- function(user.seed){
   if (!is.null(user.seed)) {
      emss <- c("Incorrect random seed specification:\n",
                "Must be single numeric value x such that abs(x) <= ",
                .Machine$integer.max)
      if (length(user.seed) != 1 | !is.numeric(user.seed)) {
         stop(emss, call. = FALSE)
      } else {
         if (abs(user.seed) > .Machine$integer.max) {
            stop(emss, call. = FALSE)
         } else {
            seed <- as.integer(user.seed) # enforce integer value
         }
      }
   } else {
      seed <- sample.int(.Machine$integer.max, 1) * sample(c(-1, 1), 1)
      seed <- as.integer(seed)
   }
   return(seed)
}


#' @title Measure and format elapsed time
#' @description Internal helper function that calculates and formats the
#'  elapsed time from a given start point, returning a human-readable string.
#' @param start.time character; `POSIXct` object representing the starting time.
#' @return A character string representing the elapsed time in "Hh Mm Ss"
#'  format.
#' @noRd
.get_time <- function(start.time) {
   TT <- difftime(Sys.time(), start.time, units = "secs")
   tt <- ceiling(as.numeric(TT))
   H <- tt %/% 3600
   rr <- tt %% 3600
   if(rr > 0) {
      M <- rr %/% 60
      rr2 <- rr %% 60
      if(rr2 > 0) {
         S <- round(rr2)
      } else {
         S <- 0
      }
   } else {
      M <- 0
      S <- 0
   }
   return(list(paste0(H, "h ", M, "m ", S, "s"), TT))
}


#' @title Set up parallel processing parameters
#' @description Internal function that validates and configures the parameters
#'  for parallel processing. Checks the number of cores requested and the
#'  future plan, making adjustments based on system capabilities and user input.
#' @param n.cores character or integer; Controls the number of cores used in
#'  parallel processes. Default option (`auto`) uses `max(1, no. free cores - 1)`.
#'  Users can also choose `max` to automatically detect and use all available
#'  cores. Overrides incompatible `fut.plan` option.
#' @param fut.plan character; String specifying the future plan.
#'   Default option (`auto`) for multiple cores is `multisession`, with
#'   `sequential` used for single core machines. Overridden by `n.cores` when
#'   incompatible option is chosen (e.g. if user chooses `n.cores` = 5 and
#'   `fut.plan` = `sequential`, `multisession` processing with 5 cores is used).
#' @param verbose logical; Determines whether progress messages are displayed.
#' @param ENV environment; Environment where the parameters will be stored.
#' @importFrom future availableCores supportsMulticore
#' @importFrom progressr handlers handler_progress handler_void
#' @return No return value. Instead, the function assigns the validated parallel
#'   processing parameters to the specified environment, or stops with an error.
#' @noRd
.par_params <- function(n.cores, fut.plan, verbose, ENV) {
   pWMSS <- NULL # capture warning messages

   # detect available cores
   N.cores <- future::availableCores()

   # check n.cores option
   if (n.cores %in% c("auto", "Auto", "AUTO")) {
      nc <- N.cores - 1
   } else {
      if (n.cores %in% c("max", "Max", "MAX")) {
         nc <- N.cores
      } else {
         if (is.numeric(n.cores)) {
            nc <- as.integer(n.cores) # force integer
            if (n.cores <= 0 | n.cores > N.cores) {
               nc <- N.cores - 1
               pWMSS <- c(pWMSS,
                          paste0("n.cores not in accepted range [1, ",
                                 N.cores, "]"))
            }
         } else {
            stop(paste0("n.cores incorrectly specified\n",
                        "must be either positive integer, 'auto' or 'max'"),
                 call. = FALSE)
         }
      }
   }

   # check fut.plan option
   valid.fp <- c("auto", "Auto", "AUTO", "sequential", "multisession",
                 "multicore")
   if (!fut.plan %in% valid.fp) {
      stop("fut.plan incorrectly specified\nMust be one of: ",
           paste(valid.fp, collapse = ", "),
           call. = FALSE)
   } else {
      if (fut.plan %in% c("auto", "Auto", "AUTO")) {
         fut.plan <- ifelse(nc == 1L, "sequential", "multisession")
      } else {
         if (fut.plan == "sequential") { # user has chosen sequential processing
            nc <- 1
         } else { # user has chosen parallel processing
            if (nc == 1) {
               fut.plan <- "sequential" # forcing sequential processing
               if (N.cores == 1) {
                  ms <- "Only single core detected; using sequential processing"
                  pWMSS <- c(pWMSS, ms)
               } else {
                  ms <- "User has set n.core = 1; using sequential processing"
                  pWMSS <- c(pWMSS, ms)
               }
            } else {
               if (fut.plan == "multicore") { # check for multicore support
                  if(!future::supportsMulticore()) {
                     fut.plan <- "multisession"
                     ms <- "multicore option not supported; using multisession"
                     pWMSS <- c(pWMSS, ms)
                  }
               }
            }
         }
      }
   }

   if (fut.plan == "sequential") {
      pMSS <- c("Sequential processing enabled [fut.plan = sequential]",
                "All processes will utilise a single core")
   } else {
      pMSS <- c(paste0("Parallel processing enabled [fut.plan = ",
                       fut.plan, "]"),
                paste0("Specific processes will utilise ", nc, " cores"))
   }

   pMSS <- paste(pMSS, collapse = "\n")
   if (!is.null(pWMSS)) pWMSS <- paste(pWMSS, collapse = "\n")

   if (verbose) { # use progress reporting
      preset.hand <- progressr::handlers(global = NA) # check for preset handlers
      if (preset.hand) { # use predefined progress handler
         prog.hand <- progressr::handlers()
      } else {
         fmt <- "[:bar] :percent in :elapsed | ETA: :eta"
         prog.hand <- progressr::handler_progress(format = fmt)
      }
   } else { # no progress reporting
      prog.hand <- progressr::handler_void()
   }

   list2env(list(n.cores = nc, fut.plan = fut.plan, prog.hand = prog.hand,
                 pWMSS = pWMSS, pMSS = pMSS), envir = ENV)
}


#' @title Process and standardise primary datasets
#' @description Internal function that standardises the `obj.info`, `set.info`,
#'   and `set.obj` data.table objects. Ensures that the object IDs and set IDs
#'   are contiguous integers, simplifing downstream processing. Also adds a
#'   `midpos` column to `obj.info`f positional information is present.
#' @param OI data.table; The gene (object) information data.table.
#' @param SI data.table; The gene set information data.table.
#' @param SO data.table; The gene to gene set mapping data.table.
#' @param pos.info logical; Is position info provided in `obj.info`?
#' @param ENV environment; Environment where the processed data frames will be
#'   stored.
#' @import data.table
#' @return No return value. Instead, the function returns standardized
#'   data frames and a vector of unique genes to the specified environment.
#' @noRd
.file_set <- function(OI, SI, SO, pos.info, ENV) {
   # copy files to avoid overwriting
   OI <- data.table::copy(OI)
   SI <- data.table::copy(SI)
   SO <- data.table::copy(SO)

   OI[, objID.orig := objID]
   if (pos.info) {
      # create unique integer for each chromosome
      OI[, chr.orig := chr]
      OI[, chr := as.integer(factor(chr.orig),
                             levels = sort(unique(chr.orig)))]
      data.table::setkey(OI, chr, midpos) # order by position
   } else {
      data.table::setkey(OI, objID.orig)
   }
   OI[, objID := .I] # create unique integer for each objID

   # ensure setIDs are contiguous and SI is ordered by setIDs
   SI[, setID.orig := setID]
   data.table::setkey(SI, setID.orig)
   SI[, setID := .I]

   # regenerate set.obj
   SO[, `:=` (setID.orig = setID, objID.orig = objID)]
   SO[, `:=` (setID = NULL, objID = NULL)]
   data.table::setkey(SO, setID.orig)
   SO[.(SI[, .(setID.orig, setID)]), setID := setID]
   data.table::setkey(SO, objID.orig)
   SO[.(OI[, .(objID.orig, objID)]), objID := objID]
   SO <- SO[order(setID, objID)]

   # determine unique genes in gene sets
   set.genes <- sort(unique(SO$objID))

   # return updated objects to function evnironment
   list2env(list(obj.info = OI, set.info = SI, set.obj = SO,
                 set.genes = set.genes), envir = ENV)
}


#' @title Restore original file formats of primary datasets
#' @description Internal function to reset the `obj.info`, `set.info`, and
#'   `set.obj` data.table objects to their original format
#' @param OI data.table; The gene (object) information data.table.
#' @param SI data.table; The gene set information data.table.
#' @param SO data.table; The gene to gene set mapping data.table.
#' @param pos.info logical; Is position info provided in `obj.info`?
#' @import data.table
#' @return No return value. Instead, the function transforms data.tables
#'   within specified environment.
#' @noRd
.file_reset <- function(OI, SI, SO, pos.info) {
   # remove keys
   data.table::setkey(OI, NULL)
   data.table::setkey(SI, NULL)
   data.table::setkey(SO, NULL)

   # reset columns
   OI[, objID := objID.orig]
   OI[, objID.orig := NULL]

   if (pos.info) {
      OI[, chr := chr.orig]
      OI[, chr.orig := NULL]
   }

   SI[, setID := setID.orig]
   SI[, setID.orig := NULL]

   SO[, setID := setID.orig]
   SO[, setID.orig := NULL]
   SO[, objID := objID.orig]
   SO[, objID.orig := NULL]
}


#' @title Calculate gene scores from SNP / binned summary statistics
#' @description Internal function that computes non-parametric residuals for
#'  gene scores, accounting for differences in overlapping statistic density.
#' @param OI data.table; The object (gene) information table.
#' @param VI data.table; The variant information table containing summary
#'  statistics computed across SNPs or intervals.
#' @param FUN character; function used to calculate gene scores. Default is
#'  `non.param`, a robust non-parametric method, that uses binned data to
#'  calculate median and median absolute deviation (MAD) to normalise scores.
#'  Alternatively, `lm.logN` applies a linear regression to the log-transformed
#'  SNP / bin counts (assumes a roughly linear relationship is appropriate). In
#'  both cases, expected scores are estimated and gene scores calculated as the
#'  residual value.
#' @param binN numeric; The minimum number of genes to include in each bin. Not
#'  used for `lm.logN` option.
#' @import data.table
#' @return list containing data.table with gene scores and separate data.table
#'  with raw maximum scores and model parameters as an attribute.
#' @noRd
.get_obj_stat <- function(OI, VI, FUN, binN) {
   # overlap gene bins and SNPs
   oiX <- data.table::copy(OI)[, .(chr, objID, startpos, endpos)] # avoid overwriting obj.info
   oiX <- oiX[VI, on = .(chr = chr, startpos <= pos, endpos >= pos), nomatch = 0]

   # get max value and record no. overlapping SNPs for each gene
   oiX <- oiX[, .(objMax = max(value), N.ovlp = .N), keyby = objID]

   if (FUN == "lm.logN") { # residuals from linear regression on log-transfromed SNP / bin count
      lm.fit <- oiX[, stats::lm(objMax ~ log(N.ovlp))]
      oiX[, objMax.res := stats::resid(lm.fit)]
      mPars <- lm.fit$coefficients
   } else { # residuals from non-paramtetic correction
      # bins genes according to no. of SNPs / bins
      ob0 <- oiX[, .(n = .N), keyby = N.ovlp] # tabulate genes by SNP / bin count
      fm <- ob0[1]$N.ovlp # start index
      fM <- ob0[.N]$N.ovlp
      n.min <- fm
      while (fm < fM) {
         ob0[N.ovlp < fm, csN := 0]
         ob0[N.ovlp >= fm, csN := cumsum(n)]
         if (any(ob0$csN >= binN)) { # combining following bins does exceed min value
            fm <- ob0[which.max(csN >= binN)]$N.ovlp + 1
            n.min <- c(n.min, fm)
         } else { # final bins do not exceed min value, merge all remaining bins
            fm <- fM
         }
      }

      nB <- length(n.min)
      oBins <- data.table::data.table(objBin = 1:nB, n.min = n.min,
                                      n.max = c(n.min[2:nB] - 1, fM))

      # assign genes to bins and compute median-based modified z score
      oiX[oBins, on = .(N.ovlp >= n.min, N.ovlp <= n.max), objBin := objBin]
      oiX[, objMed := stats::median(objMax), by = objBin]
      oiX[, objMed.res := objMax - objMed]
      oiX[, MAD := stats::median(abs(objMed.res)), by = objBin]
      oiX[, objMax.res := 0.6475 * objMed.res / MAD]
      mPars <- merge(oBins,
                     oiX[, .(.N, Med = objMed[1], MAD = MAD[1]), by = objBin],
                     by = "objBin") # keep bins and relevant parameter information
   }

   mp <- oiX[, .(objID, objMax, N.ovlp)] # retain summary
   attr(mp, "fit.params") <- mPars # add model information as attribute
   return(list(oiX, mp))
}


#' @title Convert physical distances to genetic distances
#' @description Internal function that converts base pair (bp) positions to
#'  centimorgan (cM) genetic distances. It uses a recombination rate map and a
#'  specified map function (e.g., Kosambi or Haldane) to perform the conversion.
#' @param OI data.table; The gene information table containing physical
#'  gene positions.
#' @param pos character; vector with position columns to convert
#' @param RR data.table; The recombination rate dataset (in hapmap format).
#' @param map.fun character; The genetic map function to use.
#'  Options are "kosambi", "haldane", or `NULL` in cases where map column is
#'  provided.
#' @import data.table
#' @return No return value. Instead an updated `OI` data.table with specified
#'  columns converted to genetic distances.
#' @noRd
.bp2cM <- function(OI, pos, RR, map.fun) {
   if (!is.null(map.fun)) {    # add genetic distance column if not present
      # get appropriate map function
      if (map.fun == "haldane") mf <- function(d) 0.5 * exp(-d / 50)
      if (map.fun == "kosambi") mf <- function(d) 0.5 * tanh(d / 50)
      # NB: rescale by 1e-4 = 100 (centiMorgans) / 1e6 (mega-basepairs)
      RR[, cM := cumsum(c(0, mf(d = rate[-.N]) * diff(pos))) / 1e4]
   } else {
      g.dist <- c("cM", "cm", "map", "Map") %in% names(RR)
      gd0 <-c("cM", "cm", "map", "Map")[g.dist] # determine user labels
      if (!"cM" %in% gd0) { # force canonical label use
         data.table::setnames(RR, old = gd0[1], new = "cM")
      }
   }

   RR <- RR[, .(SP.bp = pos[-.N], EP.bp = pos[-1],
                SP.cM = cM[-.N], EP.cM = cM[-1]), by = chr]
   RR[, bp2cM := (EP.cM - SP.cM) / (EP.bp - SP.bp)]
   data.table::setkey(RR, chr, SP.bp, EP.bp)

   # programmatically generate new genetic distance columns
   for (posI in pos) { # loop over all positions
      data.table::setnames(OI, old = posI, new = "p0")
      OI[.(RR), on = .(chr = chr, p0 >= SP.bp, p0 < EP.bp),
         (posI) := SP.cM + (p0 - SP.bp) * bp2cM]
      OI[, p0 := NULL]
   }
}


#' @title Consolidate merged gene sets
#' @description Internal function that resolves multiple merges of gene sets.
#'  Takes two intermediate merge results and produces a consolidated mapping of
#'  the original set IDs to their final merged set ID.
#' @param I data.table; An intermediate merge result table.
#' @param J data.table; Another intermediate merge result table.
#' @import data.table
#' @return A data.table with two columns, `setID` (the original set ID) and
#'  `setID.new` (the final set ID).
#' @noRd
.final_merge <- function(I, J) {
   # find sets that have been merged multiple times
   mgm <- merge(I, J, by.x = "setID.new", by.y = "setID", all.y = TRUE)
   K <- mgm[, .(setID = ifelse(is.na(setID), setID.new, setID),
                TEMP = setID.new.y)]
   out <- merge(I, K, by = "setID", all.x = TRUE)
   return(out[, .(setID, setID.new = ifelse(is.na(TEMP), setID.new, TEMP))])
}


#' @title Identify subnetworks of linked gene sets
#' @description Internal function that identifies groups of gene sets with
#'  shared genes above a specified similarity threshold. Models the gene
#'  sets as a graph and uses connected components to define subnetworks.
#' @param SO data.table; The gene to gene set mapping data.table.
#' @param min.sim numeric. The minimum proportion of shared genes required
#'  for unification into a single gene set.
#' @import data.table
#' @importFrom igraph components graph_from_data_frame groups
#' @return An `integer` vector where each value represents the subnetwork ID for
#'  each gene set.
#' @noRd
.get_subnetworks <- function(SO, min.sim) {
   # no. of shared genes between sets
   SxG <- SO[SO, on = .(objID = objID),
             allow.cartesian = TRUE][, .N, keyby = .(setID, i.setID)]
   sN <- SO[, .(setN = .N), keyby = setID]
   # find cases where |set A intersection set B | / |set A| exceeds min.sim
   SxG <- SxG[sN, on = .(setID = setID)][N / setN >= min.sim, .(setID, i.setID)]
   # find cases where this is reciprocated for both set A and set B
   SxG <- merge(SxG, SxG, by.x = c("setID", "i.setID"),
                by.y = c("i.setID", "setID"))
   igg <- igraph::graph_from_data_frame(SxG, directed = FALSE) # convert to graph
   sgg <- igraph::groups(igraph::components(igg)) # find connected subgraphs
   out <- data.table::data.table(setID = unlist(sgg),
                                 setG = rep(1:length(sgg),
                                            times = sapply(sgg, length)))
   if (class(SO$setID) != class(out$setID)) { # enforce class similarity
      out[, setID := as(setID, Class = class(SO$setID))]
   }
   return(out)
}


#' @title Merge gene sets based on a similarity threshold
#' @description Internal function that iteratively merges gene sets sharing a
#'  proportion of genes greater than a specified minimum similarity (`min.sim`).
#'  It uses a subnetwork-based approach to identify and consolidate sets,
#'  ensuring that the merging process is consistent and complete.
#' @param SI data.table; The gene set information data.table.
#' @param SO data.table; The gene to gene set mapping data.table.
#' @param min.sim numeric; The minimum proportion of shared genes required
#'  for unification into a single gene set.
#' @param ENV environment; The environment where the merged data tables will
#'  be stored.
#' @import data.table
#' @return No value is returned. Instead, the function modifies the
#'  provided environment (`ENV`) by adding the updated `set.info` and `set.obj`
#'  tables, as well as generating an object tracking the merged sets.
#' @noRd
.merge_sim_sets <- function(SI, SO, min.sim, ENV) {
   SI[, setN := SO[, .N, keyby = setID]$N] # get gene set size

   # identify subnetworks of gene sets with shared genes
   SI0 <- .get_subnetworks(SO = SO, min.sim = min.sim)
   merged.sets <- list()
   R <- 0

   # repeat until all gene sets fulfill gene sharing requirements
   while (max(SI0$setG) < nrow(SI)) {
      R <- R + 1
      SI[SI0, on = .(setID = setID), setG := setG] # update setIDs

      # rename merged gene sets using name from set with most genes
      SI[, `:=` (Nm = .N, setID.new = setID[which.max(setN)]), by = setG]
      merged.sets[[R]] <- SI[Nm > 1, .(setID.new, setID)][order(setID.new)]

      # update files
      SO[SI, on = .(setID = setID), setID := setID.new] # update set names
      SO <- SO[, .(objID = unique(objID)), by = setID] # ensure only one copy of each gene per set

      SI <- SI[setID == setID.new] # remove redundant gene sets
      SI[, `:=` (setG = NULL, Nm = NULL, setID.new = NULL)] # remove unused columns
      SI[, setN := SO[, .N, keyby = setID]$N] # update set sizes

      # check that new sets fulfill merging requirements
      SI0 <- .get_subnetworks(SO = SO, min.sim = min.sim)
   }

   # prepare output
   if (R == 0) {
      SO.xx <- NULL
   } else { # compile final merged sets
      SO.xx <- Reduce(.final_merge, merged.sets)
      if ("setName" %in% names(SI)) {
         SI[SO.xx, on = .(setID = setID.new), setName := paste0(setName, "!")]
      }
   }

   list2env(list(set.info = SI, set.obj = SO, set.info.merged = SO.xx, R = R),
            envir = ENV)
}


#' @title Convert Mahalanobis distances to probability weights
#' @description Internal function that converts a matrix of Mahalanobis
#'  distances into a matrix of probability weights. The conversion is based on
#'  a specified distance decay function.
#' @param wt matrix; A matrix of Mahalanobis distances, where each row is a
#'   vector containing inter-gene distances relative to a focal gene.
#' @param kern.scale numeric; The scaling parameter for the kernel function.
#' @param kern.func character; The kernel function to use:
#'  Options are 'inverse', 'gaussian', or 'exponential'.
#' @importFrom Rfast rowsums
#' @return A matrix of probability weights.
#' @noRd
.md2p <- function(wt, kern.scale, kern.func) {
   if (kern.func == "gaussian") { # Gaussian weights
      wt <- exp(-kern.scale * wt ^ 2)
   } else {
      if (kern.func == "exponential") { # Exponential weights
         wt <- exp(-kern.scale * wt)
      } else { # inverse weights
         wt <- (wt + 1e-6) ^ -kern.scale
      }
   }
   wt <- wt / Rfast::rowsums(wt) # convert weights to probabilities
   return(wt)
}


#' @title Cap and re-normalize probability weights
#' @description Internal function that converts weight matrix to probabilities
#'  and ensures no individual probability exceeds a specified maximum value
#'  (maxP). Iteratively caps and rescales the remaining probabilities in
#'  each row until all values are at or below the maximum, while maintaining the
#'  sum of 1.
#' @param mm0 matrix; The probability matrix to be capped.
#' @param maxP numeric; The maximum probability weight allowed.
#' @importFrom Rfast rowsums
#' @return A matrix where all probabilities are at or below `maxP` and each
#'  row sums to 1.
#' @noRd
.cap_probs <- function(mm0, maxP) {
   wS <- Rfast::rowsums(mm0) # get row weight sums
   wMax <- maxP * wS # maximum weight in each row
   wGT0 <- which(mm0 > wMax, arr.ind = TRUE) # check that capping is required
   if (nrow(wGT0) > 0) {
      row0 <- sort(unique(wGT0[, 1])) # rows requiring capping
      wGT <- wGT0 # track capped values
      wS0 <- wS # track changing weights
      while (length(row0) > 0) { # recursively update probabilities until none > cap
         mm0[wGT0] <- 0 # temporarily set capped weigths to 0
         rv <- table(wGT0[, 1]) * wS[row0] * maxP # revised weight at capped values
         wS0[row0] <- wS0[row0] - rv # revised row weight to be redistributed to uncapped values
         wS.ac <- Rfast::rowsums(mm0[row0, , drop = FALSE]) # actual weights in uncapped values
         mm0[row0, ] <- mm0[row0, ] * wS0[row0] / wS.ac # rescale to preserve weight sums

         # check revised values for those exceeding cap and track values
         wGTi <- which(mm0 > wMax, arr.ind = TRUE)
         wGT <- rbind(wGT, wGTi)
         wGT0 <- wGTi
         row0 <- sort(unique(wGT0[, 1]))
      }
      mm0[wGT] <- wMax[wGT[, 1]] # insert capped values
   }
   mm0 <- mm0 / wS # convert to probability weights
   return(mm0)
}


#' @title Fit local quadratic regression to centred confounder covariates
#' @description Internal function that estimates the coefficients for the local
#'  quadratic regression. Initially fits a full quadratic model (a combination
#'  of constant linear, and quadratic terms), then successively fits a linear
#'  and constant model if model fitting fails.
#' @param d1 numeric matrix; The centred covariates used in the first order
#'  (linear) regression terms.
#' @param wt numeric vector: The weights used in the local regression.
#' @param os0 numeric vector; The gene scores.
#' @param n.cv numeric vector; The number of confounder covariates.
#' @import foreach
#' @return A numeric vector with model coefficients.
#' @noRd
.fit_local_quad_reg <- function(d1, wt0,  object_scores, n.cv) {
   # get second degree terms
   d2 <- foreach::foreach(j = 1:n.cv) %do% {
      d1[, j:n.cv] * d1[, j]
   }
   d2 <- cbind(1, d1, do.call(cbind, d2)) # create design matrix
   dX <- d2 * wt0 # weighted design matrix
   LHS <- crossprod(dX, d2) # left-hand side (LHS) of the normal equation: (X'WX)
   RHS <- crossprod(dX, object_scores) # right-hand side (RHS) of the normal equation: (X'WY)
   # attempt to solve 2nd order polynomial equation
   theta0 <- tryCatch(solve(LHS, RHS), error = function(e) e)
   if (methods::is(theta0, "error")) { # try linear model
      ri <- 1:(n.cv + 1)
      theta0 <- tryCatch(solve(LHS[ri], RHS[ri, ri]), error = function(e) e)
      if (methods::is(theta0, "error")) { # constant model
         theta0 <- c(RHS[1], rep(NA, ncol(d1) - 1))
      } else {
         theta0 <- c(theta0, rep(NA, ncol(d1) - n.cv - 1))
      }
   }
   return(theta0)
}


#' @title Combine null model outputs
#' @description Internal function designed to iteratively merge the
#'  outputs from a `foreach` loop used for generating a null distribution and
#'  GPD tail fits.
#' @param A list or matrix; The first output from a null generation iteration.
#' @param B list or matrix; The second output to be merged with `A`.
#' @return A list containing the combined object from the iteration.
#'  Includes the ECDF `e0` and GDP fits `g0` for evaluated gene set sizes.
#' @noRd
.combine_results <- function(A, B) {
   # iteratively combine raw value lists
   if (is.null(A$ecdf_values$qX)) { # list form in parallel inner loop
      A$ecdf_values <- Map(function(x, y) c(x, list(y)), A$ecdf_values, B$ecdf_values) # aggregate raw value lists
      # aggregate tails
      AB <- cbind(A$gpd_params, B$gpd_params)
      nT <- ncol(AB) / 2
      rTH <- Rfast::rowMedians(AB)
      rID <- which(AB >= rTH, arr.ind = TRUE)
      tt <- Rfast::rowsums(AB >= rTH) # check for more than nT entries
      if (any(tt > nT)) {
         wTT <- which(tt > nT)
         for (wT in wTT) {
            xOUT <- which(AB[wT, ] == rTH[wT])[1:(tt[wT] - nT)] # identify values to remove
            rID <- rID[!(rID[, 1] == wT & rID[, 2] %in% xOUT), ] # remove excess rows
         }
      }
      A$gpd_params <- matrix(AB[rID[order(rID[, 1]), ]], byrow = TRUE, ncol = nT)
   } else { # raw form in sequentail outer loop (bootstraps)
      A$ecdf_values$qX <- Map(rbind, A$ecdf_values$qX, B$ecdf_values$qX) # combine ECDFs
      A$ecdf_values$qE <- Map(rbind, A$ecdf_values$qE, B$ecdf_values$qE) # combine ECDF errors
      A$gpd_params <- rbind(A$gpd_params, B$gpd_params) # bind tails
   }
   if (!is.null(A$rescaling_factors)) { # collect means for set-wise rescaling
      A$rescaling_factors <- A$rescaling_factors + B$rescaling_factors
   }
   return(A)
}


#' @title Estimate Generalized Pareto Distribution (GPD) parameters
#' @description Internal function that fits a Generalized Pareto Distribution
#'  to the tail of a null distribution. Facilitates the calculation of p-values
#'  below the minimum resolution achieved by the permutations.
#' @param x0 numeric; A vector of null exceedences.
#' @importFrom ismev gpd.fit
#' @return A list containing the estimated GPD parameters and associated
#'  measurement errors.
#' @noRd
.fit_gpd_tail <- function(x0) {
   gpd.fit0 <- tryCatch(
      ismev::gpd.fit(xdat = x0, threshold = min(x0), show = FALSE),
      error = function(e) e,
      warning = function(w) w
   )
   if (methods::is(gpd.fit0, "error") | methods::is(gpd.fit0, "warning")) {
      return(rep(NA, 4))
   } else { # gpd fitting successful
      return(c(gpd.fit0$mle, gpd.fit0$se ^ 2))
   }
}


#' @title Merge raw values and extract empirical quantiles
#' @description Internal function that combines a list of raw numeric vectors
#'  and returns requested quantiles with an error band estimated by a small
#'  symmetric probability perturbation.
#' @param tX list; List of numeric vectors to combine.
#' @param eQnt numeric; Vector of probabilities at which to estimate quantiles.
#' @param eps numeric; Small perturbation for error estimation.
#' @return A list containing a numeric vector of estimated quantiles (`qX`) and
#'  a numeric vector of error widths computed at `eps` distances either side of
#'  each quantile.
#' @noRd
.get_quantiles <- function (tX, empirical_quantiles, eps) {
   # Combine all raw values into a single vector
   all_vals <- unlist(tX, use.names = FALSE)
   # Extract quantiles using base R quantile function
   qX <- stats::quantile(all_vals, probs = empirical_quantiles, names = FALSE)
   # quantify error bounds
   qL <- stats::quantile(all_vals, probs = empirical_quantiles - eps, names = FALSE)
   qH <- stats::quantile(all_vals, probs = empirical_quantiles + eps, names = FALSE)
   return(list(qX = qX, qE = qH - qL))
}


#' @title Compute minimum achievable p-values and z cutoffs from GPD parameters
#' @description Internal function. Given a GPD fit (`gI`) and a cutoff factor
#'  (`gpd.cutoff`), returns the minimum attainable p-value in the tail region,
#'  the corresponding z threshold, and an adjustment factor to ensure monotonic
#'  tail probabilities.
#' @param gI numeric; A one-row object containing GPD parameters.
#' (`shape`, `scale`) and associated quantities used in GPD p value calculation.
#' @param gpd.cutoff numeric; Proportion of null used in GPD estimation.
#' @return A numeric matrix with adjustment factor applied to tail probabilities
#'  for monotonicity (`adj.p`), the value at which the minimum tail probability
#'  occurs (min.z), and the minimum possible tail probability (min.p).
#' @noRd
.get_gpd_minimums <- function(gI, gpd.cutoff, q.bnd) {
   gz <- (1:1e6L / 1e4L) # assume gpd scaled values
   if (gI$shape == 0) {
      pU <- exp(-gz)
   } else {
      pU <- suppressWarnings(with(gI, exp((-1 / shape) * log1p(shape * gz))))
   }
   # minimum possible p-value
   pMin <- min(pU[pU > 0], na.rm = TRUE)
   min.p <- pMin * gpd.cutoff

   # z value at minimum
   zMin <- which.max(pU == pMin)
   min.z = gI$exc.z + gz[zMin] * gI$scale

   # determine adjustment factor to ensure monotonic probabilities
   diff.z <- gI[, (cut.z - exc.z) / scale]
   if (gI$shape == 0) {
      pU <- exp(-diff.z)
   } else {
      pU <- with(gI, exp((-1 / shape) * log1p(shape * diff.z)))
   }
   adj.p <- q.bnd / pU # (q.bnd * gpd.cutoff) / (pU * gpd.cutoff)
   return(cbind(adj.p = adj.p, min.z = min.z, min.p = min.p))
}


#' @title Calculate the lower-tail p-values combining ECDF and GPD tail
#' @description Internal function that calculates p-values for gene set scores
#'  using an empirical cumulative distribution function (ECDF) for the main body
#'  of the distribution, switching to a fitted Generalized Pareto Distribution
#'  (GPD) model for the extreme lower tail values.
#' @param ss numeric; ector of gene set scores.
#' @param n integer; Index selecting the set-specific ECDF and GPD parameters.
#' @param g0 data.table; Contains GPD parameters (`shape`, `scale`) and
#'  associated quantities used in GPD p value calculation for each gene set size.
#' @param e0 list; ECDF functions for each gene set size.
#' @import data.table
#' @return A numeric vector of p-values, corresponding to the input scores.
#' @noRd
.get_pvalue_lower <- function(ss, n, gpd_params, ecdf_values) {
   pU <- ecdf_values[[n]](ss) # get ECDF quantile
   gI <- gpd_params[n]
   if (any(ss <= gI$cut.z)) { # check for set scores in GPD tail region
      wp <- which(ss <= gI$cut.z)
      if (any(ss <= gI$min.z)) { # some values below absolute minimum p
         wp0 <- which(ss <= gI$min.z)
         pU[wp0] <- gI$min.p
         wp <- wp[!(wp %in% wp0)] # update wp
      }
      q.gpd <- with(gI, (ss[wp] - exc.z) / scale) # place on gpd scale
      if (gI$shape == 0) {
         pU[wp] <- exp(-q.gpd) * gI$adj.p
      } else {
         pU[wp] <- with(gI, exp((-1 / shape) * log1p(shape * q.gpd)) * adj.p)
      }
   }
   return(pU)
}


#' @title Calculate the upper-tail p-values combining ECDF and GPD tail
#' @description Internal function that calculates p-values for gene set scores
#'  using an empirical cumulative distribution function (ECDF) for the main body
#'  of the distribution, switching to a fitted Generalized Pareto Distribution
#'  (GPD) model for the extreme upper tail values.
#' @param ss numeric; ector of gene set scores.
#' @param n integer; Index selecting the set-specific ECDF and GPD parameters.
#' @param g0 data.table; Contains GPD parameters (`shape`, `scale`) and
#'  associated quantities used in GPD p value calculation for each gene set size.
#' @param e0 list; ECDF functions for each gene set size.
#' @import data.table
#' @return A numeric vector of p-values, corresponding to the input scores.
#' @noRd
.get_pvalue_upper <- function(ss, n, gpd_params, ecdf_values) {
   pU <- 1 - ecdf_values[[n]](ss) # get ECDF quantile
   gI <- gpd_params[n]
   if (any(ss >= gI$cut.z)) {
      wp <- which(ss >= gI$cut.z)
      if (any(ss >= gI$min.z)) { # some values below absolute minimum p
         wp0 <- which(ss >= gI$min.z)
         pU[wp0] <- gI$min.p
         wp <- wp[!(wp %in% wp0)] # update wp
      }
      q.gpd <- with(gI, (ss[wp] - exc.z) / scale) # place on gpd scale
      if (gI$shape == 0) {
         pU[wp] <- exp(-q.gpd) * gI$adj.p
      } else {
         pU[wp] <- with(gI, exp((-1 / shape) * log1p(shape * q.gpd)) * adj.p)
      }
   }
   return(pU)
}


#' @title Expected set-size covariance under finite-population cumulative sums
#' @description Constructs the Brownian-bridge-like correlation induced by
#'  cumulative sums without replacement on a finite population of size `n.genes`,
#'  evaluated at specific set sizes (`x`). Diagonal is set to 1 and symmetry
#'  enforced explicitly.
#' @param x numeric; Vector of assessed set sizes (strictly between 0 and
#' `n.genes`).
#' @param n.genes integer; Number of genes used in permutations.
#' @return Numeric matrix `K` of size length(x) by length(x) with unit diagonal.
#' @noRd
.est_ss_cov <- function(x, n.genes) {
   denom <- tcrossprod(sqrt(x * (n.genes - x)))
   num_min <- outer(x, x, pmin) # min(i,j)
   num_max <- outer(x, x, pmax) # max(i,j)
   num <- num_min * (n.genes - num_max) # numerator
   K <- num / denom
   diag(K) <- 1 # enforce unit diagonal
   K <- (K + t(K)) * 0.5 # enforce symmetry
   return(K)
}


#' @title Construct whitened least-squares smooth for mgcv
#' @description Internal constructor for a custom `bs` = `wls` smooth that fits
#'  in a whitened equation space. Consumes, via `xt`, the precomputed whitened
#'  design `xW` and the original inner smooth `smI`. Reparameterisation and
#'  centering are disabled to preserve GLS geometry. Penalty matrices and
#'  scaling are carried explicitly.
#' @param object list; Smooth specification (from mgcv).
#' @param data model frame; Unused; required by signature.
#' @param knots numeric; Unused; required by signature.
#' @method smooth.construct wls.smooth.spec
#' @export
#' @importFrom mgcv smooth.construct
#' @return A smooth object with class `c("wls.smooth", class(xt$smI))`
#' @details
#'  The following elements are expected in `object$xt`:
#' \itemize{
#'   \item `xW`: whitened basis matrix (rows pre-multiplied by whitener).
#'   \item `smI`: inner smooth from `mgcv::smoothCon`.
#'   \item `S.sW`: penalty scaling pre-aligned to `xW`.
#' }
#' The constructor sets `sm$X <- xW`, disables additional reparameterisation,
#' and assigns penalties `D`, `S`, and `S.scale` to match `xW`.
#' @noRd
smooth.construct.wls.smooth.spec <- function(object, data, knots) {
   xt <- object$xt
   sm <- xt$smI # unwhitened inner smooth from smoothCon
   sm$X <- xt$xW # use the whitened design for fitting
   sm$smI <- xt$smI # keep the original inner smooth for prediction
   sm$side.constrain <- FALSE # do not impose centering
   sm$repara <- FALSE # do not rotate basis
   sm$null.space.dim <- 0 # tell mgcv there is no null space to remove
   sm$D <- xt$smI$D # carry the original difference-penalty structure
   sm$S <- xt$smI$S # carry penalty matrices
   sm$S.scale <- xt$S.sW # control rescaling
   sm$sp <- xt$smI$sp # carry starting sp
   attr(sm,  "qrc") <- NULL # drop qrc only on the constructed smooth to avoid warning
   attr(sm$X, "qrc") <- NULL # drop any qrc attr on design
   class(sm) <- c("wls.smooth", class(xt$smI))
   return(sm)
}


#' @title Predict matrix for whitened least-squares smooth (internal)
#' @description Internal function. Prediction method for the custom `bs` = `wls`
#'  smooth. Delegates to the next class (e.g., `pspline.smooth`) so that the
#'  prediction basis matches the parameterisation used at fit time.
#' @param object smooth object; Construct with class `wls.smooth`.
#' @param data data.frame; Contains the predictor variable(s).
#' @method Predict.matrix wls.smooth
#' @export
#' @importFrom mgcv Predict.matrix
#' @return A prediction matrix compatible with the fitted coefficients.
#' @noRd
Predict.matrix.wls.smooth <- function(object, data) {
   cls <- class(object)
   class(object) <- setdiff(cls, "wls.smooth")
   mgcv::Predict.matrix(object, data)
}


#' @title Covariance-aware penalized GLS smoother
#' @description Internal function that fits a Gaussian additive model in a
#'  whitened equation space to account for known process covariance across set
#'  sizes. The function builds a covariance matrix that combines the known
#'  autocovariance (`K`) structure with the measurment error (`tau`) from GPD
#'  parameter or quantile estimation. Computes a pivoted Cholesky whitener to
#'  whiten the response, the rows of the inner basis, and a parametric intercept
#'  column. Rescales the penalty strength to the whitened design, and calls
#'  `mgcv::gam` with a custom smooth that uses the whitened basis at fit time
#'  while predicting on the original set-size scale.
#' @param y numeric; Vector of parameter/quantile estimates at assessed set
#'  sizes (averaged over bootstraps).
#' @param x numeric; Vector of assessed set sizes.
#' @param sm0 smooth object; Inner smooth returned by `mgcv::smoothCon` built
#'  on `x`. Must provide `sm0$X` (unwhitened design), `sm0$bs.dim` (basis
#'  dimension), and penalties `sm0$S`, `sm0$D`, `sm0$S.scale`.
#' @param tau numeric; Vector of measurement variances for `y`. (averaged over
#'  bootstraps).
#' @param K numeric; Covariance matrix across set sizes that captures the
#'  correlation created by cumulative summing and FPC-standardisation.
#' @param n.th integer; Number of assessed set sizes.
#' @return A fitted `mgcv::gam` object.
#' @importFrom mgcv gam s
#' @details
#'  This function performs GLS by pre-whitening the equation:
#'  \deqn{
#'  W y = (W X)\beta + (W \mathbf{1})\alpha + \varepsilon, \quad \varepsilon \sim N(0, I),
#'  }
#'  then fits with `mgcv::gam()` using a custom smooth `bs` = `wls` that consumes:
#'  \itemize{
#'   \item the whitened basis `xW` (rows of `sm0$X` left-multiplied by \eqn{W}),
#'   \item the original inner smooth `smI` = `sm0` (for prediction on the natural
#'   set-size scale),
#'   \item a penalty scale `S.sW` adjusted to the whitened design via
#'   \code{sm0$S.scale * sum(xW^2) / sum(sm0$X^2)}.
#'  }
#' @noRd
.parallel_smooth <- function(y, x, sm0, tau, K, n.th) {
   Sigma <- K + diag(tau + rep(1e-9, n.th)) # include autocovariance and measurement error
   R <- chol(Sigma, pivot = TRUE) # Cholesky whitening matrix
   piv <- attr(R, "pivot") # allow column pivoting
   pivU <- order(piv) # undo pivot
   yW <- backsolve(R, y[piv], transpose = TRUE)[pivU] # whiten response
   xW <- backsolve(R, sm0$X[piv, ], transpose = TRUE)[pivU, ] # whiten basis from initial smooth
   oneW <- backsolve(R, rep(1, n.th), transpose = TRUE)[pivU] # whiten intercept
   S.sW <- sm0$S.scale * (sum(xW ^ 2) / sum(sm0$X ^ 2)) # rescale
   # Fit GAM with custom whitened smooth: pass xW and S.sW for fitting and sm0 for prediction
   fit <- mgcv::gam(yW ~ 0 + oneW + s(x, bs = "wls", k = sm0$bs.dim,
                                      xt = list(xW = xW, smI = sm0,
                                                S.sW = S.sW)),
                    data = data.frame(x = x, yW = yW), method = "REML")
   return(fit)
}


#' @title Construct exponentially spaced lag bins for inter-gene autocovariance
#' @description Internal utility function to construct optimally sized, exponentially spaced
#'  genomic lag bins for summarising inter-gene autocovariance. The binning scheme
#'  prioritises resolution at small genomic distances while progressively
#'  coarsening at larger lags. Bins are adaptively merged to ensure a minimum
#'  number of gene pairs per bin on each chromosome.
#' @param maxL numeric; Maximum genomic lag (in base pairs) to be considered
#'  when constructing bins.
#' @param minL numeric; Minimum non-zero genomic lag used to partition 0-lag
#'  covariance values.
#' @param cgm.bin integer; Minimum required number of inter-gene pairs per bin
#'  and chromosome. Bins containing fewer observations are sequentially merged
#'  with neighbouring bins until this threshold is met.
#' @param OI list; List of data objects, one per chromosome, each containing at
#'  least the chromosome identifier (`chr`) and a vector of gene midpoints
#'  (`pM`).
#' @return Numeric vector giving the final lag-bin boundaries (in base pairs).
#'  The vector includes the lower boundary of each retained bin and the upper
#'  boundary of the final bin.
#' @import data.table
#' @details
#'  Initial bin boundaries are constructed on an exponential scale between
#'  `minL` and `maxL`, yielding between 100 and 1000 candidate bins depending on
#'  the genomic range. Pairwise inter-gene distances are then binned by
#'  chromosome, and adjacent bins are merged iteratively from small to large
#'  lags to satisfy the minimum pair-count requirement (`cgm.bin`) across all
#'  chromosomes. This adaptive binning strategy concentrates information at
#'  small genomic distances—where autocovariance structure is strongest—while
#'  maintaining statistical stability at larger lags.
#' @noRd
.get_cgm_bins <- function (maxL, minL, cgm.bin, OI) {
   n.chr <- length(OI)
   min.log <- log10(1e5 * minL) # minimum possible lag range
   max.log <- log10(5e7 * minL) # maximum possible lag range
   nBins <- floor(100 + 900 * (log10(maxL) - min.log) / (max.log - min.log)) # range between 100 to 1000 bins
   p <- log(minL / maxL) / log(1 / nBins) # calculate required power exponent
   cBins <- maxL * (0:nBins  / nBins) ^ p # calculate initial bin boundaries

   # determine lag-bin counts
   cN <- lapply(OI, function(x) {
      xD <- dist(x$midpos)
      xD <- xD[which(xD <= maxL)] # values within range

      # extract and bin covariances
      bins <- cut(xD, cBins, labels = 1:nBins, right = FALSE)

      data.table::data.table(chr = x[1]$chr, bin = 1:nBins,
                             N = unname(c(table(bins))))
   })
   cN <- do.call(rbind, cN) # merge chromosome results into single table

   # sequentially aggregate bins to meet minimum size requirement
   bin0 <- 1:nBins
   if (nrow(cN[bin != 1L & N < cgm.bin]) > 0L) { # at least one bin below required size
      fm <- 2 # start index: ignore 0-lag bins
      while (fm < nBins) {
         csN <- cN[bin %in% fm:nBins, .(csN = cumsum(N)), by = chr] # cumulative sum
         if (any(csN >= cgm.bin)) { # combining following bins does exceed min value
            mB <- max(csN[, which.max(csN >= cgm.bin), by = chr]$V1)
            fm0 <- fm + mB - 1
         } else { # final bins do not exceed min value, merge all remaining bins
            fm0 <- nBins
         }
         bin0[fm:fm0] <- fm
         fm <- fm0 + 1
      }
   }
   bin0 <- unique(bin0)

   return(c(cBins[bin0], tail(cBins, 1)))
}


#' @title Exponential LD covariance with overlap correction
#' @description Internal covariance function combining an exponential
#'  recombination-based decay with a locally acting overlap correction. The LD
#'  component models background linkage-disequilibrium-driven covariance as a
#'  function of genetic distance, while the overlap term captures short-range
#'  excess covariance induced by overlapping or embedded gene intervals.
#' @param h numeric; Vector of genetic lag distances (e.g. centiMorgans)
#'  at which the covariance is evaluated.
#' @param ovlp numeric; Vector of gene-pair overlap measures derived from
#'  physical (base-pair) gene coordinates.
#' @param beta1 numeric; Marginal variance (amplitude) of the LD component.
#' @param beta2 numeric; Strength of the overlap-induced covariance component.
#' @param scale numeric; Exponential decay scale (range) governing LD decay
#'  in genetic distance units.
#' @return Numeric vector of model-implied covariances evaluated at the
#'  supplied lags and overlap values.
#' @details
#'  The covariance is defined as:
#'  \deqn{
#'    C(h) = \beta_1 \exp(-h / r_1) + \beta_2 \, \mathrm{ovlp} \, \exp(-h / r_2),
#'  }
#'  where the first term captures genome-wide LD-mediated dependence as a
#'  function of recombination distance, and the second term provides a
#'  locally acting adjustment for gene overlap effects. The overlap term
#'  is gated by the overlap metric and decays rapidly with lag, ensuring
#'  that annotation-driven covariance does not distort long-range LD
#'  inference.
#' @noRd
.covariance_function <- function(h, ovlp, beta1, beta2, scale) {
   ld_cov   <- beta1 * exp(-h / scale) # LD-driven covariance
   ovlp_cov <- beta2 * ovlp # Local overlap correction
   return(ld_cov + ovlp_cov)
}


#' @title Composite likelihood fit of exponential covariance model
#' @description Internal objective function for fitting a two-component
#'  exponential covariance model to binned empirical covariances of gene
#'  scores using a heteroskedastic Gaussian quasi–composite likelihood.
#'  The model combines an LD-driven exponential decay in genetic distance
#'  with a locally acting overlap correction that captures excess covariance
#'  induced by overlapping or embedded gene intervals.
#' @param par numeric; Vector of model parameters to esitmate. Elements include:
#'  \describe{
#'    \item{beta1}{Amplitude of the LD-driven covariance component.}
#'    \item{beta2}{Strength of the overlap-induced covariance component.}
#'    \item{scale}{Exponential decay scale for LD effects.}
#'  }
#' @param fixed numeric; Global scale parameter for chromosome models
#' @param cov_emp numeric; Vector of empirical covariances for each lag bin.
#' @param se numeric; Vector of standard errors associated with `cov_emp`,
#'  typically derived from robust covariance estimators.
#' @param ovlp numeric; Vector of overlap values corresponding to each bin.
#' @param h numeric; Vector of lag bin midpoints (genetic distance units).
#' @param w numeric; Vector of non-negative bin weights, typically proportional
#'  to pair density.
#' @return Numeric scalar giving the weighted negative composite
#'  log-likelihood.
#' @details
#'  The loss corresponds to a heteroskedastic Gaussian quasi-likelihood:
#'  \deqn{
#'    \hat C(h) \sim N(C_\theta(h), \sigma_h^2),
#'  }
#'  where \eqn{C_\theta(h)} is the model-implied covariance and \eqn{\sigma_h^2}
#'  is supplied via `se^2`. Bin weights scale the contribution of each lag bin,
#'  reflecting differing information content. The overlap component provides
#'  targeted correction at short lags, preventing distortion of the inferred LD
#'  decay scale.
#' @noRd
.fit_cgm_function <- function(par, fixed = NULL, cov_emp, se, ovlp, h, w) {
   # Extract parameters
   b1 <- par[1]
   b2 <- par[2]
   lam <- ifelse(is.null(fixed), par[3], fixed)

   # Predicted covariance: Matérn + linear gene interval terms
   cov_pred <- .covariance_function(h = h, ovlp = ovlp, beta1 = b1, beta2 = b2, scale = lam)
   cov_pred <- pmax(cov_pred, 1e-12) # Numerical stabilisation

   # Heteroskedastic variances from Huber-M SEs
   se2 <- se ^ 2
   var.emp <- pmax(se2, quantile(se2, 0.05)) # protect against tiny error margins

   # Composite likelihood contribution
   loss <- (cov_emp - cov_pred) ^ 2 / var.emp + log(var.emp)

   # Weighted composite likelihood
   sum(w * loss)
}


#' @title Estimate chromosome-level covariance parameters via empirical Bayes
#' @description Internal helper function to estimate chromosome-level covariance
#'  amplitude parameters using a Gaussian random-effects model fitted by REML.
#'  Chromosome-specific LD and overlap amplitudes are treated as random effects
#'  around genome-wide means, providing partial pooling.
#' @param cPar numeric matrix; Initial chromosome-level parameter estimates.
#'  Rows correspond to model parameters (`beta.LD`, `beta.ovlp`) and columns
#'  to chromosomes.
#' @param gPar numeric vector; Global starting values for the fixed effects.
#'  Must contain genome-wide means for `beta.LD`, `beta.ovlp`, and the fixed
#'  LD decay scale.
#' @param cN0 data.table; Empirical covariance summaries by chromosome and lag.
#'  Must contain columns `chr`, `pos.lag`, `cov.hubM.mu`, `cov.hubM.se`,
#'  and `cov.ovlp`. Rows with `chr == "All"` are excluded internally.
#' @param corr logical; Default is `TRUE`, whereby a full 2×2 random-effects
#'  covariance matrix allowing correlation between LD and overlap amplitudes. If
#'  `FALSE` (default), random effects are assumed independent.
#' @param maxit integer; Maximum number of REML iterations.
#' @return A list with components:
#' \describe{
#'   \item{mu}{Estimated genome-wide mean amplitudes.}
#'   \item{beta.var}{Estimated between-chromosome covariance matrix of amplitudes.}
#'   \item{beta.chr}{Empirical Bayes (BLUP) chromosome-specific estimates.}
#'   \item{iter}{Number of iterations to convergence.}
#' }
#' @details
#'  The model assumes a fixed exponential LD decay scale and a fixed overlap
#'  geometry. Observation-level uncertainty is accounted for via the empirical
#'  covariance standard errors (`cov.hubM.se`), which enter the REML estimation
#'  through a generalized least-squares formulation. Additional lag-bin weights
#'  from the covariance-fitting stage are not required in this step unless
#'  alternative weighting schemes are explicitly desired.
#' @import data.table
#' @noRd
.empirical_bayes_estimate <- function(cPar, gPar, cN0, corr = FALSE, maxit = 1e3) {
   # extract values
   cN0 <- data.table::copy(cN0)[chr != "All"]
   nChr <- ncol(cPar)
   cov_emp <- t(matrix(cN0$cov.hubM.mu, ncol = nChr)) # observed covariances
   X <- list(t(matrix(exp(-cN0$pos.lag / gPar[3]), ncol = nChr)),
             t(matrix(cN0$cov.ovlp, ncol = nChr))) # design matrices
   cN0[, se2 := cov.hubM.se^2]
   Vr <- t(matrix(cN0[, pmax(se2, quantile(se2, 0.05)), by = chr]$V1,
                      ncol = nChr)) ## observation variances
   Xv <- lapply(1:nChr, function(cI) {
      Xc <- cbind(X[[1]][cI, ], X[[2]][cI, ])
      Xv <- Xc / Vr[cI, ]
      list(Xvc = Xv, A = crossprod(Xv, Xc))
   })

   # initial values
   mu <- gPar[1:2]
   tau2 <- rep(0.01, 2)
   rho <- 0

   for (it in 1:maxit) {
      tau2.old <- tau2
      rho.old  <- rho
      res <- cov_emp - mu[1] * X[[1]] - mu[2] * X[[2]] # residual after fixed effects

      # build G
      G <- matrix(c(tau2[1], rep(rho * sqrt(prod(tau2)), 2), tau2[2]), nrow = 2)
      Ginv <- solve(G + diag(1e-9, 2)) # ensure invertibility

      # joint BLUP update
      uHat <- sapply(1:nChr, function(j) {
         b <- crossprod(Xv[[j]]$Xvc, res[j, ])
         solve(Ginv + Xv[[j]]$A, b)
      })

      ## update fixed effects using raw BLUPs
      mu.c <- rowMeans(uHat)
      mu <- mu + mu.c
      uHat <- uHat - mu.c # ensure estimates are centred

      # REML variance update
      tau2 <- rowSums(uHat ^ 2) / (nChr - 1)
      if (corr) { # update correlation if requested
         rho <- cov(uHat[1, ], uHat[2, ]) / sqrt(prod(tau2))
         rho <- max(-0.999, min(0.999, rho))
      }

      # check convergence
      tauCv <- max(abs(tau2 - tau2.old)) < 1e-6
      rhoCv <- ifelse(corr, abs(rho - rho.old) < 1e-6, TRUE)
      if (tauCv & rhoCv) break
   }

   # calculate REML for optimal parameters
   res <- cov_emp - mu[1] * X[[1]] - mu[2] * X[[2]]
   G <- matrix(c(tau2[1], rep(rho * sqrt(prod(tau2)), 2), tau2[2]), nrow = 2)
   Ginv <- solve(G + diag(1e-9, 2)) # ensure invertibility
   dG <- determinant(G, logarithm = TRUE)$modulus
   logLik <- sum(sapply(1:nChr, function(cI){
      aI  <- Xv[[cI]]$A
      uhI <- uHat[, cI]
      vrI <- Vr[cI, ]
      dGaI <- determinant(Ginv + aI, logarithm = TRUE)$modulus
      qf <- sum(res[cI, ] ^ 2 / vrI) - crossprod(uhI, crossprod(Ginv + aI, uhI))
      -(sum(log(vrI)) + qf + dG + dGaI) / 2
   }))

   # return results
   betaChr <- lapply(1:2, function(j) (mu[j] + uHat[j, ]))
   names(mu) <- names(betaChr) <- c("beta.LD", "beta.ovlp")
   betaVar <- matrix(c(tau2[1], rep(rho * sqrt(prod(tau2)), 2), tau2[2]),
                     nrow = 2)
   rownames(betaVar) <- colnames(betaVar) <- c("beta.LD", "beta.ovlp")

   return(list(mu = mu, beta.var = betaVar, beta.chr = betaChr,
               logLik = logLik, iter = it))
}


#' @title Calculate autocovariance for null gene sets
#' @description Internal function to efficiently compute the autocovariance for
#'  all gene sets in a null model. Uses data.table to efficiently extract and
#'  sum covariances for all gene pairs in in each null gene set (which are used
#'  for rescaling null scores).
#' @param csj integer; A vector of object IDs.
#' @param rsX matrix; The covariance matrix.
#' @param SS0 data.table; A table of gene sets and their gene IDs.
#' @param ac0 data.table; A table of gene-gene autocovariances.
#' @import data.table
#' @return A data.table of autocovariances (`C`) for each null gene set (`sID`)
#'   and by successive set sizes (`sN`)
#' @noRd
.get_covariance_null <- function(csj, rsX, set_scores, autocov_data) {
   SSi <- set_scores[, oID := c(csj)]
   ACx <- autocov_data[SSi, on = .(A = oID), nomatch = 0, allow.cartesian = TRUE]
   ACxx <- ACx[SSi, on = .(B = oID, sID = sID), nomatch = 0]
   ACxx[, sN0 := ifelse(sN > i.sN, sN, i.sN)] # get last appearance of gene and adjust for current minimum
   RSi <- ACxx[, .(CV = sum(CV)), keyby = .(sID, sN0)] # calculate sum of covariance and sort
   rsX[cbind(RSi$sN0, RSi$sID)] <- RSi$CV # reassign covariances to matrix
   return(rsX)
}


#' @title Calculate autocovariance for gene sets
#' @description Internal function to efficiently compute the autocovariance for
#'  gene sets. Uses data.table to efficiently extract and sum covariances for
#'  all gene pairs in in each gene set (which are used for rescaling scores).
#' @param SSi data.table; A table of gene sets and their gene IDs.
#' @param ac0 data.table; A table of gene-gene autocovariances.
#' @import data.table
#' @return A data.table of autocovariances (`C`) for each null gene set (`sID`)
#' @noRd
.get_covariance <- function(SSi, autocov_data) {
   # extract covariances
   ACx <- autocov_data[SSi, on = .(A = oID), nomatch = 0, allow.cartesian = TRUE]
   ACxx <- ACx[SSi, on = .(B = oID, sID = sID), nomatch = 0]
   # calculate sum of covariance and sort
   RSi <- ACxx[, .(CV = sum(CV)), keyby = sID]
   return(RSi)
}


#' @title Correct gene set scores for shared genes
#' @description Internal function performs a sequential gene pruning algorithm
#'  to correct for shared genes among gene sets. It iteratively identifies the
#'  most significant gene set, removes its genes from all other sets, and
#'  re-calculates the scores and p-values for the affected sets.
#' @param OBS numeric; A vector of observed set scores.
#' @param EXP list; A list of empirical CDF functions for each gene set size.
#' @param GPD data.table; A table of GPD parameters for each gene set size.
#' @param acX data.table; The autocovariance table.
#' @param pX numeric; A vector of p-values for each set.
#' @param nX integer; A vector of gene counts for each set.
#' @param osX numeric; A vector of observed gene scores.
#' @param PI data.table; A table of gene to set mappings.
#' @param fpc numeric; A vector of scaling factors for each gene set size.
#' @param mss numeric; The minimum gene set size.
#' @import data.table
#' @return A `matrix` summarizing the pruning results, including ranks, set
#'  sizes, scores, and p-values for the top sets.
#' @noRd
.prune_gene_sets <- function(OBS, EXP, GPD, acX, pX, nX, osX, PI, finite_pop_correction, min_set_size, gP) {
   # keep track top gene sets and their genes
   n.set.rem <- max(PI$A) # no. gene sets remaining
   n.obj.rem <- length(osX) # no. genes remaining
   si.out <- matrix(0, nrow = n.set.rem, ncol = 6)
   rescaled <- !is.null(acX)
   if (rescaled) { # track autocovariance
      cX <- rep(0, n.set.rem)
      cx0 <- .get_covariance(SSi = PI[A == B, .(sID = A, oID = X)], autocov_data = acX)
      cX[cx0$sID] <- cx0$CV
   }

   K <- 0
   while (n.set.rem > 1) { # continue while more than one valid set remains
      K <- K + 1 # increment

      # determine most significant gene set
      ts <- which(pX == min(pX)) # don't use which.min to allow for ties
      if (length(ts) > 1) { # use largest gene set in case of ties
         ts <- ts[which.max(nX[ts])]
      }

      # record top gene set values
      n.ts <- nX[ts] # set size
      p.ts <- pX[ts] # p value
      if (rescaled) {
         ss.ts <- OBS[ts] / sqrt((n.ts + cX[ts]) * finite_pop_correction[n.ts]) # set score
      } else {
         ss.ts <- OBS[ts] / sqrt(n.ts * finite_pop_correction[n.ts]) # set score
      }

      # determine genes and sets to update
      nXY0 <- PI[.(ts)] # shared genes with top scoring gene set
      obj.out <- nXY0[A == B]$X # genes to remove
      ss.u <- unique(nXY0$B) # gene sets sharing genes with top scoring gene set
      n.new <- nX[ss.u] - nXY0[, .N, by = B]$N # decrement to gene set size
      nX[ss.u] <- n.new # update gene set size

      # remove sets with insufficient genes after removing shared genes
      ws.out <- n.new < min_set_size
      set.out <- ss.u[ws.out]
      pX[set.out] <- 1 # invalidate removed sets from minimum p values in subsequent iterations
      nXY <- nXY0[!.(set.out), on = "B"] # determine gene sets with sufficient genes to modify

      # revise scores
      if (nrow(nXY) > 0) {
         ss.i <- unique(nXY$B) # gene sets to update
         ss.j <- unique(nXY$X) # genes to remove

         # calculate decrement to summed scores
         if (length(ss.i) == 1) { # only one gene set with shared genes
            obs.dec <- sum(osX[ss.j]) # decrement to sum
         } else {
            obs.dec <- sapply(split(osX[nXY$X], nXY$B), sum)
         }

         obs.new <- OBS[ss.i] - obs.dec # revised sum
         OBS[ss.i] <- obs.new # update sums
         nU <- n.new[!ws.out]

         if (rescaled) { # add autocovariance to variance
            ss.k <- intersect(ss.i, which(cX > 0))
            ssi0 <- PI[.(ss.k, ss.k), .(sID = A, oID = X)] # need to consider all genes in each updated gene set
            ACx <- acX[ssi0, on = .(A = oID), nomatch = 0, allow.cartesian = TRUE]
            ACxx <- ACx[ssi0, on = .(B = oID, sID = sID), nomatch = 0]
            ACxx <- ACxx[A %in% ss.j | B %in% ss.j] # only keep interactions involving removed genes
            cx0 <- ACxx[, .(CV = sum(CV)), keyby = sID]
            cX[cx0$sID] <- cX[cx0$sID] - cx0$CV
            sso <- obs.new / sqrt((nU + cX[ss.i]) * finite_pop_correction[nU]) # revised set score
         } else {
            sso <- obs.new / sqrt(nU * finite_pop_correction[nU]) # revised set score
         }

         # estimate p values
         for(sUi in split(data.frame(ss.i, sso, nU), nU)) {
            pX[sUi$ss.i] <- with(sUi, gP(ss = sso, n = nU[1], g0 = GPD, e0 = EXP))
         }
      }

      # update gene and gene set sizes
      PI <- PI[!.(set.out)][!(B %in% set.out | X %in% obj.out)] # remove unused sets and genes
      n.set.rem <- n.set.rem - length(set.out)
      n.obj.rem <- n.obj.rem - length(obj.out)
      si.out[ts, ] <- c(K, n.obj.rem, n.set.rem, n.ts, ss.ts, p.ts) # record results
   }

   si.out
}


#' @title Compute FDR bin counts
#' @description This internal function calculates the number of observed and
#'  expected p-values within a series of predefined bins. These counts are
#'  used to estimate the false discovery rate (FDR) and the proportion of true
#'  null hypotheses.
#' @param O data.table; The table of observed p-values.
#' @param E data.table; The table of p-values from the null distribution.
#' @param p.bins numeric; A vector of boundaries for the p-value bins.
#' @param n.fdr numeric; The number of null replications used to generate the
#'  expected p-values.
#' @import data.table
#' @return A data.table; containing the observed and expected counts (`OBS` and
#'  `EXP`) for each p-value bin.
#' @noRd
.fdr_bin_counts <- function(O, E, p.bins, n.fdr) {
   O[, FDR.bin := cut(setScore.pr.p, p.bins)]
   rp <- O[, .(OBS = .N), by = FDR.bin]

   n.rs <- nrow(O) / n.fdr # expected count rescaling factor
   E[, FDR.bin := cut(setScore.pr.p, p.bins)]
   e.bin <- E[, .N, by = c("FDR.rep", "FDR.bin")]
   e.bin[, N.tot := sum(N), by = FDR.rep]
   vp <- e.bin[!is.na(FDR.bin), .(EXP = sum(N / N.tot) * n.rs), by = FDR.bin]

   #merge and fill in empty bins
   out <- merge(rp, vp, by = "FDR.bin", all = TRUE)
   out[is.na(OBS), OBS := 0]
   out[is.na(EXP), EXP := 0]

   return(out)
}


#' @title Estimate pi0 for FDR correction
#' @description This internal function estimates `pi0`, the proportion of true
#'  null hypotheses, for false discovery rate (FDR) correction. It uses an
#'  iterative algorithm to find a stable estimate based on the comparison of
#'  observed and expected p-value counts in specific bins.
#' @param n.fdr.sets numeric; The total number of gene sets used for FDR
#'  estimation.
#' @param SSO.pr data.table; The table of observed pruned set scores.
#' @param FDR.pr data.table; The table of pruned set scores from the null
#'  distribution.
#' @param tolerance numeric; The convergence tolerance for the pi0 estimation.
#' @param n.fdr numeric; The number of null replications used.
#' @param env environment; The environment where the estimated `pi0` and the bin
#'  counts table will be stored.
#' @import data.table
#' @return No value returned. Instead, `pi0` and `pi0.dt` objects are returned
#'  to the environment (`env`).
#' @noRd
.estimate_pi0 <- function(n.fdr.sets, SSO.pr, FDR.pr, tolerance, n.fdr, env) {
   n.cuts <- min(floor(n.fdr.sets), 100) + 1 # number of bins up to max 100
   pi0.p.bins <- seq(0, 1, length.out = n.cuts)

   # create percentile bins and calculate observed and expected values
   pi0.dt <- .fdr_bin_counts(O = SSO.pr[!is.na(Rank)],
                             E = FDR.pr[!is.na(Rank)],
                             p.bins = pi0.p.bins, n.fdr = n.fdr)

   # determine obs < exp up to first failure and evaluate probabilities in prior bins
   obs.excess <- pi0.dt[, which(cumsum(OBS < EXP) == 1) - 1][1]
   pi0 <- 1 # initial value
   if(obs.excess > 0) {
      track.pi0 <- pi0
      k <- 0; TT <- TRUE
      pi0.dt[, EXP.orig := EXP]
      while (isTRUE(TT)) { # run to convergence
         k <- k + 1
         true.pos <- pi0.dt[1:obs.excess, sum(OBS - EXP) / n.fdr.sets]
         pi0 <- 1 - true.pos
         pi0.dt[, EXP.scaled := pi0 * EXP.orig]
         obs.excess <- pi0.dt[, which(cumsum(OBS < EXP) == 1) - 1][1]
         track.pi0 <- c(track.pi0, pi0)
         TT <- abs(track.pi0[k] - pi0) > tolerance # test convergence
      }
   }

   list2env(list(pi0 = pi0, pi0.dt = pi0.dt), envir = env)
}


#' @title Determine optimal ggplot facet dimensions
#' @description Internal function calculates the optimal number of rows and
#'  columns for a `ggplot2` facet plot. Considers the total number of panels, a
#'  maximum number of facets per page, and a set of heuristic criteria to
#'  minimize unused space, keep the layout balanced, and limit the number of
#'  pages.
#' @param n.panel numeric; The total number of panels (e.g., gene sets) to be
#'  displayed across all pages.
#' @param max.facets numeric; The maximum number of facets to display on a
#'  single page.
#' @import data.table
#' @return A list containing three numeric values: `n.row` (the optimal number
#'  of rows per page), `n.col` (the optimal number of columns per page), and
#' `n.pg` (the total number of pages required).
#' @noRd
.facet_opt <- function(n.panel, max.facets) {
   max.cols <- min(7, n.panel)
   OPT <- min(n.panel, max.facets)
   ffX <- outer(seq_len(ceiling(OPT / 2)), seq_len(max.cols))
   pg.dim <- which(ffX <= max.facets, arr.ind = TRUE)
   pg.dim <- as.data.frame(pg.dim)
   data.table::setDT(pg.dim)
   data.table::setnames(pg.dim, c("n.row", "n.col"))
   pg.dim <- pg.dim[n.col >= n.row]
   pg.dim[, n.facet := n.col * n.row]
   pg.dim[, n.page := ceiling(n.panel / n.facet)]

   if (nrow(pg.dim) > 1) { #  more than one option
      pg.dim[, unused := n.page * n.facet - n.panel] # no. of unused facets
      pg.dim <- pg.dim[unused < n.col] # remove options with missing final row
      pg.dim[, f.dist := max.facets - n.facet] # distance between max.facets and facets per page
      pg.dim[, pg.wt := log(n.page / (n.row * n.col))] # measure of page no. to facet weight
      pg.dim[, p1 := unused / ifelse(var(unused) == 0, 1, stats::sd(unused))]
      pg.dim[, p2 := f.dist / ifelse(var(f.dist) == 0, 1, stats::sd(f.dist))]
      pg.dim[, p3 := pg.wt / ifelse(var(pg.wt) == 0, 1, stats::sd(pg.wt))]
      pg.dim[, WT := p1 + p2 + p3]
      pg.dim <- pg.dim[WT == min(WT)]

      if (nrow(pg.dim) > 1) { # still more than one option
         pg.dim <- pg.dim[abs(n.row - n.col) == min(abs(n.row - n.col))] # similar rows and columns
      }
   }
   return(list(n.row = pg.dim$n.row, n.col = pg.dim$n.col,
               n.pg = pg.dim$n.page))
}


#' @title Convert function arguments to a plot title
#' @description Internal function to create a title string from a list of
#'  function arguments for plotting.
#' @param args list; The list of arguments.
#' @import foreach
#' @return A character string formatted as a title.
#' @keywords internal
#' @noRd
.args2title <- function(args) {
   args <- args[grep("group|seed|input", names(args), invert = T)] # remove seed info
   args <- args[!sapply(args, is.null)]

   s.args <- split(args, ceiling(seq_along(args) / 4)) # split into chunks of 4
   out <- foreach::foreach(ARG = s.args) %do% {
      arg0 <- sapply(ARG, paste, collapse = ",")
      paste(names(ARG), "=", arg0, collapse = "; ")
   }
   return(paste(c("Parameters:", out), collapse = "\n"))
}

# Deprecated function wrappers for backward compatibility

#' @rdname verbose_msg
#' @keywords internal
#' @noRd
.vrb <- function(...) {
   .verbose_msg(...)
}

#' @rdname combine_results
#' @keywords internal
#' @noRd
.cmb <- function(...) {
   .combine_results(...)
}

#' @rdname fit_local_quad_reg
#' @keywords internal
#' @noRd
.fit_lqr <- function(...) {
   .fit_local_quad_reg(...)
}

#' @rdname fit_gpd_tail
#' @keywords internal
#' @noRd
.fit_gpd <- function(...) {
   .fit_gpd_tail(...)
}

#' @rdname get_quantiles
#' @keywords internal
#' @noRd
.get_quant <- function(...) {
   .get_quantiles(...)
}

#' @rdname get_gpd_minimums
#' @keywords internal
#' @noRd
.get_gpd_mins <- function(...) {
   .get_gpd_minimums(...)
}

#' @rdname get_pvalue_lower
#' @keywords internal
#' @noRd
.get_p_lt <- function(...) {
   .get_pvalue_lower(...)
}

#' @rdname get_pvalue_upper
#' @keywords internal
#' @noRd
.get_p_ut <- function(...) {
   .get_pvalue_upper(...)
}

#' @rdname get_covariance_null
#' @keywords internal
#' @noRd
.get_cov0 <- function(...) {
   .get_covariance_null(...)
}

#' @rdname get_covariance
#' @keywords internal
#' @noRd
.get_cov <- function(...) {
   .get_covariance(...)
}

#' @rdname prune_gene_sets
#' @keywords internal
#' @noRd
.prune_sets <- function(...) {
   .prune_gene_sets(...)
}

#' @rdname estimate_pi0
#' @keywords internal
#' @noRd
.est_pi0 <- function(...) {
   .estimate_pi0(...)
}

#' @rdname empirical_bayes_estimate
#' @keywords internal
#' @noRd
.emp_bayes_est <- function(...) {
   .empirical_bayes_estimate(...)
}

#' @rdname covariance_function
#' @keywords internal
#' @noRd
.cov_fun <- function(...) {
   .covariance_function(...)
}

#' @rdname fit_cgm_function
#' @keywords internal
#' @noRd
.cgm_fun_fit <- function(...) {
   .fit_cgm_function(...)
}

#' @rdname parallel_smooth
#' @keywords internal
#' @noRd
.par_smooth <- function(...) {
   .parallel_smooth(...)
}

#' @rdname check_file_paths
#' @keywords internal
#' @noRd
.path_check <- function(...) {
   .check_file_paths(...)
}

#' @rdname unpack_plr_object
#' @keywords internal
#' @noRd
.plR_unpack <- function(...) {
   .unpack_plr_object(...)
}

#' @rdname check_arguments
#' @keywords internal
#' @noRd
.arg_check <- function(...) {
   .check_arguments(...)
}
