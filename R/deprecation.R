#' Deprecation helpers for polylinkR tidyverse naming transition
#'
#' @description
#' Internal helper functions for managing the transition from legacy naming
#' conventions to tidyverse-style naming in polylinkR v0.6.0. These functions
#' provide soft deprecation with informative warning messages.
#'
#' @keywords internal
#' @noRd

#' Deprecate a function name
#'
#' Creates a deprecated alias for a function that calls the new function
#' with an informative warning.
#'
#' @param old_name Character string of the deprecated function name
#' @param new_name Character string of the new function name
#' @param ... Arguments passed to the new function
#' @return Result from calling the new function
#'
#' @examples
#' \dontrun{
#' # Old function definition:
#' plR_read <- function(...) {
#'   .deprecate_function("plR_read", "read_polylinkr_data", ...)
#' }
#' }
#'
#' @keywords internal
.deprecate_function <- function(old_name, new_name, ...) {
  .Deprecated(new_name, package = "polylinkR",
              msg = paste0("'", old_name, "' is deprecated. ",
                          "Use '", new_name, "' instead. ",
                          "See ?", new_name, " for details."))
  
  # Get the new function from the package namespace
  new_fun <- get(new_name, envir = asNamespace("polylinkR"))
  
  # Call the new function with all provided arguments
  new_fun(...)
}


#' Deprecate a parameter name
#'
#' Checks if a deprecated parameter name was used and issues a warning,
#' returning the value mapped to the new parameter name.
#'
#' @param arg_list List of arguments from the calling function
#' @param old_param Character string of the deprecated parameter name
#' @param new_param Character string of the new parameter name
#' @return The value for the parameter (from old or new name), or NULL
#'
#' @examples
#' \dontrun{
#' # In function definition:
#' read_polylinkr_data <- function(input_path = NULL, input.path = NULL) {
#'   # Handle deprecated parameter
#'   if (!is.null(input.path) && is.null(input_path)) {
#'     input_path <- .deprecate_param(
#'       list(input.path = input.path, input_path = input_path),
#'       "input.path", "input_path"
#'     )
#'   }
#'   # ... rest of function
#' }
#' }
#'
#' @keywords internal
.deprecate_param <- function(arg_list, old_param, new_param) {
  if (!is.null(arg_list[[old_param]])) {
    if (!is.null(arg_list[[new_param]])) {
      # Both provided - use new one but warn
      warning("Both '", old_param, "' and '", new_param, 
              "' provided. Using '", new_param, "'.",
              call. = FALSE)
      return(arg_list[[new_param]])
    } else {
      # Only old param provided - deprecate and return
      .Deprecated(new_param, package = "polylinkR",
                  msg = paste0("Parameter '", old_param, "' is deprecated. ",
                              "Use '", new_param, "' instead."))
      return(arg_list[[old_param]])
    }
  }
  # Return new param value (may be NULL)
  return(arg_list[[new_param]])
}


#' Map deprecated parameter names to new names
#'
#' Takes a list of arguments and maps any deprecated parameter names
#' to their new equivalents, issuing deprecation warnings.
#'
#' @param arg_list List of arguments (usually from match.call or list(...))
#' @param param_map Named character vector where names are old parameters
#'   and values are new parameters
#' @return Modified argument list with deprecated parameters renamed
#'
#' @examples
#' \dontrun{
#' # Define parameter mapping
#' param_map <- c(
#'   "input.path" = "input_path",
#'   "obj.info.path" = "object_info_path",
#'   "n.perm" = "n_permutations"
#' )
#' 
#' # In function:
#' args <- .map_deprecated_params(args, param_map)
#' }
#'
#' @keywords internal
.map_deprecated_params <- function(arg_list, param_map) {
  for (old_name in names(param_map)) {
    new_name <- param_map[[old_name]]
    if (!is.null(arg_list[[old_name]])) {
      if (!is.null(arg_list[[new_name]])) {
        warning("Both '", old_name, "' and '", new_name, 
                "' provided. Using '", new_name, "'.",
                call. = FALSE)
      } else {
        .Deprecated(new_name, package = "polylinkR",
                    msg = paste0("Parameter '", old_name, "' is deprecated. ",
                                "Use '", new_name, "' instead."))
        arg_list[[new_name]] <- arg_list[[old_name]]
        arg_list[[old_name]] <- NULL
      }
    }
  }
  arg_list
}


#' Create a deprecated function alias
#'
#' Helper to create a deprecated wrapper around a new function.
#' This is used to maintain backward compatibility while warning users.
#'
#' @param old_name Character string - name of the deprecated function
#' @param new_name Character string - name of the new function
#' @param new_fun Function object - the new function to wrap
#' @return A function that calls the new function with deprecation warning
#'
#' @keywords internal
.make_deprecated_alias <- function(old_name, new_name, new_fun) {
  function(...) {
    .Deprecated(new_name, package = "polylinkR",
                msg = paste0("'", old_name, "' is deprecated and will be removed ",
                            "in a future version. Use '", new_name, "' instead."))
    new_fun(...)
  }
}


#' Deprecation message helper
#'
#' Creates a standardised deprecation message for functions or parameters.
#'
#' @param old_name Character string - deprecated name
#' @param new_name Character string - replacement name
#' @param what Character string - "function" or "parameter"
#' @return Character string with formatted deprecation message
#'
#' @keywords internal
.deprecation_msg <- function(old_name, new_name, what = "function") {
  paste0("The ", what, " '", old_name, "' is deprecated as of polylinkR v0.6.0. ",
         "Please use '", new_name, "' instead. ",
         "The old name will be removed in v1.0.0.")
}
