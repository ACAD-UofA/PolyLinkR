# Environment compatible with polylinkR:::.arg_check(f = "read", ENV = e).
read_arg_env <- function(...) {
  e <- new.env(parent = asNamespace("polylinkR"))
  defs <- list(
    input.path = NULL,
    obj.info.path = NULL,
    set.info.path = NULL,
    set.obj.path = NULL,
    var.info.path = NULL,
    rec.rate.path = NULL,
    min.set.n = 2L,
    max.set.n = Inf,
    group = NULL,
    map.fun = "kosambi",
    obj.buffer = "auto",
    obj.stat.fun = "non.param",
    bin.size = 250L,
    obj.in = NULL,
    obj.out = NULL,
    set.in = NULL,
    set.out = NULL,
    set.merge = 0.95,
    rem.genes = FALSE,
    verbose = TRUE
  )
  ovr <- list(...)
  for (nm in names(defs)) {
    e[[nm]] <- defs[[nm]]
  }
  for (nm in names(ovr)) {
    e[[nm]] <- ovr[[nm]]
  }
  e
}
