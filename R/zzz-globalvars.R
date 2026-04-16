## Suppress R CMD check NOTEs for data.table NSE column bindings
##
## These are column names used in data.table expressions throughout the package.
## Declaring them here prevents "no visible binding" warnings and NOTES.
##
## NOTE: This file must be last in the collate order for globalVariables to work,
## hence the "zzz" prefix.
utils::globalVariables(c(

  # data.table column names used across plR_read / permute / rescale / prune
  ".", "A", "ARG", "B", "C", "CV", "EP.bp", "EP.cM", "EXP", "EXP.orig",

  "EXP.scaled", "F.i", "FDR.bin", "FDR.rep", "FUNCTION", "ID", "J", "MAD",
  "info.messages", "N", "N.ovlp", "N.tot", "Nm", "OBS", "R", "RP", "Rank", "SE",
  "SP.bp", "SP.cM", "TEMP", "V1", "VP", "warning.messages", "WT", "X",

  # plR object and attribute variables
  "ac", "adj.p", "alt", "bin", "bin.max", "bin.min", "bin.size", "bootN",
  "bp2cM", "cM", "cgm.bin", "cgm.range", "cgm.wt.max", "chr", "chr.orig",
  "coord", "cov.emp", "cov.hubM.se", "cov.info", "cov.labels", "cov.mat",
  "cov.mean", "cov.names", "cov.ovlp", "csN", "cut.z",

  # more variables from check output
  "eI", "eJ", "eR", "eS", "ecdf.rs", "ecdf.std", "ei", "ej", "emp.bayes",
  "endpos", "endpos.base", "est.pi0", "exc.z",

  "f.dist", "f.i", "file.paths",

  "gI", "gU", "get.objStat", "gpd.cutoff", "gpd.rs", "gpd.std", "gpdI",

  "i", "i.A", "i.B", "i.J", "i.sN", "i.setID", "j", "k",

  "kern.bound", "kern.func", "kern.scale", "kern.wt.max", "lag.wts",

  "m", "map.fun", "max.facets", "max.set.n", "md.meth", "midpos", "midpos.bp",
  "min.p", "min.rho", "min.set.n", "min.z", "mu",

  "n", "n.boot", "n.chr", "n.col", "n.cov", "n.facet", "n.fdr", "n.genes",
  "n.max", "n.page", "n.par", "n.perm", "n.row", "n.set.genes", "n.sets",
  "no.deconf", "no.share",

  "oID", "obj.buffer", "obj.in", "obj.info", "obj.out", "obj.stat.fun",
  "objBin", "objID", "objID.A", "objID.B", "objID.new", "objID.orig",
  "objMax", "objMax.res", "objMed", "objMed.res", "objStat", "objStat.res",
  "objStat.std", "osI", "osJ", "output.path",

  "p.i", "p0", "p1", "p2", "p3", "param.messages", "param.warnings", "par", "perm.path",
  "permute", "pg.wt", "plot.all", "plot.name", "plr.seed", "pos", "pos.info",
  "pos.lag", "prog.hand",

  "rate", "rec.rate", "rem.genes", "rescale", "rho",

  "sI", "sID", "sJ", "sN", "sN0", "se2", "set.genes", "set.in", "set.info",
  "set.info.merged", "set.merge", "set.obj", "set.out", "setG", "setID",
  "setID.new", "setID.new.y", "setID.orig", "setN", "setName",
  "setScore.pr.p", "setScore.rs", "setScore.rs.p", "setScore.std",
  "setScore.std.p", "startpos", "startpos.base",

  "tail", "tolerance",

  "unused", "user.ac",

  "value", "var.info", "vb.i",

  "x", "x0", "xN"
))