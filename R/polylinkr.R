#===========================================================================
# Function: get.null(x)
# permutation algorithm for PolyLinkR used to generate null distribution
#
# - set.info: data.frame file with fields:
#                  setID, setName, ...
# - set.obj : data.frame with fields:
#                  setID, objID
# - obj.info: data.frame with fields:
#                  objID, objName, objStat, (objBin), (objSNPcnt), ...
# - n.cores: integer specifying the number of cores to use in computation
# - emp.nruns: number of iterations of permatution algorithm to comput null
# - NN: subset of iterations on which to summarize results (improves computation speed); default = 1000
#
#===========================================================================

#permutation function
#' Compute null distribution and probability estimation of gene score enrichment in pathways or gene sets.
#'
#' This is the core function of PolyLinkR, it runs the permuation algorithm that
#' creates a random mapping of genomic scores while preserving the inherent linkage
#' disequilibrium amongt the different genomic regions.
#' This process performed iteratively to generate a null distribution for testing
#' enrichment in biological pathways or gene sets.
#'
#'
#' @param set.info data.frame: four required fields; gene IDs, associated scores, chromosome, and start position
#' @param obj.info data.frame: two required fields; genomic regions and associated scores
#' @param set.obj data.frame: two required fields; genomic regions and associated scores
#' @param n.cores integer: number of cores to run in parallel
#' @param emp.nruns integer: number of iterations used to compute null
#' @param NN integer: subset of interations for summary statistic calculations; default = 1000
#'
#' @return Returns a data.table.
#'
#' @seealso \code{\link{PolyLinkR_SetInfo}}, \code{\link{PolyLinkR_SetObj}}
#'
#'
#' @examples
#' output = polylinkr(obj.info = Anatolia_EF_CLR, set.info = PolyLinkR_SetInfo, set.obj = PolyLinkR_SetObj,
#'              n.cores = 8, emp.nruns = 1000, NN = 10000)
#'
#' @export
polylinkr <- function(set.info, obj.info, set.obj,
                         n.cores=NA, emp.nruns=10000, NN=1000){

  #check for common errors
  if(emp.nruns<=0){ 
    emp.nruns <- max(NN, 10000)
    warning(paste("WARNING: iterations incorrectly entered, resetting to", emp.nruns),
            immediate.=T)
  }
  if(NN>emp.nruns){
    NN <- min(emp.nruns, 1000)
    warning("WARNING: summary statistic evalution block size exceeds total number of iterations",
            immediate.=T)
    warning("summary statistic evalution set to arg.min(1000, emp.nruns)",
            immediate.=T)
  }
  
  #set cores
  dc <- detectCores()
  if(n.cores<=0 | n.cores>dc | is.na(n.cores)) n.cores <- dc-1
  
   print(paste("Running enrichment test..."))

   #permutation function
   permute.data <- function(obj.info, n.chr, n.genes, gene.pos, chr.ord.now){
      # 1. string together chromosomes
      chr.ord.now <- sample(1:n.chr, n.chr)
      new.ord <- unlist(gene.pos[chr.ord.now])
      r1 <- obj.info[new.ord, objStat]

      # 2. rotate scores
      rotate.now <- sample(2:n.genes, 1)
      r1[c(rotate.now:n.genes, 1:(rotate.now-1))]
   }

   #sumstat calculation function
   sum.stat <- function(set.obj, ID, perm){
      # 3. compute new sumstat scores
      mm.e <- merge(data.table(objID=ID, objStat=perm), set.obj, by="objID")[order(setID)]
      mm.e[, lapply(.SD, sum), .SDcols=grep("objStat", names(mm.e)), by="setID"]
   }

   #p.value calculation function
   compute.p.val <- function(obs, exp){
      # 4. compute p values
      e <- exp[, setID:=NULL]
      rowSums(obs <= e)
   }

   #timer function
   timer <- function(time.point, emp.nruns, I){
      if((I-1000) %% time.point == 0 & I>1000){ # estimate time remaining
         d.now <- as.numeric(difftime(Sys.time(), s0, units="secs"))
         est.time <- (emp.nruns-I+1000)*(d.now/(I-1000))
         print(paste("Completed", I-1000, "iterations"))
         mins <- est.time %/% 60
         secs <- round(est.time %% 60, 0)
         if(mins>=60){
            hours <- mins %/% 60
            mins <- mins %% 60
            print(paste("Estimated time remaining:", hours, "h", mins, "m", secs, "s"))
         }else{
            if(mins==0){
               print(paste("Estimated time remaining:", secs, "s"))
            }else{
               print(paste("Estimated time remaining:", mins, "m", secs, "s"))
            }
         }
      }
   }

   #break up emp.nruns into specific size iteration chunks for computation
   #(improves run time; 1000 seems optimal for ~10k genes)
   get.blocks <- function(emp.nruns, block.size){
      if(emp.nruns<block.size){
         rb <- emp.nruns
      }else{
         onek.blocks <- rep(block.size, emp.nruns %/% block.size)
         remainder <- emp.nruns %% block.size
         if(remainder>0){
            rb <- c(onek.blocks, remainder)
         }else{
            rb <- onek.blocks
         }
      }
      return(rb)
   }

   #convert to datatable format (if not already data.table)
   if(!is.data.table(set.info)) data.table::setDT(set.info)
   if(!is.data.table(set.info)) data.table::setDT(obj.info)
   if(!is.data.table(set.info)) data.table::setDT(set.obj)

   #remove na's
   no.scores <- obj.info[is.na(objStat), objID]
   obj.info <- obj.info[!(objID %in% no.scores)]
   set.obj <- set.obj[!(objID %in% no.scores)]

   #set up parallel backend for foreach %dopar%
   doParallel::registerDoParallel(cores=n.cores)

   # order obj.info table (needed to compare rotated values against)
   obj.info <- obj.info[order(chr, startpos)]
   ID <- obj.info$objID

   #get relevant variables
   n.genes <- obj.info[, .N]
   n.paths <- set.info[, .N]
   n.chr <- obj.info[, length(unique(chr))]

   #determine matrix position of each gene for each chromosome
   gene.pos <- foreach::foreach(i=1:n.chr) %do% obj.info[, which(chr==i)]

   # compute observed sumstat scores
   mm.o <- merge(set.obj, obj.info[, .(objID, objStat)], by="objID")
   mm.n <- mm.o[order(setID), .(N=length(unique(objID))),
                by=c("setID")]
   m.obs <- mm.o[order(setID), .(SumStat=sum(objStat, na.rm=T)),
                 by=c("setID")]

   #housekeeping
   rm(mm.o)
   gc()

   # compute expected sumstat scores
   run.blocks <- get.blocks(emp.nruns, NN)
   I=0
   sig.tests <- rep(0, n.paths)
   s0 <- Sys.time()
   for(l in run.blocks){
      I=I+l
      timer(NN, emp.nruns, I)
      pp <- foreach(i=1:l) %dopar% {
         permute.data(obj.info, n.chr, n.genes, gene.pos)
      }
      perm <- matrix(unlist(pp), ncol=l, byrow=FALSE)
      m.exp <- sum.stat(set.obj, ID, perm)
      sig.tests <- sig.tests + rowSums(m.obs[, SumStat] <= m.exp[, -1])
      #compute.p.val(obs=m.obs$SumStat, exp=m.exp)
      rm(perm, m.exp)
      gc()
   }

   ##-----------------------------------------##
   ## Compute p and q values
   ##-----------------------------------------##

   p.vals <- sig.tests/emp.nruns
   q.vals <- qvalue::qvalue(p.vals, pi0.method="smoother")

   data.table::data.table(set.info[order(setID), .(setID, setName)],
              setScore=m.obs$SumStat, setSize=mm.n$N,
              setP=p.vals, setQ=q.vals$qvalues)[order(setP)]

}

#read input data function
#' Read in polylink data.
#'
#' Reads in data. Optionally, can merge gene sets that share more than a specified proportion of genes, or remove gene sets with less/more than a specified number of genes.
#'
#'
#' @param in.path character: pathway to input files
#' @param population character: label used to identify input files. Three input files are required, specified by 'setInfo', 'setObj', and 'objInfo'.
#' @param minsetsize integer: minimum number of genes required in gene set
#' @param maxsetsize integer: maximum number of genes required in gene set
#' @param min.sim double: sets sharing at least this proportion of genes will be merged
#' @param n.cores integer: number of cores to run in parallel
#'
#' @return Returns a list with three elements.
#'
#' @export
#===========================================================================
# Function: ReadSetObjTables(in.path, set.info.file,set.obj.file,
#                            obj.info.file)
# Read in all required gene (object) and gene set (set) tables
#
# - in.path      : path to directory with input files. 
#                  default = local folder './'
# - population   : label specifying which dataset to load. Must be specified. 
#                  Will search for the following files:
#                 - setInfo: tab seperated file with fields:
#                  setID, setName, ...
#                 - setObj : tab seperated file with fields:
#                  setID, objID
#                 - objInfo: tab seperated file with fields:
#                  objID, objName, objStat, chr, start, ...
# - minsetsize   : exclude gene sets with size below minsetsize
#                  (default = 10)
# - maxsetsize   : exclude gene sets with size above maxsetsize
#                  (default = 1000)
# - min.sim      : minimum proportion of shared genes used to 
#                  concatenate gene sets (default = 1, i.e. no concatenation)
# - n.cores      : no. of computation cores to use (default = no. cores - 1)
#
# These files must contain headers, IDs can be strings
# Internal numeric IDs will be assigned to objects and sets to improve
# further computations
#===========================================================================
ReadSetObjTables<-function(in.path="./", population=NA, 
                           minsetsize=10, maxsetsize=1000,
                           merge.set.prop=1, n.cores=NA){
  
  #Reading in data
  print(paste("Reading data for", population))
  
  #determine input files
  ll <- list.files(in.path)
  pf <- ll[grep(population, ll)]
  
  if(merge.set.prop>1 | merge.set.prop<=0){
    merge.set.prop <- min(emp.nruns, 1)
    warning("WARNING: minimim (0) or maximum (1) gene set similarity exceeded",
            immediate.=T)
    warning("merge.set.prop set to 1",
            immediate.=T)
  }

  # Read in information on gene sets
  set.info <- fread(file.path(in.path, pf[grep("SetInfo", pf)]))
  set.obj <- fread(file.path(in.path, pf[grep("SetObj", pf)]))
  obj.info <- fread(file.path(in.path, pf[grep("ObjInfo", pf)]))

  # Cleaning data
  #merge similar sets
  if(merge.set.prop<1){
    print(paste0("Merging gene sets with >", merge.set.prop, " similarity"))
    r <- MergeSimilarSets(SI=set.info, SO=set.obj,
                          min.sim=merge.set.prop, 
                          n.cores=n.cores)
    set.info <- r$set.info
    set.obj <- r$set.obj
  }
  
  # Remove resized gene sets that have too many/too few genes
  if(!is.na(minsetsize) | !is.na(maxsetsize)){
    setN <- set.obj[, .N, by=setID][N>minsetsize & N<maxsetsize, setID]
    set.info <- set.info[setID %in% setN]
    set.obj <- set.obj[setID %in% setN]
  }
  
  # Create new field with setName and setSource
  # to tell apart sets with the same name
  dups <- r
  if(length(dups)>0){
    print("Relabeling duplicated gene set names")
    dup.set <- set.info[, .N, by=setName][N>1, setName]
    for(d.now in dup.set){
      set.info[setName==d.now, setName:=paste0(setName, " set", 1:.N)]
    }
  }
  
  #add in setIDs for unmerged sets
  set.obj[is.na(setID.merged), setID.merged:=setID]
  
  return(list(set.info=set.info, set.obj=set.obj, obj.info=obj.info))
}


#===========================================================================
# Function: MergeSimilarSets(set.info, set.obj)
# Merge gene sets that have more than 95% similarity
#
# -SI : dataframe with fields setID, setName, ...
# -SO  : dataframe with fields setID, objID
# -min.sim  : minimum proportion of gene sharing before merging
# -n.cores  : number of cores to use
#===========================================================================


MergeSimilarSets<- function(SI, SO, min.sim=0.95, n.cores=NA){

  # Get similarity matrix
  # Which sets are min.sim proportion similar (two way)
  # Choose the one with largest original set to keep
  # Remove rest, but keep link in set.info.old
  
  #set cores
  dc <- detectCores()
  if(n.cores<=0 | n.cores>dc | is.na(n.cores)) n.cores <- dc-1
  
  #create set.obj matrix
  SO[, X:=1]
  ss <- dcast(SO, objID~setID, value.var="X",
              fun.aggregate=length)
  s.mat <- as.matrix(ss, nrow=nrow(ss), rownames="objID")
  s.mat.t <- t(s.mat)
  
  NP <- ncol(s.mat)
  #determine number of shared genes
  np <- 100
  if(NP>np){
    #speed up computation by using multiple cores and blocks
    #split matrix into blocks of length 100 (subset gene sets)
    #set up parallel backend for foreach %dopar%
    doParallel::registerDoParallel(cores=n.cores)
  
    sim.mat <- foreach(i=1:ceiling(NP/np), .combine=rbind) %dopar% {
      print(paste0("Quantified gene overlap in ", i*np, " of ", NP, " gene sets"))
      
      s.mat.t[(((i-1)*np)+1):min(i*np, NP), ] %*% s.mat
    }
  }else{
    sim.mat <- t(s.mat) %*% s.mat
    print(paste0("Quantified gene overlap in ", NP, " gene sets"))
  }
  
  set.n <- diag(sim.mat)
  
  #proportion shared
  p.mat <- sim.mat/set.n
  
  #check for similarity prop > min.sim in both pathways
  ff <- foreach(i=seq_len(NP), .combine=rbind) %dopar% {
    data.table(P=i, X=which(p.mat[i, ]>min.sim))
  }
  
  xx <- rbind(data.table(Z=1, ff), data.table(Z=2, P=ff$X, X=ff$P))
  xx[, PX:=paste(P, X)]
  keep <- xx[X!=P, .N, by="PX"][N>1, PX]
  
  set.m <- xx[Z==1 & PX %in% keep, .(P, X)]
  
  #concatenate pathways
  pc <- list()
  cc=0
  for(i in set.m[, unique(P)]){
    PC <- c(i, set.m[P==i, unique(X)])
    if(!any(PC %in% unlist(pc))){
      cc=cc+1
      pc[[cc]] <- PC
    }
  }
  
  if(length(pc)>0){
    so.new <- foreach(i=pc, .combine=rbind) %dopar% {
      so.now <- SO[setID %in% i]
      #rename according to largest set 
      #take random if equivalent
      new.setID <- so.now[, .N, by=setID][N==max(N), setID][1]
      data.table(setID=new.setID, objID=so.now[, unique(objID)], 
                 setID.merged=so.now[, paste0(unique(setID), collapse=",")])
    }
    
    SO.out <- rbind(SO[!(setID %in% unlist(pc)), .(setID, objID, setID.merged=NA)], 
                    so.new)[order(setID, objID)]
    SI.out <- SI[setID %in% SO.out$setID]
    SI.out <- SI.out[setID %in% unlist(pc), setName:=paste0(setName, "*")]
    
    print(paste0("Merged ", length(pc), " gene sets >", min.sim, " similarity"))

  }else{
    print(paste0("No gene sets >", min.sim, "similarity"))
  }

  return(list(set.obj=SO.out, set.info=SI.out))
}


