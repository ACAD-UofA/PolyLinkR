
#===========================================================================
# Function: read.setobj.tables(in.path, population, minsetsize, maxsetsize,
#                            merge.set.prop, n.cores)
# JD: changed function name following naming convention
# Read in all required gene (object) and gene set (set) tables
#
# JD: changed function to allow for input of set.info, obj.info and set.obj
# instead of reading obj's from files
# this allows for some proprocessing before calling this function, such
# as binning obj.stat and rescaling obj.stat
#
# - in.path        : path to directory with input files. 
#                    default = local folder './'
# - population     : label specifying which dataset to load. Must be specified. 
#                    Will search for the following files:
#                    - setInfo: tab seperated file with fields:
#                      setID, setName, ...
#                    - setObj : tab seperated file with fields:
#                      setID, objID
#                    - objInfo: tab seperated file with fields:
#                      objID, objName, objStat, chr, start, ...
# - minsetsize     : exclude gene sets with size below minsetsize
#                    (default = 10)
# - maxsetsize     : exclude gene sets with size above maxsetsize
#                    (default = 1000)
# - merge.set.prop : minimum proportion of shared genes used to concatenate 
#                    gene sets (default = 0.95, 1 = no concatenation)
# - n.cores        : no. of computation cores to use (default = no. cores - 1)
# JD: non.test.sets is not used anymore, right?
# - non.test.sets  : names of gene sets that will provide genes to the null 
#                    computations, but will not be directly tested for 
#                    enrichment
#
# These files must contain headers, IDs can be strings
# 
# Note: numeric IDs are internally assigned to objects and sets to improve 
# computation
#===========================================================================

read.setobj.tables <- function(in.path="./", set.info=NULL, 
                               obj.info=NULL, set.obj=NULL,
                               population=NA, minsetsize=10, 
                               maxsetsize=1000, merge.set.prop=0.95,
                               obj.in.set=T, n.cores='default'){
  

  if (is.null(set.info)|is.null(obj.info)|is.null(set.obj)){
    # Reading in data
    cat(paste0("Reading data for ", population, "\n"))
    # Determine input files
    ll <- list.files(in.path)
    PF <- ll[grep(population, ll)]
    
    if(length(PF)<3){
      miss.f <- setdiff(1:3, grep(".ObjInfo|.SetInfo|.SetObj", PF))
      ff <- paste(paste0(population, 
                         c(".ObjInfo", ".SetInfo", ".SetObj")[miss.f]), 
                  collapse=" & ")
      stop(paste0(ff, "were not detected. Check path and ensure correct ",
                  "file names are used."))
    }
    # Read in information on gene sets
    set.info <- data.table::fread(file.path(in.path, PF[grep("SetInfo", PF)]))
    set.obj <- data.table::fread(file.path(in.path, PF[grep("SetObj", PF)]))
    obj.info <- data.table::fread(file.path(in.path, PF[grep("ObjInfo", PF)]))
  } else {
    if(!is.data.table(set.info)) data.table::setDT(set.info)
    if(!is.data.table(obj.info)) data.table::setDT(obj.info)
    if(!is.data.table(set.obj)) data.table::setDT(set.obj)
  }
 
  if(merge.set.prop>1 | merge.set.prop<=0){
    merge.set.prop <- min(n.perm, 1)
    warning("WARNING: minimim or maximum (0,1] gene set similarity exceeded",
            immediate.=T)
    warning("merge.set.prop set to 1 (i.e. no merging performed)",
            immediate.=T)
  }

  # Cleaning data
  # Merge similar sets
  if(merge.set.prop<1){
    cat(paste0("Merging gene sets with >", merge.set.prop, " similarity\n"))
    r <- merge.similar.sets(SI=set.info, SO=set.obj,
                            min.sim=merge.set.prop, 
                            n.cores=n.cores, bfs=T)
    set.info <- r$set.info
    set.obj <- r$set.obj
  }else{
    set.obj[, setID.merged:=setID]
  }
  
  # Remove resized gene sets that have too many/too few genes
  cat(paste("Removing gene sets with <", minsetsize, "or >", 
            maxsetsize, "genes\n"))
  if(!is.na(minsetsize) | !is.na(maxsetsize)){
    setN <- set.obj[, .N, by=setID][N>=minsetsize & N<=maxsetsize, setID]
    set.info <- set.info[setID %in% setN]
    set.obj <- set.obj[setID %in% setN]
  }
  # Removing and orphan genes
  # JD: better to keep full obj.info list?
  # added for now parameter obj.in.set (default: T)
  if (obj.in.set) obj.info <- obj.info[objID %in% set.obj[, objID]]
  
  # Create new field with setName to distinguish sets with the same name
  dups <- set.info[, .N, by=setName][N>1, setName]
  if(length(dups)>0){
    cat("Relabeling duplicated gene set names\n")
    set.info[setName %in% dups, setName:=paste0(setName, ": set", 1:.N), 
             by=setName]
  }
  
  # Add in setIDs for unmerged sets
  set.obj[is.na(setID.merged), setID.merged:=setID]
  
  # Record full gene set size (including genes with missing data)
  set.info <- data.table::merge.data.table(set.info, set.obj[, .N, by=setID], 
                                           by="setID")
  
  #------------------------------------------------------#
  ## Prepare input for permutations
  #------------------------------------------------------#
  
  # Remove genes with no information
  no.scores <- obj.info[is.na(objStat), objID]
  obj.info <- obj.info[!(objID %in% no.scores)]
  
  # Create unique integer for each chromosome
  obj.info[, chr.orig:=chr]
  obj.info[, chr:=as.integer(factor(sort(chr)))]
  
  # Reorder gene IDs
  obj.info[, objID.orig:=objID]
  obj.info[order(chr, startpos, endpos, objName), objID:=1:.N]
  
  # Reorder set IDs
  set.info[, setID.orig:=setID]
  set.info[order(setName), setID:=1:.N]
  
  # Regenerate set x object table
  set.obj.tmp <- data.table::merge.data.table(set.obj[, .(setID.orig=setID, 
                                                          objID.orig=objID)], 
                                              obj.info[, .(objID, objID.orig)], 
                                              by="objID.orig")
  
  set.obj <- data.table::merge.data.table(set.obj.tmp, 
                                          set.info[, .(setID, setID.orig)], 
                                          by="setID.orig")
  
  cat(paste("Finished data entry:", nrow(set.info), "gene sets and", 
            set.obj[, uniqueN(objID)], "genes remaining\n"))
  
  return(list(set.info=set.info[order(setID)], obj.info=obj.info[order(objID)],
              set.obj=set.obj[order(setID, objID), .(setID, objID, setID.orig, 
                                                     objID.orig)], 
              minsetsize=minsetsize, population=population))
}



#===========================================================================
# Function: merge.similar.sets(set.info, set.obj)
# Merge gene sets that have more than 95% similarity
#
# -SI : dataframe with fields setID, setName, ...
# -SO : dataframe with fields setID, objID
# -min.sim : float defining the minimum proportion of gene sharing before merging; default = 0.95
# -n.cores : integer specifying number of cores to use; default = 'default' (max cores - 1)
#===========================================================================

merge.similar.sets <- function(SI, SO, min.sim=0.95, n.cores='default', bfs=T){

  # Get similarity matrix
  # Which sets are min.sim proportion similar (two way)
  # Choose the one with largest original set to keep
  # Remove rest, but keep link in set.obj.merged
  
  # JD: added pc.bfs function to find connected components
  pc.bfs <- function(graph){
    # Implementation of Breadth-First-Search (BFS) using adjacency matrix
    # to find connected components in graph
    # Returns list with components, name of each element is first node
    pc<-list()
    cc<-0
    ix.do<-which(rowSums(graph)>1)
    visited = rep(FALSE, nrow(graph))
    for (k in ix.do){
      #for (k in seq_along(visited)){
      if (!visited[k]){
        k.char<-as.character(k)
        queue <- c(k)
        visited[k]<-T
        cc<-cc+1
        pc[[cc]]<-c(k)
        # While there are nodes left to visit...
        while(length(queue) > 0) {
          node = queue[1] # get...
          queue = queue[-1] # ...and remove next node
          ix<-graph[node,] & !visited
          add2q <- which(ix)
          pc[[cc]]<-c(pc[[cc]],add2q)
          visited <- visited | ix
          queue <- c(queue,add2q)
        }
      }
    }
    return (pc)
  }
  
  #parallize processing
  library(doFuture)
  N.cores <- future::availableCores() 
  
  if(n.cores %in% c("max", "MAX", "Max", "default", "Default", "DEFAULT") | is.numeric(n.cores)){
    if(n.cores %in% c("max", "MAX", "Max")){
      n.cores <- N.cores
      cat(paste0("Using all available (", N.cores, ") cores\n"))
    }
    if(n.cores %in% c("default", "Default", "DEFAULT")){
      n.cores <- N.cores-1
      cat(paste("Using", N.cores-1, "cores\n"))
    }
    if(is.numeric(n.cores)){
      n.cores <- as.integer(n.cores)
      if(n.cores<0 | n.cores>N.cores){
        n.cores <- N.cores-1
        e.mss <- paste("No. cores misspecified: using", n.cores, "cores")
        warning(e.mss, immediate.=T)
      }
    }
  }else{
    n.cores <- N.cores-1
    e.mss <- paste("No. cores misspecified: using", n.cores, "cores")
    warning(e.mss, immediate.=T)
  }
  
  doFuture::registerDoFuture()
  future::plan(multiprocess, workers=n.cores) #uses multicore if available
  
  #generate pathway x gene matrix
  
  # JD: assumption is made that setID and objID are integers, 
  # which is not always the case (e.g. Ensembl gene IDs)
  # better to convert to integers 1:N
  # To avoid errors: set.obj should have setIDs without 'gaps'
  # otherwise diag will contain zeros and division by diag will give NAs
  # ===========================================================================

  # JD: using factors
  SO$setID.f<-as.factor(SO$setID)
  SO$objID.f<-as.factor(SO$objID)
  PxG <- Matrix::sparseMatrix(i=as.integer(SO$setID.f), 
                              j=as.integer(SO$objID.f), x=1L)
  # done using factors
  # ===================================================
  
  #PxG <- Matrix::sparseMatrix(i=SO$setID, j=SO$objID, x=1L)
  NP <- nrow(PxG)
  
  #compute proportion shared genes
  sim.mat <- tcrossprod(PxG)

  # JD: I got an error with the following line: 'problem too large' 
  # problem solved when using factors (see above)
  p.mat <- (sim.mat/diag(sim.mat)) > min.sim
  
  # IJ.pairs <- which((triu(t(p.mat)) + triu(p.mat)) == 2)
  # set.m <- data.table(P=IJ.pairs %% NP, 
  #                     X=ceiling(IJ.pairs / NP))
  # set.m[P==0, P:=NP]
  
  # if (! bfs){
  #   # JD: better to use 'arr.ind':
  #   set.m<-data.table(which(triu(t(p.mat)) + triu(p.mat) == 2, arr.ind = T))
  #   colnames(set.m)<-c("P","X")
  #   set.m <- set.m[P!=X] # remove diagonal
  #   
  #   #concatenate pathways
  #   pc <- list()
  #   cc=0
  #   for(i in set.m[, unique(P)]){
  #     PC <- c(i, set.m[P==i, unique(X)])
  #     if(!any(PC %in% unlist(pc, use.names=F))){
  #       cc=cc+1
  #       pc[[cc]] <- PC
  #     }
  #   }
  # } else {
    #alternative: bfs search
    adj.mat <- t(p.mat) + p.mat == 2
    pc <- pc.bfs(adj.mat)
  #}
  cat(paste0("Quantified gene overlap in ", NP, " gene sets\n"))
  
  # ============================================
  # JD: if using factors, convert them back now
  # 
  for (i in seq_along(pc)){
    if (is.numeric(SO$setID)){
      pc[[i]]<-as.numeric(levels(SO$setID.f)[pc[[i]]])
    } else {
      pc[[i]]<-levels(SO$setID.f)[pc[[i]]]
    }
  } 
  SO$setID.f<-NULL
  SO$objID.f<-NULL
  #=============================================
  
  if(length(pc)>0){
    so.new <- foreach(i=pc, .combine=rbind,
                      .export=c("SO", "new.setID")) %dopar% {
      so.now <- SO[setID %in% i]
      #rename according to largest set; take first if equivalent sized sets
      # JD: if we would sort on set size in the beginning taking max would not be 
      # necessary
      new.setID <- so.now[, .N, by=setID][N==max(N), setID][1]
      data.table(setID=new.setID, objID=so.now[, unique(objID)], 
                 setID.merged=so.now[, paste0(unique(setID), collapse=",")])
    }

    SO.out <- rbind(SO[!(setID %in% unlist(pc, use.names=F)), 
                       .(setID, objID, setID.merged=NA)], 
                    so.new)[order(setID, objID)]
    SI.out <- SI[setID %in% SO.out$setID]
    SI.out <- SI.out[setID %in% unlist(pc, use.names=F), 
                     setName:=paste0(setName, "*")]
    SI.merged <- so.new[, .(setID.merged=unique(setID.merged)), by=setID]
    
    cat(paste0("Merged ", length(pc), " gene sets >", min.sim, " similarity\n"))
  }else{
    SO.out <- SO
    SI.out <- SI
    SI.merged <- NA
    cat(paste0("No gene sets >", min.sim, "similarity\n"))
  }
  
  return(list(set.obj=SO.out, set.info=SI.out, set.obj.merged=SI.merged))
}


#===========================================================================
# Function: polylinkr(set.info, obj.info, set.obj, n.cores=NA, linked=T,
#                     minsetsize=10, n.perm.p=1000000, perm.block.size=1000,
#                     prune=T, n.perm.pruning=10000, n.pruned.sets=100,
#                     n.FDR=100, FDR.block.size=10, perm.mat=NULL)
# core function used generate null distribution, perform pruning 
# and calculate corrected p values
#
# - set.info: data.frame or data.table file with fields:
#                  setID, setName, ...
# - set.obj : data.frame or data.table with fields:
#                  setID, objID
# - obj.info : data.frame or data.table with fields:
#                  objID, objName, objStat, (objBin), (objSNPcnt), ...
# - n.cores : character or integer specifying the number of cores to use; 
#   options include 'max' = all cores, 'default' = max cores - 1, 
#   or integer value; default = 'default'
# - linked : logical indicating if linkage or non-linkage based permutations 
#   should be used to compute the null; default = TRUE
# - n.perm.p : integer indicating the number of iterations of permutation 
#   algorithm to use to compute null used in p-value estimation; default = 100000
# - perm.block.size : integer defining the number of iterations on which to 
#   summarize results during p-value estimation (improves computation speed); 
#   default = 1000
# - prune : logical indicating if pruning step should be used; default = FALSE; 
#   if TRUE is selected, tests will be performed on the complete complement of 
#   genes for each gene set
# - n.perm.pruning : if pruning is used, integer defining the number of 
#   iterations during pruning step which to summarize results. if using a 
#   limited number of processors (<10), this should be smaller than the amount 
#   used to calculate the p-values (i.e. n.perm.p); default = 10000, ignored if 
#   prune == FALSE
# - n.pruned.sets : if pruning is used, integer specifying the subset of 
#   number of gene sets to compute p-values for; default = 100, 
#   ignored if prune == FALSE
# - n.FDR : if pruning is used, number of random pruning steps to run to 
#   calculate the FDR-corrected p-values (note that the number of permuted 
#   set scores used in the FDR correction is n.FDR * n.pruned.sets); 
#   default = 100, ignored if prune == FALSE
# - precise.p.method : if pruning is used, determines the method to used for 
#   null - options: 'full' = all genes, 'pruned' = pruned genes removed, 
#   'both' = both methods used; default = 'both'
# - perm.mat : optional permutation matrix values provided by user; default = NULL
#===========================================================================

polylinkr <- function(set.info, obj.info, set.obj, n.cores="default", linked=T, 
                      minsetsize=10, n.perm.p=100000, perm.block.size=1000, 
                      prune=T, n.perm.pruning=5000, n.FDR=100, n.pruned.sets=50, 
                      precise.p.method="both", est.pi0=TRUE, perm.mat=NULL, 
                      startT=NA)
{
  
  ##===================================================================##
  ## internal functions
  ##===================================================================##
  
  # records cumulative run time at specific points
  get.time <- function(t) { # time conversion
    tt <- as.numeric(difftime(Sys.time(), t, units = "secs"))
    H <- tt %/% 3600
    rr <- tt %% 3600
    if(rr > 0) {
      M <- rr %/% 60
      rr2 <- rr %% 60
      if(rr2 > 0) {
        S <- round(rr2)
      }else{
        S <- 0
      }
    }else{
      M <- 0
      S <- 0
    }
    return(paste0(H, "h ", M, "m ", S, "s"))
  }
  
  # unique linked null permutations
  make.perm.mat <- function(n.chr, n.genes, N){
    chr.ord <- c(t(rowRanks(matrix(sample(1:(n.chr*N)), 
                                   ncol=n.chr))))
    rot <- floor(runif(N, min=2, max=n.genes))
    list(CHR.ORD=chr.ord, ROT=rot)
  }
  
  # enumerate linked null permutations
  permute.data <- function(score.list, rot.list, chr.ord.now, rot.now){
    all.perm.chr <- unlist(score.list[chr.ord.now], use.names=F)
    all.rot <- unlist(rot.list[rot.now], use.names=F)
    all.perm.chr[all.rot]
  }
  
  # improve runtime by separating iterations into fixed-sized blocks
  get.blocks <- function(n.perm, block.size, n.chr) {
    if(n.perm < block.size){
       rb <- cbind(1, n.perm)
       cob <- cbind(1, n.perm*n.chr)
    } else {
       ss <- unique(c(seq(0, n.perm, block.size), n.perm))
       rb <- cbind(ss[1:(length(ss)-1)]+1, ss[2:length(ss)])
       cob  <- cbind(ss[1:(length(ss)-1)]*n.chr+1, 
                     ss[2:length(ss)]*n.chr)
    }
    return(list(RB=rb, COB=cob, N=nrow(cob)))
  }
  
  # update parameters during pruning step
  pruning.param <- function(set.obj.now, top.sets.now, 
                            nobj=n.genes, nset=n.sets){
    #get genes from top ranked set
    gene.rm <- set.obj.now[setID %in% top.sets.now, objID] 
    set.obj.now[objID %in% gene.rm, KEEP:=0]
    #update removed gene list
    obj.out<- set.obj.now[KEEP==0, sort(unique(objID))]
    #determine sets with too few genes
    set.out <- set.obj.now[, .(N=sum(KEEP)), by=setID][N<minsetsize, unique(setID)]
    return(list(so.now=set.obj.now, 
                set.n.remaining=nset-length(set.out),
                obj.n.remaining=nobj-length(obj.out), 
                set.out=set.out, obj.out=obj.out))
  }
  
  # estimate pi0 using histogram method
  pi0.estimator <- function(pp, n.EXP, n.OBS, n.bins, tolerance=10^-4){
    #create percentile bins
    p.bins <- pp[, seq(min(P), max(P), length.out=n.bins+1)] 
    pp[, P.bin:=cut(P, p.bins, include.lowest=T)]
    #calculate observed and expected values
    fdr.dt <- data.table::as.data.table(cbind(pp[REP==1, table(P.bin)], 
                                              pp[REP!=1, table(P.bin)]),
                                        keep.rownames=TRUE)
    data.table::setnames(fdr.dt, c("P.bin", "OBS", "EXP"))
    fdr.dt[, EXP.scaled:=(EXP/n.EXP)*n.OBS]
    obs.excess <- fdr.dt[, which(cumsum(OBS < EXP.scaled)==1)-1][1] # determine obs > exp up to first failure
    
    pi0 <- 1 # set initial condition
    if(obs.excess>0){ # calculate pi0 using histrogram method
      #estimate pi0 using iterative procedure
      track.pi0 <- pi0
      k=0; TT=TRUE
      while(isTRUE(TT)){ #run to covergence
        k=k+1
        true.pos <- fdr.dt[1:obs.excess, sum(OBS-EXP.scaled)/n.OBS]
        pi0 <- 1 - true.pos
        fdr.dt[, EXP.scaled:=pi0*EXP.scaled]
        obs.excess <- fdr.dt[, which(cumsum(OBS < EXP.scaled)==1)-1][1]
        track.pi0 <- c(track.pi0, pi0)
        # test convergence
        TT <- abs(track.pi0[k] - pi0) > tolerance
      }
    }
    
    return(pi0)
  }
    
  # calculate q values
  q.estimator <- function(pp, pi0, n.EXP, n.OBS){ 
    # cumulative number of rejected hypotheses (account for non-unique p-values)
    RP.star <- pp[REP==1, cumsum(table(P))]
    # bins to evaluate expected number of true nulls
    emp.p.bins <- pp[REP==1, c(0, sort(unique(P)))] 
    # count expected true nulls if all nulls are true
    VP.star <- cumsum(table(pp[REP!=1, cut(P, emp.p.bins, right=TRUE)]))
    # compute expected true nulls after scaling by pi0 and number of tests
    VP.star.scaled <- pi0*(VP.star/n.EXP)*n.OBS 
    FDR.vect <- VP.star.scaled / RP.star
    FDR.vect[FDR.vect==0] <- 1/n.EXP # reset values == 0
    FDR.vect[FDR.vect>1] <- 1 # reset values larger than 1
    
    #compute q values from FDR values
    M <- min(FDR.vect)
    I <- max(which(FDR.vect==M))
    K <- 1
    q.vect <- c(rep(M, I), rep(NA, length(FDR.vect)-I))
    while(I<length(FDR.vect)){# reset FDRs for larger p values with smaller FDR
      K <- I+1
      M <- min(FDR.vect[K:length(FDR.vect)])
      I <- max(which(FDR.vect==M))
      q.vect[K:I] <- M
    }
    
    # modify q value vector to include repeated p values
    if(length(q.vect)<n.OBS){ 
      q.it <- pp[REP==1, unname(table(P))]
      q.vect <- unlist(lapply(1:length(q.it), function(x) rep(q.vect[x], q.it[x])))
    }
    
    #output
    return(list(q=data.table::data.table(pp[REP==1][order(P), .(RANK, P)], 
                                         Q=q.vect), pi0=pi0))
  }


  ##===================================================================##
  ##PART 1: clean data and run checks
  ##===================================================================##
  
  # Overall starting time
  if(is.na(startT)){
    startT <- Sys.time()
  }
  cat("[STEP 1] Performing data checks...\n")
  
  if(!is.data.table(set.info)) data.table::setDT(set.info)
  if(!is.data.table(obj.info)) data.table::setDT(obj.info)
  if(!is.data.table(set.obj)) data.table::setDT(set.obj)
  
  # add in mandatory column if data is not read in by read.setobj.tables
  if(!("setID.orig" %in% names(set.info))){
    set.info[, setID.orig:=setID]
    set.obj[, setID.orig:=setID]
  }
  if(!("objID.orig" %in% names(obj.info))){
    obj.info[, objID.orig:=objID]
    set.obj[, objID.orig:=objID]
  }
  
  # set scores to integers (reduce memory pressure and speed up calculations)
  if(obj.info[, !is.integer(objStat)]){
    max.sum.int <- .Machine$integer.max
    max.sum.fl <- sum(obj.info[order(-objStat),objStat][1:set.info[, max(N)]])
    #JD: better floor instead of round
    prec <- floor(log10(max.sum.int / max.sum.fl))
    obj.info[, objStat.orig:=objStat]
    obj.info[, objStat:=as.integer(round(objStat*10^prec), 0)]
  }else{
    prec <- NA
  }
  
  #get no. paths and genes
  n.genes <- obj.info[, .N]
  n.sets <- set.info[, .N]
  
  #---------------------------------#
  ## check iteration options
  #---------------------------------#
  
  # check for violations in parameter settings
  n.check <- sapply(list(n.perm.p, n.FDR, perm.block.size, n.perm.pruning,
                         n.pruned.sets), function(x) is.numeric(x) & x>0)
  
  if(!all(n.check)){
    num.param <- c("n.perm.p", "n.FDR", "perm.block.size", 
                   "n.perm.pruning", "n.pruned.sets")[which(!n.check)]
    stop(paste("Operation halted:", paste(num.param, collapse=" & "), 
               ifelse(length(num.param)==1, "is", "are"),
               "non-numeric or less than 0. Re-enter or use defaults."))
  }
  
  # permutation parameters
  if(prune){ #FDR computation (if selected)
    if(n.perm.pruning > n.perm.p){
      n.perm.pruning <- n.perm.p
      warning(paste0("Number of pruning iterations exceeds null iterations ",
                     "for p-value estimation. Resetting to ", n.perm.pruning, 
                     " iterations"), immediate.=T)
    }
    
    if(n.pruned.sets > n.sets){
      n.pruned.sets <- n.sets
      warning(paste0("Number of gene sets evaluated in pruning step exceeds ",
                     "total number of gene sets. Resetting to ", n.pruned.sets,
                     " iterations"), immediate.=T)
    }
    
    # in case non-integer values are used
    n.FDR <- as.integer(n.FDR)
    n.perm.pruning <- as.integer(n.perm.pruning)
    n.pruned.sets <- as.integer(n.pruned.sets)
  }
  
  if(perm.block.size <= 0 | perm.block.size > n.perm.p | 
     perm.block.size > n.perm.pruning){
    perm.block.size <- min(n.perm.pruning, 1000)
    warning(paste0("null evalution block size incorrectly entered or ",
                   "size exceeds total number of iterations"), immediate.=T)
    warning(paste0("initial null evalution set to blocks of ", perm.block.size,
                  " iterations"), immediate.=T)
  }
  
  # in case non-integer values are used
  n.perm.p <- round(n.perm.p, 0)
  perm.block.size <- round(perm.block.size, 0)
  
  #---------------------------------#
  ## set up parallel back-end
  #---------------------------------#
  
  if(prune){ # allow parallel processing if pruning used
    library(doFuture)
    
    #set up permutation block size
    N.cores <- future::availableCores()

    if(n.cores %in% c("max", "MAX", "Max", "default", "Default", "DEFAULT") | 
       is.numeric(n.cores)){
      if(n.cores %in% c("max", "MAX", "Max")){
        n.cores <- N.cores
        cat(paste0("Using all available (", N.cores, ") cores\n"))
      }
      if(n.cores %in% c("default", "Default", "DEFAULT")){
        n.cores <- N.cores-1
        cat(paste("Using", N.cores-1, "cores\n"))
      }
      if(is.numeric(n.cores)){
        n.cores <- as.integer(n.cores)
        if(n.cores<0 | n.cores>N.cores){
          n.cores <- N.cores-1
          e.mss <- paste("No. cores misspecified: using", n.cores, "cores")
          warning(e.mss, immediate.=T)
        }
      }
    } else {
      n.cores <- N.cores-1
      e.mss <- paste("No. cores misspecified: using", n.cores, "cores")
      warning(e.mss, immediate.=T)
    }
  
    doFuture::registerDoFuture()
    #uses multicore if available
    future::plan(multiprocess, workers=n.cores) 
    #increase default maximum size of objects exported to functions
    #JD: I propose to increase the maxsize to 1500
    #options(future.globals.maxSize=1000*1024^2) 
    options(future.globals.maxSize=1500*1024^2) 
  }
  
  cat("Completed checking data\n")
  cat(paste("Passed", n.sets, "gene sets &", n.genes, "unique genes\n"))
  if(n.genes-set.obj[, uniqueN(objID)]!=0){
    cat(paste0("NOTE: ", n.genes-set.obj[, uniqueN(objID)], " genes are not ",
               "in any gene set, but will be retained for computing null\n"))
  }
  cat(paste0(" *** Time elapsed: ", get.time(startT), " *** \n\n"))
  
  
  ##===================================================================##
  # PART 2: Set up permutations
  ##===================================================================##
  
  if(prune){
    cat("[STEP 2] Pruning requested: preparing inputs for permutation tests",
        "and FDR correction...\n")
  }else{
    cat("[STEP 2] Pruning not requested: running enrichment test...\n")
  }
  
  # generate pathway x gene matrix
  PxG <- Matrix::sparseMatrix(i=set.obj$setID, j=set.obj$objID, 
                              x=1L, dims=c(n.sets, n.genes))
  
  #----------------------------------------------------------#
  ## generate observed and FDR scores (if pruning chosen)
  #----------------------------------------------------------#

  # calculate set scores for observed data and for FDR correction
  # If no pruning only generate observed scores
  # If pruning is used compute single matrix of observed and FDR scores
  
  # generate parameters for blocks of null computations
  n.chr <- obj.info[, length(unique(chr))]
  run.blocks <- get.blocks(n.perm=ifelse(prune, n.perm.pruning, n.perm.p), 
                           perm.block.size, n.chr)
  
  if(linked){
    # if using sampling with linked genes
    # generate list of gene positions and scores by chromosome, 
    # and unique gene rotations
    pos.list <- split(obj.info[order(chr, startpos, endpos), objID], 
                        sort(obj.info$chr))
    score.list <- split(obj.info[order(chr, startpos, endpos), objStat], 
                        sort(obj.info$chr))
    rot.list <- append(list(1:n.genes),
                       lapply(2:n.genes, function(x) c(x:n.genes, 1:(x-1))))
    if(prune){
      # create additional perms to account for any redundancy
      fdr.mat <- make.perm.mat(n.chr, n.genes, n.FDR) 
      # observed values + randomised values
      obs.fdr.mat <- list(CHR.ORD=c(1:n.chr, fdr.mat$CHR.ORD),
                          ROT=c(1, fdr.mat$ROT)) 
    } else {
      obs.fdr.mat <- list(CHR.ORD=1:n.chr, ROT=1)
    }
    PRM <- permute.data(score.list, rot.list,
                        chr.ord.now=obs.fdr.mat$CHR.ORD, 
                        rot.now=obs.fdr.mat$ROT)# permuted gene scores
  } else {
    # generate fully random gene score sampling
    library(dqrng)
    dqrng::dqRNGkind("Xoroshiro128+") # fast random sampling
    fdr.mat <- dqrng::dqsample.int(.Machine$integer.max, 1) # set random seed
    dqrng::dqset.seed(fdr.mat)
    if(prune){
      obs.fdr.mat <- c(seq.int(1, n.genes), # observed values + randomised FDR values
                       dqrng::dqsample.int(n.genes*n.FDR) %/% 
                       as.integer(n.FDR)) # create additional perms to account for any reduncancy
      obs.fdr.mat[obs.fdr.mat==0] <- as.integer(n.genes) # permuted gene positions
    } else {
      obs.fdr.mat <- seq.int(1, n.genes) # permuted gene positions
    }
    PRM <- obj.info$objStat[obs.fdr.mat] # permuted gene scores
  }
  
  # observed set score matrix: 1st entry = observed, 
  # remainder for FDR estimation (if pruning used)
  SS <- matrix(PRM, ncol=ifelse(prune, n.FDR+1, 1))
  m.obs <- matrix(as.integer((PxG %*% SS)@x), 
                  ncol=ifelse(prune, n.FDR+1, 1)) 
  rm(PRM, obs.fdr.mat) # clean up large objects not subsequently used
  gc()
  
  #---------------------------------#
  ## generate null scores
  #---------------------------------#
  
  # generate full permutation matrix for null computations
  if(is.null(perm.mat)){ # generate null parameters
    if(linked){
      cat("Generating null scores using linkage-based permutations\n")
      perm.mat <- make.perm.mat(n.chr, n.genes, N=n.perm.p)
    } else {
      cat("Generating null scores using standard permutations (no linkage)\n")
      #keep seed to replicate sampling
      perm.mat <- dqsample.int(.Machine$integer.max, 
                               n.perm.p %/% perm.block.size) 
    }
  } else {
    cat("Permutations performed using user-specified parameters\n")
    if(linked){ # check that user matrix matches input data
      if(n.perm.p != length(perm.mat$ROT)){
        stop("Operation halted: user-specified parameters are not compatible",
             " with input data\n")
      }
    } else {
      if(n.perm.p %/% perm.block.size != length(perm.mat)){
        stop("Operation halted: user permutation parameters are not compatible",
             " with input data\n")
      }
    }
  }
  
  # Compute permuted scores
  # NOTE: only a subset of the full permutation matrix is used if pruning step 
  # performed
  if(prune){
    PRM.pos <- list()
    PRM <- list()
  } else {
    cat(paste("Computing enrichment tests using", n.perm.p, "permutations\n"))
    # JD: maybe a bit faster: initialize failed.test as vector(#tests)?
    #failed.tests <- integer(n.sets)
    failed.tests <- 0 # counter to collect number of failed tests
  }
  
  progressr::with_progress({ # allow updating during processing
    prog <- progressr::progressor(steps=run.blocks$N)
    for(l in 1:run.blocks$N){
      prog(sprintf("STEP=%i", l))
      if(linked){ # linkage-based permutations
        rot.now <- run.blocks$RB[l, ]
        cob.now <- run.blocks$COB[l, ]
        rot.now.i <- perm.mat$ROT[rot.now[1]:rot.now[2]]
        cob.now.i <- perm.mat$CHR.ORD[cob.now[1]:cob.now[2]]
        if(prune){ # keep gene positions for pruning step
          PRM.pos[[l]] <- permute.data(pos.list, rot.list,
                                       chr.ord.now=cob.now.i, 
                                       rot.now=rot.now.i) # permuted gene positions
          PRM[[l]] <- permute.data(score.list, rot.list,
                                   chr.ord.now=cob.now.i, 
                                   rot.now=rot.now.i)
        } else { # estimate p values on the fly
          score.mat <- PxG %*% matrix(permute.data(score.list, rot.list, 
                                                   chr.ord.now=cob.now.i, 
                                                   rot.now=rot.now.i), 
                                      ncol=perm.block.size)
          failed.tests <- failed.tests + rowSums(score.mat > c(m.obs))
        }
      } else { # standard permutations
        dqset.seed(perm.mat[l])
        PRM.pos.temp <- dqsample.int(n.genes*perm.block.size) %/% 
                          as.integer(perm.block.size)
        PRM.pos.temp[PRM.pos.temp==0] <- as.integer(n.genes)
        if(prune){ # keep gene positions for pruning step
          PRM.pos[[l]] <- PRM.pos.temp
          PRM[[l]] <- obj.info$objStat[PRM.pos.temp]
        } else { # estimate p values
          score.mat <- PxG %*% matrix(obj.info$objStat[PRM.pos.temp], 
                                      ncol=perm.block.size)
          failed.tests <- failed.tests + rowSums(score.mat > c(m.obs)) 
        }
        rm(PRM.pos.temp) # clean up
        gc()
      }
    }
  })
  if(prune){
    cat("Completed set score computation\n")
  }else{
    cat("Completed gene set enrichment tests\n")
  }
  
  ##===================================================================##
  # PART 3: determine top ranking pathways for observed and FDR
  ##===================================================================##
  
  if(prune){
    #generate null score matrix (set scores x permutations)
    score.mat <- matrix((PxG %*% matrix(unlist(PRM, use.names=F), 
                                        ncol=n.perm.pruning))@x, 
                        ncol=n.perm.pruning)
  
    #determine top gene sets for pruning step
    rr1 <- matrixStats::rowRanks(-m.obs)
    rr2 <- matrixStats::rowRanks(cbind(-m.obs, -score.mat))[, 1:(n.FDR+1)]
    p.vals <- (rr2 - rr1 + 1) / (n.perm.pruning+1)
    rm(rr1, rr2)
    gc()
    
    #if equal highest rank, keep pathway with most genes
    min.p <- matrixStats::colMins(p.vals)
    top.sets <- foreach(i=1:length(min.p), .combine=c) %do% {
      pk <- which(p.vals[, i] == min.p[i])
      if(length(pk)>1){
        pk <- set.info[setID %in% pk][order(-N), setID][1]
      }
      pk
    }

    #convert permutation lists to vectors
    PRM <- unlist(PRM, use.names=F)
    PRM.pos <- unlist(PRM.pos, use.names=F)
    
  } else {
    p.vect <- (failed.tests + 1) / (n.perm.p + 1)
    q.vect <- qvalue::qvalue(p.vect)
    sio <- data.table::data.table(set.info[order(setID),
                                           .(setName, setID=setID.orig, N)], 
                       N.pruned=NA, 
                       setScore=c(m.obs)/ifelse(is.na(prec), 1, 10^prec), 
                       setScore.pruned=NA, P=p.vect, 
                       Q=q.vect$qvalues)[order(P)]
    OUT <- list(set.info=sio, set.obj=set.obj, pi0.est=q.vect$pi0, 
                null.set=NA, permutation.mat=perm.mat, 
                permutation.type=c("random", "linked")[as.numeric(linked)+1])
  }
  cat(paste0(" *** Time elapsed: ", get.time(startT), " *** \n\n"))
  
  
  if(prune){ # run pruning step if requested
    ##===================================================================##
    # PART 4: run pruning step
    ##===================================================================##
  
    cat("[STEP 3] Running pruning step (this may take some time)...\n")
    cat(paste("Generating", n.perm.pruning, "null permutations for", 
              n.FDR, "FDR iterations and observed data\n"))
    cat(paste("The top", n.pruned.sets, 
              "gene sets will be evaluated in each case\n"))
    
    with_progress({ # allow updating during parallel processing
      prog <- progressor(steps=n.FDR+1)
      pruned.sets <- foreach::foreach(f=1:(n.FDR+1),
                                      .export=c("m.obs", "p.vals", "top.sets", 
                                                "set.obj", "perm.mat", "PRM", 
                                                "PRM.pos", "SS")) %dopar% {
        
        #get input datasets
        ts.now <- top.sets[f]
        
        #set up dataset for pruning 
        SP <- pruning.param(set.obj.now=data.table(set.obj, KEEP=1L), 
                            top.sets.now=ts.now)
        
        #run pruning step until criteria met
        #i.e. more than 1 path remaining and up to n.pruned.sets iterations
        set.n.remaining <- SP$set.n.remaining
        I=1L
        
        #keep track of top sets and their genes
        top.sets.dt <- data.table::data.table(REP=f, RANK=I,
                                              set.info[setID %in% ts.now, 
                                                       .(setName, setID, 
                                                         setID.orig, N, 
                                                         N.pruned=N)], 
                                              setScore=m.obs[ts.now, f], 
                                              setScore.pruned=m.obs[ts.now, f],
                                              P.raw=p.vals[ts.now, f], 
                                              P.raw.pruned=p.vals[ts.now, f],
                                              P.full.null=NA, P.pruned.null=NA)
        
        set.obj.dt <- SP$so.now[setID==ts.now, .(setID, objID, RANK=1)]
        
        while(set.n.remaining>1 & I<n.pruned.sets){ # iterate until criteria not met
          I=I+1 # update iteration counter
          
          #create new PxG matrix removing unused genes
          PxG.now <- PxG %*% Matrix::Diagonal(n=n.genes, x=1L)[, -SP$obj.out]
          #create new observed scores
          m.obs.now <- PxG.now %*% SS[-SP$obj.out, f]
          
          #generate updated set scores
          PRM.now <- PRM[!(PRM.pos %in% SP$obj.out)]
          score.mat.now <- PxG.now %*% matrix(PRM.now, nrow=SP$obj.n.remaining)
          
          #determine top path and update top set datatable
          p.now <- (rowSums(score.mat.now > m.obs.now@x)+1) / (n.perm.pruning+1)
          p.rank <- rank(p.now, ties.method="average")
          ts.now <- setdiff(which(p.rank == min(p.rank[-SP$set.out])), 
                            SP$set.out)
          #if more than 1 top set, take largest set
          if(length(ts.now)>1){
            top.sets.now <- SP$so.now[setID %in% ts.now, .(N=sum(KEEP)), 
                                      by=setID]
            ts.now <- top.sets.now[N==max(N), setID][1]
          }
          #update top sets and their genes
          tp.dt.now <- data.table::data.table(REP=f, RANK=I, 
                                              set.info[setID %in% ts.now, .(setName, setID, setID.orig, N)], 
                                              N.pruned=SP$so.now[setID %in% ts.now, sum(KEEP)], 
                                              setScore=m.obs[ts.now, f], 
                                              setScore.pruned=m.obs.now[ts.now],
                                              P.raw=p.vals[ts.now, f], 
                                              P.raw.pruned=p.now[ts.now],
                                              P.full.null=NA, P.pruned.null=NA)
          top.sets.dt <- rbind(top.sets.dt, tp.dt.now)
          set.obj.dt <- rbind(set.obj.dt,
                              SP$so.now[setID==ts.now & KEEP==1, 
                                        .(setID, objID,  RANK=I)])
          
          #update loop criteria and parameters
          SP <- pruning.param(set.obj.now=SP$so.now, top.sets.now=ts.now)
          set.n.remaining <- SP$set.n.remaining
        }
        #update user
        prog(sprintf("STEP=%i", f))

        list(set.info.pruned=top.sets.dt, set.obj.pruned=set.obj.dt)
      }
    })
    
    #clean up
    rm(PRM, PRM.pos, PxG, p.vals, m.obs, fdr.mat, score.mat)
    if(linked & !(precise.p.method %in% c("full", "both"))){
      rm(rot.list, score.list)
    }
    gc()
    cat("Completed pruning step\n")
    cat(paste0(" *** Time elapsed: ", get.time(startT), " *** \n\n"))
  
    # only recalculate if higher precision requested
    if(n.perm.p > n.perm.pruning){ 
      ##===================================================================##
      # PART 5: calculate high precision p-values
      ##===================================================================##
      
      cat(paste("[STEP 4] Recalculating p-values using", n.perm.p, 
                "iterations...\n"))
      
      ##===================================================================##
      # PART 5.1: calculate p-values using all genes in null
      ##===================================================================##
      
      if(precise.p.method %in% c("full", "both")){
        
        cat(paste("Computing null with all genes\n"))
        
        #compute multiple iterations in blocks
        tot.pruned.sets <- (n.FDR+1)*n.pruned.sets
        if(tot.pruned.sets > 1000){
          block.size <- 1000 %/% n.pruned.sets # compute 1000 gene sets per block
          block.bins <- cut(1:(n.FDR+1), seq(0, n.FDR+1, by=block.size))
          block.bins[is.na(block.bins)] <- rev(levels(block.bins))[1] # include unbinned iterations into last bin
          blocks <- split(1:(n.FDR+1), block.bins) # iterations in each block
        }else{ # compute observed and FDR iterations in a single block
          blocks <- list(1:(n.FDR+1))
        }
        run.blocks <- get.blocks(n.perm.p, perm.block.size, n.chr) #set up block computation for high precision p-values
        
        with_progress({ # allow updating during parallel processing
          prog <- progressor(steps=length(pruned.sets))
          #calculate p values in blocks
          p.vals.full <- foreach::foreach(z=blocks, .combine=rbind,
                                          .export=c("pruned.sets", "n.perm.p", 
                                                    "perm.mat", "n.genes", "SS",
                                                    "score.list", "rot.list")) %dopar% {
                                        
            pp.now <- pruned.sets[z] # get block of results
            
            #-----------------------------------#
            ## generate set matrices
            #-----------------------------------#
            
            block.size.now <- length(pp.now)
            set.info.now <- NULL
            set.obj.now <- NULL
            for(I in pp.now){
              sin <- I$set.info.pruned
              set.info.now <- rbind(set.info.now, sin)
              son <- I$set.obj.pruned
              son[, REP:=sin[, unique(REP)]]
              set.obj.now <- rbind(set.obj.now, son)
            }
            
            m.obs.now <- set.info.now[order(REP, RANK), setScore.pruned]
            
            #-----------------------------------#
            ## generate pathway x gene matrix
            #-----------------------------------#
            
            # create temporary unique set ID
            set.info.now[, setID.now:=paste0(REP, ".", RANK)] 
            ord <- set.info.now[order(REP, RANK), unique(setID.now)]
            set.info.now[, setID.now:=as.integer(factor(setID.now, levels=ord))]
            set.obj.now <- merge(set.obj.now, set.info.now[, .(REP, RANK, setID.now)], 
                                 by=c("REP", "RANK"))
            n.sets.now <- set.info.now[, uniqueN(setID.now)]
            
            PxG.now <- sparseMatrix(i=set.obj.now$setID.now, j=set.obj.now$objID,
                                    x=1L, dims=c(n.sets.now, n.genes))
            
            #---------------------------------#
            ## generate null scores
            #---------------------------------#
            
            #generate updated set scores -- run on full set of permutations
            score.mat <- foreach(l=1:run.blocks$N, .combine="+") %do% {
              if(linked){
                rr <- run.blocks$RB[l, ]
                cc <- run.blocks$COB[l, ]
                
                cc.n <- perm.mat$CHR.ORD[cc[1]:cc[2]]
                rr.n <- perm.mat$ROT[rr[1]:rr[2]]
                PRM <- permute.data(score.list, rot.list,
                                    chr.ord.now=cc.n, rot.now=rr.n)
              }else{
                dqset.seed(perm.mat[l])
                perm.mat.now <- dqsample.int(n.genes*perm.block.size) %/%
                                as.integer(perm.block.size)
                perm.mat.now[perm.mat.now==0] <- as.integer(n.genes)
                PRM <- obj.info$objStat[perm.mat.now]
              }
              
              sc.m <- PxG.now %*% matrix(PRM, ncol=perm.block.size)
              rowSums(m.obs.now < sc.m)
            }
            #input p-values
            set.info.now[, P.full.null:=(score.mat+1)/(n.perm.p+1)] 
            set.info.now[, setID.now:=NULL]
            #set.obj.now[, setID.now:=NULL]
            
            #update user
            prog(sprintf("STEP=%i", z))
            
            set.info.now
          }
        })
        
        #create output
        set.info.out <- split(p.vals.full, p.vals.full$REP)
        pruned.sets <- foreach::foreach(J=1:(n.FDR+1)) %do% {
          list(set.info.pruned=set.info.out[[J]],
               set.obj.pruned=pruned.sets[[J]]$set.obj.pruned)
        }
        
        if(linked){
          rm(rot.list, score.list) # remove unused large datasets
          gc()
        }
        cat("Completed p-value calculation using all genes for null\n")
        cat(paste0(" *** Time elapsed: ", get.time(startT), " *** \n\n"))
      }
      
      ##===================================================================##
      # PART 5.2 calculate p-values without pruned genes in null
      ##===================================================================##
      
      if(precise.p.method %in% c("pruned", "both")){
        
        cat(paste("Computing null excluding pruned genes\n"))
        
        with_progress({ # allow updating during parallel processing
          prog <- progressor(steps=length(pruned.sets))
          p.vals <- foreach::foreach(z=1:length(pruned.sets), .combine=rbind,
                                     .export=c("pruned.sets", "n.perm.p", 
                                               "perm.mat", "n.genes")) %dopar% {
            
            pp.now <- pruned.sets[[z]] #iterate over observed and FDR ranked sets
            
            #enumerate linked null permutations
            OUT <- foreach(j=1:nrow(pp.now$set.info.pruned)) %do% { #iterate over each set
              so.now <- pp.now$set.obj.pruned[RANK==j]
              si.now <- pp.now$set.info.pruned[j]
  
              #account for reduced genes after pruning
              if(j==1){
                n.genes.remaining <- n.genes
                oi.now <- obj.info[order(chr, startpos, endpos)]
                scores.now <- obj.info$objStat
              }else{
                so.N <- pp.now$set.obj.pruned[RANK %in% 1:(j-1)][, unique(objID)]
                n.genes.remaining <- n.genes - length(so.N)
                oi.now <- obj.info[-so.N][order(chr, startpos, endpos)]
                scores.now <- oi.now$objStat
              }
              #set specific genes
              n.genes.now <- si.now$N.pruned
              
              if(linked){ # using linkage-based permutations
                #create gene pos and chrom order lists specific to current gene set
                oi.now[, chrI:=1:.N, by=chr]
                oi.now[, chrL:=max(chrI), by=chr]
                oi.new <- oi.now[objID %in% so.now$objID]
                oi.new[, chrJ:=as.integer(c(rep(0, .N-1), unique(chrL))), by=chr]
                
                pos.list.now <- split(oi.new$chrI, oi.new$chr)
                chr.list.now <- split(oi.new$chrJ, oi.new$chr)
                if(oi.new[, uniqueN(chr)] < n.chr){ # add back missing chromosomes
                  miss.chr <- oi.now[chr %in% setdiff(1:n.chr, oi.new$chr), 
                                     .(chrL=unique(chrL)), by=chr]
                  for(mm in 1:nrow(miss.chr)){
                    m <- miss.chr[mm, chr]
                    pos.list.now <- append(pos.list.now, list(NA), after=m-1)
                    chr.list.now <- append(chr.list.now, list(miss.chr[mm, chrL]), after=m-1)
                  }
                }
                
                #generate random chromosome assemblies
                pl.now <- unlist(pos.list.now[perm.mat$CHR.ORD], use.names=F)
                cl.now <- unlist(chr.list.now[perm.mat$CHR.ORD], use.names=F)
                
                pl.mat <- matrix(pl.now, ncol=n.perm.p)
                cl.mat <- rbind(rep(0, n.perm.p), matrix(cl.now, ncol=n.perm.p))
                cl.mat <- colCumsums(cl.mat)[-nrow(cl.mat), ]
                pos.mat <- pl.mat + cl.mat 
                
                #generate random rotations
                npm <- nrow(pos.mat)
                PRM.pos.temp <- (pos.mat - matrix(rep(perm.mat$ROT-1, npm), 
                                                  nrow=npm, byrow=T)) %% n.genes.remaining
                PRM.pos.temp[PRM.pos.temp==0] <- n.genes.remaining
                score.mat <- colSums(matrix(scores.now[PRM.pos.temp], 
                                            nrow=npm), na.rm=T)
              }else{ # using standard permutations
                if(z==1){
                  dqset.seed(perm.mat[1]) # reuse first random seed
                }
                PRM.pos.temp <- dqsample.int(n.genes.remaining*n.perm.p, 
                                             size=n.genes.now*n.perm.p) %% 
                                as.integer(n.genes.remaining)
                PRM.pos.temp[PRM.pos.temp==0] <- as.integer(n.genes.remaining)
                score.mat <- colSums(matrix(scores.now[PRM.pos.temp], 
                                            nrow=n.genes.now))
              }
              
              #calculate p-values
              p.now <- (sum(si.now$setScore.pruned < score.mat)+1) / (n.perm.p+1)
              
              #update scores
              si.now[, P.pruned.null:=p.now]
              rm(PRM.pos.temp)
              if(linked){
                rm(score.mat, pos.mat, cl.mat, pl.mat)
              }
              gc()
              si.now
            }
            #update user
            prog(sprintf("STEP=%i", z))
            do.call(rbind, OUT)
          }
        })
        cat("Completed p-value calculation excluding pruned genes from null\n")
        cat(paste0(" *** Time elapsed: ", get.time(startT), " *** \n\n"))
      }
      
    }else{ 
      # if no. of pruning steps is >= no. steps used for precise p values reuse the p values from the pruning step
      e.mss1 <- paste0("No. of iterations for p-value recalculation less than",
                       "iterations used in pruning step\n")
      e.mss2 <- "q values will be calculated using results from pruning step"
      warning(paste0(e.mss1, e.mss2), immediate.=T)
      p.vals <- foreach(ps.now=pruned.sets, .combine=rbind) %do% {
        OUT <- ps.now$set.info.pruned
        OUT[, p.pruned.null:=p.raw.pruned]
      }
    }
    
    
    ##===================================================================##
    # PART 6: calculate q-values
    ##===================================================================##
    
    if(n.perm.p > n.perm.pruning){
      cat("[STEP 5] Computing q-values and preparing output...\n")
    }else{
      cat("[STEP 4] Computing q-values and preparing output...\n")
    }
  
    #create outputs
    set.obj.temp <- data.table::merge.data.table(pruned.sets[[1]]$set.obj.pruned[, .(setID, objID)],
                                                 set.obj, by=c("setID", "objID"))
    set.info.temp <- p.vals[REP==1]
    
    #n.EXP <- ~ n.FDR * n.pruned.sets # total number of null hypothesis tests
    n.EXP <- p.vals[REP!=1, .N]
    #n.OBS <- n.FDR # total number of observed hypothesis tests
    n.OBS <- set.info.temp[, .N]
    n.bins <- min(50, n.OBS) #bin p values; number of tested gene sets up to maximum 50 bins
    
    # force pi0 == 1 if any of the following conditions is not observed:
    #1. the total number of evaluated FDR gene sets exceeds 200
    #2. >= 50% of total genes sets are evaluated OR >= 50 genes sets are evaluated
    if(!(n.EXP > 200 & (n.OBS/n.sets >= 0.5 | n.OBS >= 50)) & est.pi0){ 
      est.pi0 <- FALSE
      e.mss <- paste0("Too few FDR realisations to provide reasonable estimate",
                      " pi0 using histogram method\n", "Setting pi0 to 1 ",
                      "(equivalent to Holm method)")
      warning(e.mss, immediate.=T)
    }
  
    pi0.full <- NA; pi0.pruned <- NA
    if(precise.p.method %in% c("full", "both") & n.perm.p > n.perm.pruning){
      p.vals.full <- data.table::data.table(p.vals, P=p.vals$P.full.null)
      if(est.pi0){
        pi0.full <- pi0.estimator(pp=p.vals.full, n.OBS=n.OBS, n.EXP=n.EXP, 
                                  n.bins=n.bins)
      }else{
        pi0.full <- 1
      }
      q1 <- q.estimator(pp=p.vals.full, pi0=pi0.full, n.OBS=n.OBS, n.EXP=n.EXP)
      set.info.temp <- data.table::merge.data.table(set.info.temp, 
                                                    q1$q[, .(RANK, Q.full.null=Q)], 
                                                    by="RANK")[order(P.full.null)]
    }else{
      set.info.temp[, Q.full.null:=NA]
    }
    if(precise.p.method %in% c("pruned", "both")){
      p.vals.pruned <- data.table::data.table(p.vals, P=p.vals$P.pruned.null)
      if(est.pi0){
        pi0.pruned <- pi0.estimator(pp=p.vals.pruned, n.OBS=n.OBS, 
                                    n.EXP=n.EXP, n.bins=n.bins)
      }else{
        pi0.pruned <- 1
      }
      q2 <- q.estimator(pp=p.vals.pruned, pi0=pi0.pruned, n.OBS=n.OBS, 
                        n.EXP=n.EXP)
      set.info.temp <- data.table::merge.data.table(set.info.temp, 
                                                    q2$q[, .(RANK, Q.pruned.null=Q)], 
                                                    by="RANK")[order(P.pruned.null)]
    }else{
      set.info.temp[, Q.pruned.null:=NA]
    }
    
    if(!is.na(prec)){ # revert set scores
      set.info.temp[, setScore:=setScore/10^prec]
      set.info.temp[, setScore.pruned:=setScore.pruned/10^prec]
    }
    
    #return runs used to estimate FDR
    null.out <- p.vals[REP!=1]
    null.out[, REP:=REP-1] # rescale replicates from 1 to n.FDR
    null.out[, setID:=setID.orig] #revert to original setIDs
    null.out[, setID.orig:=NULL]
    
    #revert IDs back to original
    OUT <- list(set.info=set.info.temp[, .(setName, setID=setID.orig, 
                                           N, N.pruned, setScore, 
                                           setScore.pruned, RANK, P.raw, 
                                           P.raw.pruned, P.full.null, 
                                           P.pruned.null, Q.full.null, 
                                           Q.pruned.null)], 
                set.obj=set.obj.temp[, .(setID=setID.orig, objID=objID.orig)], 
                pi0.est=data.table::data.table(Test=c("Full", "Pruned"),
                                               pi0=c(pi0.full, pi0.pruned)), 
                null.set=null.out, permutation.mat=perm.mat, 
                permutation.type=c("random", "linked")[as.numeric(linked)+1])
    
    cat("Completed q-value calculation step\n")
    cat(paste0(" *** Time elapsed: ", get.time(startT), " *** \n\n"))
  }
  
  #return output
  cat("Done.\n")
  return(OUT)
}
  

#===========================================================================
# Function: polylinkr.wrapper(in.path="./", population=NA, minsetsize=10, 
#                             maxsetsize=1000, merge.set.prop=1, set.info, 
#                             obj.info, set.obj, n.cores="default", linked=T, 
#                             n.perm.p=100000, perm.block.size=1000, prune=T, 
#                             n.perm.pruning=5000, n.FDR=100, n.pruned.sets=50, 
#                             precise.p.method="both", est.pi0=TRUE, 
#                             perm.mat=NULL)
# wrapper function for data input and gene set enrichment testing
#===========================================================================

polylinkr.wrapper <- function(in.path="./", population=NA, minsetsize=10, 
                              maxsetsize=1000, merge.set.prop=1, obj.in.set=T,
                              set.info=NULL, obj.info=NULL, set.obj=NULL, 
                              n.cores="default", linked=T, n.perm.p=100000, 
                              perm.block.size=1000, prune=T, 
                              n.perm.pruning=5000, n.FDR=100, 
                              n.pruned.sets=50, precise.p.method="both", 
                              est.pi0=TRUE, perm.mat=NULL, save=F,
                              out.path="./"){
  
  #call in libraries
  library(data.table)
  library(foreach)
  library(matrixStats)
  library(parallel)
  library(doFuture)
  library(qvalue)
  library(Matrix)
  library(progressr)
  
  # records cumulative run time at specific points
  get.time <- function(t) { # time conversion
    tt <- as.numeric(difftime(Sys.time(), t, units = "secs"))
    H <- tt %/% 3600
    rr <- tt %% 3600
    if(rr > 0) {
      M <- rr %/% 60
      rr2 <- rr %% 60
      if(rr2 > 0) {
        S <- round(rr2)
      }
    }
    return(paste0(H, "h ", M, "m ", S, "s"))
  }
  
  startT <- Sys.time() #initiate timer
  cat("[Part 0] Reading Data...\n")
  
  dd <- read.setobj.tables(in.path=in.path, set.info=set.info, 
                           obj.info=obj.info, set.obj=set.obj,
                           population=population, minsetsize=minsetsize, 
                           maxsetsize=maxsetsize, merge.set.prop=0.95, 
                           obj.in.set=obj.in.set, n.cores=n.cores)
  
  cat("Completed data entry\n")
  cat(paste0(" *** Time elapsed: ", get.time(startT), " *** \n\n"))

  #pruning, with linkage, fast
  pl <- polylinkr(set.info=dd$set.info, obj.info=dd$obj.info, 
                  set.obj=dd$set.obj, n.cores, linked, minsetsize=dd$minsetsize, 
                  n.perm.p, perm.block.size, prune, n.perm.pruning, n.FDR, 
                  n.pruned.sets, precise.p.method, est.pi0, perm.mat, 
                  startT=startT)
  
  #add merged set info to output
  pl$merged.set.info <- dd$set.info.merged
  
  # save to file
  if (save) save(pl, file=file.path(out.path,"pl.R"))
  return(pl)
  
}
  

#===========================================================================
# Function: cluster.genes(set.obj, set.info, use.recomb.rate=F, recomb.rate.path=NA, 
#                         fictive.gene.size=1, n.cores="max")
# create fictive genes for each gene set and generate clustering metrics
#
# - use.recomb.rate: logical, if TRUE genes will be binned based on genetic distance (requires genetic distance files), if FALSE genes will be binned based on physical distance; default = FALSE
# - recomb.rate.path: path to files with recombination rate estimates; users can specify the path there own file with the following fields
#         chr=chromosome name, pos=physical position, cM=inferred genetic position relative to start of chromosome. Alternatively, the if the option 'Bherer2017' is chosen, then data with sex-averaged human rates from https://github.com/cbherer/Bherer_etal_SexualDimorphismRecombination is loaded
# - fictive.gene.size: size of genomic regions used to create fictive genes, if use.recomb.rate = TRUE units are cM, otherwise units are Mb; default = 1 [cM, Mb]
# - n.cores: character or integer specifying the number of cores to use; options include 'max' = all cores, 'default' = max cores - 1, or integer value; default = 'default'
# 
#===========================================================================

cluster.genes <- function(set.obj, set.info, obj.info, use.recomb.rate=F, 
                          recomb.rate.path=NA, fictive.gene.size=1, 
                          n.cores="max"){
  
  #------------------------------------------------------------------#
  ##set up for parallize processing
  #------------------------------------------------------------------#
  
  library(doFuture)
  N.cores <- future::availableCores() 
  
  if(n.cores %in% c("max", "MAX", "Max", "default", "Default", "DEFAULT") | 
     is.numeric(n.cores)){
    if(n.cores %in% c("max", "MAX", "Max")){
      n.cores <- N.cores
      cat(paste0("Using all available (", N.cores, ") cores\n"))
    }
    if(n.cores %in% c("default", "Default", "DEFAULT")){
      n.cores <- N.cores-1
      cat(paste("Using", N.cores-1, "cores\n"))
    }
    if(is.numeric(n.cores)){
      n.cores <- as.integer(n.cores)
      if(n.cores<0 | n.cores>N.cores){
        n.cores <- N.cores-1
        e.mss <- paste("No. cores misspecified: using", n.cores, "cores")
        warning(e.mss, immediate.=T)
      }
    }
  }else{
    n.cores <- N.cores-1
    e.mss <- paste("No. cores misspecified: using", n.cores, "cores")
    warning(e.mss, immediate.=T)
  }
  
  doFuture::registerDoFuture()
  future::plan(multiprocess, workers=n.cores) #uses multicore if available
  
  if(use.recomb.rate){
    
    #------------------------------------------------------------------#
    ##create genetic map for each gene by interpolating gene distance
    #------------------------------------------------------------------#
    
    cat("Creating genetic map of gene positions...\n")
    
    if(recomb.rate.path=="Bherer2017"){
      load("~/Dropbox/R_projects/PolyLinkR/PolyLink_extension/human_sexavg_RR_Bherer2017.rda") ##NOTE THAT THIS PATH WILL NEED TO BE UPDATED FOR PACKAGE
    }else{
      rr.dt <- data.table::fread(recomb.rate.path)
    }
    unq.chr <- sort(rr.dt[, unique(chr)])
    gene.map <- foreach::foreach(chr.now=unq.chr, .combine=rbind) %do% { 
      #use genetic distances
      map.now <- rr.dt[chr==chr.now]
      genes.now <- obj.info[chr==chr.now][order(startpos)]
      genes.now[, midpos:=(startpos+endpos)/2]
      
      rr <- approx(map.now$pos, map.now$cM, xout=genes.now[, midpos])$y 
      
      #correct for non-interpolated regions outside of estimated RR regions
      ends1 <- which(is.na(rr))
      ends2 <- which(!is.na(rr))
      if(length(ends1)>0){
        if(1 %in% ends1){
          rr[1:(min(ends2)-1)] <- 0
        }
        if(nrow(genes.now) %in% ends1){
          rr[max(ends2):nrow(genes.now)] <- max(map.now$cM)
        }
      }
      
      data.table::data.table(genes.now, MapPos=rr)
    }
  }else{ # use physical distances
    gene.map <- obj.info[order(chr, startpos)]
  }
  
  #------------------------------------------------------------------#
  ##create fictive genes for each pathway
  #------------------------------------------------------------------#
  
  if(use.recomb.rate){
    cat("Creating fictive genes using genetic map...\n")
    cat(paste("For each gene set, fictive genes will be created for successive", 
              fictive.gene.size, "cM intervals\n"))
  }else{
    cat("Creating fictive genes using physical map...\n")
    cat(paste("For each gene set, fictive genes will be created for successive", 
              fictive.gene.size, "Mb intervals\n"))
    fictive.gene.size <- fictive.gene.size * 10^6 # convert to megabases
  }
  
  paths <- set.info[, unique(setID)]
  progressr::with_progress({ # allow updating during processing
    prog <- progressr::progressor(steps=length(paths))
    #loop over all paths
    path.fict <- foreach::foreach(p=paths, .combine=rbind,
                          .export=c("paths", "gene.map", "set.obj")) %dopar% { 
      prog(sprintf("STEP=%i", p))
      set.obj.now <- set.obj[setID==p]
      dd.now <- gene.map[objID %in% set.obj.now$objID]
      pp <- dd.now[, .N, by=chr]
      
      if(any(pp$N>1)){ # create fictive genes for each chromosome
        CC <- 0
        M <- foreach(chr.now=pp[, chr], .combine=rbind) %do% {
          if(use.recomb.rate){
            x <- dd.now[chr==chr.now, MapPos]
          }else{
            x <- dd.now[chr==chr.now, (endpos+startpos)/2]
          }
          xx <- as.matrix(dist(x)) < fictive.gene.size # get pairwise distances
          #xx[upper.tri(xx)] <- NA
          x1 <- ceiling(which(xx) / ncol(xx))
          x2 <- which(xx) %% nrow(xx); x2[x2==0] <- nrow(xx)
          XX <- data.table::data.table(X1=x1, I=x2)
          min.X1 <- XX[, min(X1), by=I]
          DD <- data.table(setID=p, objID=dd.now[chr==chr.now, objID], 
                     objID.fictive=min.X1$V1+CC)
          CC <- max(DD$objID.fictive)
          DD
        }
      }else{ # all genes on separate chromosomes
        M <- data.table(setID=p, objID=dd.now[, objID], 
                        objID.fictive=dd.now[, objID])
      }
      M
    }
  })
  
  #return results
  path.fict[, objID.fictive:=paste0(setID, ".", objID.fictive)]
  set.obj.out <- data.table::merge.data.table(path.fict, set.obj, 
                                              by=c("setID", "objID"))
  set.info.out <- data.table::merge.data.table(set.obj.out[, .(N.fictive=uniqueN(objID.fictive)), 
                                                           by=setID],
                                               set.info, by="setID")
  set.info.out[, setClustering:=(N-N.fictive)/N]
  
  return(list(set.info=set.info.out, set.obj=set.obj.out))
}


#===========================================================================
# Function: compare.clustering(set.obj, set.info, use.recomb.rate=F, recomb.rate.path=NA, fictive.gene.size=1, 
#                              n.cores="max", n.perm.p=100000, perm.block.size=1000, bin.size='default',
#                              make.plot=T, plot.path="./")
#
# test if there is any impact of clustering amongst gene sets
# - make.plot: logical, if TRUE a plot comparing the p-values of the two methods is saved, no plotting is performed if FALSE is chosen
# - n.bins: character or integer specifying the number of gene sets to analyse according to their clustring metric; 'default' = or user specified; default = 'default'
# - plot.path: character indicating path to save plot; default = current working directory
# - plot.label: character indicating the the figure file name; default = NA
# - q.sig: numeric value determining the q-value significance threshold used in plots; default = 0.05
#===========================================================================

compare.clustering <- function(set.obj, set.info, obj.info, use.recomb.rate=F, 
                               recomb.rate.path=NA, fictive.gene.size=1, 
                               n.cores="max", n.perm.p=100000, 
                               perm.block.size=1000,  make.plot=T, 
                               n.bins='default', plot.path="./", 
                               plot.label=NA, q.sig=0.05){
  
  # records cumulative run time at specific points
  get.time <- function(t) { # time conversion
    tt <- as.numeric(difftime(Sys.time(), t, units = "secs"))
    H <- tt %/% 3600
    rr <- tt %% 3600
    if(rr > 0) {
      M <- rr %/% 60
      rr2 <- rr %% 60
      if(rr2 > 0) {
        S <- round(rr2)
      }
    }
    return(paste0(H, "h ", M, "m ", S, "s"))
  }
  
  startT <- Sys.time() #initiate timer
  cat("Generating clustering metics...\n")
  
  #cluster genes
  cl.g <- cluster.genes(set.obj, set.info, obj.info,
                        use.recomb.rate=use.recomb.rate, 
                        recomb.rate.path=recomb.rate.path, 
                        fictive.gene.size=fictive.gene.size, n.cores=n.cores)
  
  cat("Completed computing clustering metric for", nrow(set.info), 
      "gene sets\n")
  cat(paste0(" *** Time elapsed: ", get.time(startT), " *** \n\n"))
  
  #run polylinkr
  cat("Performing permutation tests using random selected genes...\n")
  pl.random <- polylinkr(set.info, obj.info, set.obj, n.cores="default", 
                         linked=F, n.perm.p=100000, perm.block.size=1000, 
                         prune=F, startT=startT)
  
  cat("Performing permutation tests using linked sets of genes...\n")
  pl.linked <- polylinkr(set.info, obj.info, set.obj, n.cores="default", 
                         linked=T, n.perm.p=100000, perm.block.size=1000, 
                         prune=F, startT=startT) 
  
  ##TESTS
  if(make.plot){
    cat("Creating plots...\n")
    
    #create successive bins of 100 sets 
    n.sets <- nrow(set.info)
    if(n.bins=="default" | n.bins > n.sets){
      if(n.sets>=500){
        n.bins <- 8
      }else{
        if(n.sets<200){
          if(n.sets<40){
            n.bins <- 2
          }else{
            n.bins <- 4
          }
        }else{
          n.bins <- 6
        }
      }
    }
    bin.size <- ceiling(n.sets/n.bins)
    bin.quants <- (0:n.bins)/n.bins
    cat("Using", n.bins, "bins containing up to", bin.size, "gene sets each\n")
    
    #merge datasets
    pl <- data.table::merge.data.table(pl.random$set.info[, .(setName, 
                                                              P.random=P, 
                                                              Q.random=Q)], 
                                       pl.linked$set.info[, .(setName, 
                                                              P.linked=P, 
                                                              Q.linked=Q)], 
                                       by="setName")
    pl.cl <- data.table::merge.data.table(pl, cl.g$set.info, by="setName")
    
    pl.cl[, clust.bin:=cut(setClustering, quantile(setClustering, bin.quants), 
                           include.lowest=T)]
    pl.cl[, `Q < 0.05`:="Random not sig.\nLinked not sig."]
    pl.cl[Q.random<0.05 & Q.linked>=0.05,
          `Q < 0.05`:="Random sig.\nLinked not sig."]
    pl.cl[Q.random>=0.05 & Q.linked<0.05,
          `Q < 0.05`:="Random not sig.\nLinked sig."]
    pl.cl[Q.random<0.05 & Q.linked<0.05,
          `Q < 0.05`:="Random sig.\nLinked sig."]
    pl.cl[, `Q < 0.05`:=factor(`Q < 0.05`,
                               levels=c("Random not sig.\nLinked not sig.",
                                        "Random sig.\nLinked not sig.",
                                        "Random not sig.\nLinked sig.",
                                        "Random sig.\nLinked sig."))]
    min.x.y <- 1/n.perm.p
    
    library(ggplot2)
    
    plot.label <- ifelse(is.na(plot.label), "linked_vs_random", 
                         paste0(plot.label, "_linked_vs_random"))
    
    #scatter plots
    ggplot(pl.cl, aes(x=P.random, y=P.linked)) +
      geom_point(aes(col=`Q < 0.05`), alpha=0.8) +
      geom_abline(size=0.5) + 
      geom_smooth(size=0.5) +
      scale_y_log10(limits=c(min.x.y, 1)) +
      scale_x_log10(limits=c(min.x.y, 1)) +
      facet_wrap(~clust.bin, nrow=2) + 
      scale_color_manual(values=c("grey50", "orange", "pink", "red"),
                         drop=FALSE) +
      theme_bw()
    
    ggsave(file.path(plot.path, paste0(plot.label, "_scatter_plot.pdf")), 
           height=max(4, 0.75*n.bins), width=max(6, 1.25*n.bins), dpi=600)
    
    #distributions
    pl.cl.comb <- rbind(data.table(Clustering=pl.cl$clust.bin, Perm="Random", 
                                   P=pl.cl$P.random, Q=pl.cl$Q.random),
                        data.table(Clustering=pl.cl$clust.bin, Perm="Linked", 
                                   P=pl.cl$P.linked, Q=pl.cl$Q.linked))
    ggplot(pl.cl.comb, aes(x=P, col=Perm)) +
      geom_density() +
      scale_x_log10() +
      facet_wrap(~Clustering, nrow=2) + 
      scale_color_manual(values=c("blue", "red")) +
      theme_bw()
    
    ggsave(file.path(plot.path, paste0(plot.label, "_density_plot.pdf")), 
           height=max(4, 0.75*n.bins), width=max(6, 1.25*n.bins), dpi=600)
  }

  #output
  return(pl.cl[, c(1, 6:9, 10, 2:5, 11), with=F])
}


##END OF FUNCTION CODE
#===========================================================================
#===========================================================================

#JD: don't forget to remove code below


##-----------------------------------------##
#TESTING
##-----------------------------------------##

library(data.table)
library(foreach)
library(matrixStats)
library(parallel)
library(doFuture)
library(qvalue)
library(Matrix)
library(progressr)

# dd <- read.setobj.tables(in.path="~/Dropbox/R_projects/PolyLink_code/Example/PolyLink/Input", 
#                        population="Anatolia_EF", minsetsize=10, 
#                        maxsetsize=1000, merge.set.prop=0.95, 
#                        n.cores=4)
# 
# #test clustering
# cl1 <- compare.clustering(set.obj=dd$set.obj, set.info=dd$set.info, obj.info=dd$obj.info,
#                           use.recomb.rate=F, recomb.rate.path=NA, fictive.gene.size=1, 
#                           n.cores="max", n.perm.p=100000, perm.block.size=1000, make.plot=T, 
#                           n.bins='default', plot.path="~/Dropbox/R_projects/PolyLinkR/PolyLink_extension/plots", 
#                           q.sig=0.05, plot.label=paste0(dd$population, "_PhysicalMap"))
# 
# cl2 <- compare.clustering(set.obj=dd$set.obj, set.info=dd$set.info, obj.info=dd$obj.info, 
#                           use.recomb.rate=T, recomb.rate.path="~/Dropbox/R_projects/PolyLinkR/PolyLink_extension/RecombRate", 
#                           fictive.gene.size=1, n.cores="max", n.perm.p=100000, perm.block.size=1000, 
#                           make.plot=T, n.bins='default', plot.path="~/Dropbox/R_projects/PolyLinkR/PolyLink_extension/plots", 
#                           q.sig=0.05, plot.label=paste0(dd$population, "_GeneticMap"))
# 
# #pruning, with linkage, fast
# pl1 <- polylinkr(set.info=dd$set.info, obj.info=dd$obj.info, set.obj=dd$set.obj, 
#                  n.cores="max", linked=T, minsetsize=dd$minsetsize, n.perm.p=100000, 
#                  perm.block.size=1000, prune=T, n.perm.pruning=5000, n.FDR=100, 
#                  n.pruned.sets=50, precise.p.method="both", est.pi0=TRUE, perm.mat=NULL, 
#                  startT=NA)
# 
# save(pl1, file="~/Dropbox/R_projects/PolyLinkR/PolyLink_extension/polylink.test.linked.Robj")
# 
# #pruning, without linkage, fast
# pl2 <- polylinkr(set.info=dd$set.info, obj.info=dd$obj.info, set.obj=dd$set.obj, 
#                  n.cores=4, linked=F, minsetsize=dd$minsetsize, n.perm.p=100000, 
#                  perm.block.size=1000, prune=T, n.perm.pruning=5000, n.FDR=100, 
#                  n.pruned.sets=50, precise.p.method="both", est.pi0=TRUE, perm.mat=NULL,
#                  startT=NA)
# 
# save(pl2, file="~/Dropbox/R_projects/PolyLinkR/PolyLink_extension/polylink.test.unlinked.Robj")
# 
# #pruning, with linkage, fast
# pl3 <- polylinkr(set.info=dd$set.info, obj.info=dd$obj.info, set.obj=dd$set.obj, 
#                  n.cores="max", linked=T, minsetsize=dd$minsetsize, n.perm.p=100000, 
#                  perm.block.size=1000, prune=F, n.perm.pruning=5000, n.FDR=100, 
#                  n.pruned.sets=50, precise.p.method="both", est.pi0=TRUE, perm.mat=NULL, 
#                  startT=NA)
# 
# save(pl3, file="~/Dropbox/R_projects/PolyLinkR/PolyLink_extension/polylink.test.linked.NoPruning.Robj")
# 
# #pruning, without linkage, fast
# pl4 <- polylinkr(set.info=dd$set.info, obj.info=dd$obj.info, set.obj=dd$set.obj, 
#                  n.cores=4, linked=F, minsetsize=dd$minsetsize, n.perm.p=100000, 
#                  perm.block.size=1000, prune=F, n.perm.pruning=5000, n.FDR=100, 
#                  n.pruned.sets=50, precise.p.method="both", est.pi0=TRUE, perm.mat=NULL,
#                  startT=NA)
# 
# save(pl4, file="~/Dropbox/R_projects/PolyLinkR/PolyLink_extension/polylink.test.unlinked.NoPruning.Robj")


