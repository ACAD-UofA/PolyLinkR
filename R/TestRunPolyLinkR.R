linked <- T
prune <- T
toyexample <- F
txt <- "_min"

in.path="~/Documents/PolyLinkR/data"
out.path="~/Documents/PolyLinkR/data"
population="HumanPops"
minsetsize=10
maxsetsize=1000
merge.set.prop=0.95
n.cores=4

if (toyexample){
  # toy example parameters
  maxsetsize = 50; 
  n.perm.p=500; perm.block.size=50 
  n.perm.pruning=250; n.FDR=25; n.pruned.sets=20
} else {
  maxsetsize = 1000; 
  n.perm.p=100000; perm.block.size=1000 
  n.perm.pruning=5000; n.FDR=100; n.pruned.sets=50
}


# first use polysel to correct for gene length bias
source('~/Documents/Polysel/R/polysel.R')
res<-ReadSetObjTables(in.path=in.path,
                      set.info.file="HumanPops_SetInfo.txt",
                      set.obj.file="HumanPops_SetObj.txt",
                      obj.info.file="HumanPops_ObjInfo.txt",
                      minsetsize=1,
                      obj.in.set=F,
                      merge.similar.sets=F)
obj.info<-res$obj.info
set.info<-res$set.info
set.obj<-res$set.obj

obj.info<-AssignBins(obj.info,fld="SNPcount", min.bin.size = 1000)
obj.info<-obj.info[order(obj.info$objID),]
obj.info<-RescaleBins(obj.info)

# return original IDs
m<-match(set.obj$setID,set.info$setID)
set.obj$setID<-set.info$setID.orig[m]
m<-match(set.obj$objID,obj.info$objID)
set.obj$objID<-obj.info$objID.orig[m]
set.info$setID<-set.info$setID.orig
obj.info$objID<-obj.info$objID.orig
set.info$setID.orig<-NULL
obj.info$objID.orig<-NULL

source('~/Documents/PolyLinkR/R/polylinkr.R')
# Testing polylinkr.wrapper with objects as given

pl<-polylinkr.wrapper(in.path=in.path, population = population,
                      set.info=set.info, 
                      obj.info=obj.info, set.obj=set.obj,
                      obj.in.set=F,
                      linked = linked,
                      n.perm.p=n.perm.p, 
                      prune = prune,
                      perm.block.size=perm.block.size, # default: 1000
                      n.perm.pruning=n.perm.pruning, # default: 5000
                      n.FDR=n.FDR, # default: 100
                      n.pruned.sets=n.pruned.sets, # default: 50
                      minsetsize = minsetsize,
                      maxsetsize = maxsetsize, 
                      merge.set.prop = 0.95, n.cores=n.cores,
                      save=T, out.path=out.path)

write.table(pl$set.info, 
            file=file.path(out.path, 
                           paste0("res_set_info_", 
                                  ifelse(linked,"linked","random"), "_", 
                                  ifelse(prune,"pruned","notpruned"),
                                  ifelse(toyexample,"_toy",""), 
                                  txt, ".txt")),
            quote=F, row.names=F, sep="\t")
