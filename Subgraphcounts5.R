#Code to get counts
#Code requires an edgelist tab separeted with node indices starting from 1 (NOT ZERO).
#To call this code use
# Rscript Subgraphcounts.R myargs.txt
#myargs.txt contains three lines giving (1) the paths to the edgelist (2) the output file and (3) orbits or noorbits.
#(4) "four" "five" (3 has to be orbits)


##############################
#read_argments from 

args=(commandArgs(TRUE))
batch_args <- read.table(args[1],header = FALSE)


elistdir= as.character(batch_args[1,1])
writeto= as.character(batch_args[2,1])
wantorbits= as.character(batch_args[3,1])
orbitsize= as.character(batch_args[4,1])

if(orbitsize=="five" & wantorbits=="noorbits") print("Only orbits are computed with five nodes.")

##############################
#elistdir="~/Desktop/aring_graph.txt"
#writeto="~/Desktop/example_subgraphs_aring.txt"
require(orca)
#
orca_orbitcounts=function(e_list_num, subgraph_size="four"){
  #If graph does not have edges
  if(length(e_list)<1){ 
    print("empty edgelist?");stop("empty edgelist?")
  }
  #If graph has edges#it now starts with 1, required by orca
  e_list=e_list_num
  if(!is.integer(e_list_num)){
    e_list=matrix(ncol=2,nrow=nrow(e_list_num))
    e_list[,1]=as.integer(e_list_num[,1]); e_list[,2]=as.integer(e_list_num[,2])
  }
  rm(e_list_num)
  if(min(e_list) <1 | min(e_list) > 1 ){ print("Wrong indexing for orca. It has to start at 1 and indices must be integers.");stop("Wrong indexing for orca.  It has to start at 1 and indices must be integers.")}
  if(subgraph_size=="four")
  orbitcounts=count4(e_list) #needs integer matrix
  if(subgraph_size=="five")
    orbitcounts=count5(e_list) #needs integer matrix
  return(orbitcounts)
}
##########/////////
orca_subgraphcounts=function(e_list){
  orbits=orca_orbitcounts(e_list,"four")
  subgraph=list(1,2:3,4,5:6,7:8,9,10:12,13:14,15)
  if(nrow(orbits) <1 ){
    print("empty edgelist?"); stop("empty edgelist?")
  }
  if(nrow(orbits) >= 1){
    counts=sapply(subgraph,function(l){if(length(l)<2) return( orbits[,l] ); rowSums(orbits[,l])  })
   }
  colnames(counts)=c("edge","2s","tri","4l","3s","sq","triedg","sqdiag","4clique")
  return(counts)  
}
###################
e_list=read.table(file = elistdir,sep = "\t")

if(wantorbits=="orbits"){
result= orca_orbitcounts(e_list,orbitsize)
}
if(wantorbits=="noorbits")
  result= orca_subgraphcounts(e_list)

#
#Write results
write.table(result,file=writeto,sep = "\t",row.names = FALSE,col.names = TRUE)


###OTHER INFORMATION##
#ABOUT COMMAND LINE CALL: TAKE THE FOLLOWING CODE AT THE BEGINING OF THE SCRIPT
#cmd <- paste(commandArgs(), collapse=" ")
#cat("How R was invoked:\n");
#cat(cmd, "\n")





