# ./Phased_consensus.R REFERENCE FATHER MOTHER BED MASKED_CNS OUTPUT_NAME
library(Biostrings)

args <- commandArgs(trailingOnly = TRUE)

ref<-readDNAStringSet(args[1], format="fasta")
father<-readDNAStringSet(args[2], format="fasta")
mother<-readDNAStringSet(args[3], format="fasta")
bed<-read.table(args[4])
mask<-readDNAStringSet(args[5], format="fasta")
mask2<-mask
startCMP<-0
endCMP<-0
file.remove(paste(args[6],".choose_blocks.txt",sep=""), showWarnings=F)

aln0 <- pairwiseAlignment(ref[[1]], mask[[1]], type="global") # Computes global alignment with Needleman-Wunsch algorithm. 
minus<-start(insertion(aln0))
startsub<-start(subject(aln0))-1
for (i in 1:length(father)){
  start=bed[i,2]+1
  end=bed[i,3]
  startCMP<-start-startsub+sum(width(insertion(aln0))[start>=minus])
  endCMP<-end-startsub+sum(width(insertion(aln0))[end>=minus])
  aln1 <- pairwiseAlignment(subseq(DNAString(as.character(pattern(aln0))),startCMP,endCMP), father[i], type="global") # Computes global alignment with Needleman-Wunsch algorithm. 
  aln2 <- pairwiseAlignment(subseq(DNAString(as.character(pattern(aln0))),startCMP,endCMP), mother[i], type="global") # Computes global alignment with Needleman-Wunsch algorithm. 
  if (score(aln1)>score(aln2)) {
    #    mask<-replaceAt(mask,IRanges((bed[i,2]+1),father[i])) #Other way to do it for newer versions of R
#    subseq(mask,start=start,end=end)<-father[i]
#    subseq(mask2,start=start,end=end)<-mother[i]
    write("father",file=paste(args[6],".choose_blocks.txt",sep=""),append=T)
  }else{
    #    mask<-replaceAt(mask,IRanges((bed[i,2]+1),mother[i])) #Other way to do it for newer versions of R
#    subseq(mask,start=start,end=end)<-mother[i]
#    subseq(mask2,start=start,end=end)<-father[i]
    write("mother",file=paste(args[6],".choose_blocks.txt",sep=""),append=T)
  }
}

#names(mask)<-args[6]
#names(mask2)<-paste(args[6],"_NOT",sep="")
#writeXStringSet(mask,paste(args[6],".fasta",sep=""))
#writeXStringSet(mask2,paste(args[6],".NOT.fasta",sep=""))
