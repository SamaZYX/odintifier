library(Biostrings)

args <- commandArgs(trailingOnly = TRUE)

blocks<-read.table(args[1], header=F)
phaseout<-read.table(args[2], header=F, fill=T)
bamzero<-readDNAStringSet(args[3], format="fasta")
bamone<-readDNAStringSet(args[4], format="fasta")
file.remove(paste(args[5],".new_choose_blocks.txt",sep=""))

n=0
n2=0
seqzero=""
seqbamzero=""
seqbamone=""
flag=0

for (i in 2:nrow(phaseout)){
  if(phaseout[i,5]==""){
    if(flag==1){
      n2=n2+1
      aln1 <- pairwiseAlignment(DNAString(seqbamzero),DNAString(seqzero), type="global")
      aln2 <- pairwiseAlignment(DNAString(seqbamone),DNAString(seqzero), type="global")
      if (score(aln1)>=score(aln2)) {
        if(blocks[n2,1]=="0"){
          write("0",file=paste(args[5],".new_choose_blocks.txt",sep=""),append=T)
        }else{
          write("1",file=paste(args[5],".new_choose_blocks.txt",sep=""),append=T)
        }
      }else{
        if(blocks[n2,1]=="0"){
          write("1",file=paste(args[5],".new_choose_blocks.txt",sep=""),append=T)
        }else{
          write("0",file=paste(args[5],".new_choose_blocks.txt",sep=""),append=T)
        }
      }
      seqzero=""
      seqbamzero=""
      seqbamone=""
      flag=0
    }
  }else{
    n=n+1
    flag=1
    seqzero=paste(seqzero,phaseout[i,5],sep="")
    seqbamzero=paste(seqbamzero,bamzero[n],sep="")
    seqbamone=paste(seqbamone,bamone[n],sep="")
  }
}

if(flag==1){
  n2=n2+1
  aln1 <- pairwiseAlignment(DNAString(seqbamzero),DNAString(seqzero), type="global")
  aln2 <- pairwiseAlignment(DNAString(seqbamone),DNAString(seqzero), type="global")
  if (score(aln1)>=score(aln2)) {
    if(blocks[n2,1]=="0"){
      write("0",file=paste(args[5],".new_choose_blocks.txt",sep=""),append=T)
    }else{
      write("1",file=paste(args[5],".new_choose_blocks.txt",sep=""),append=T)
    }
  }else{
    if(blocks[n2,1]=="0"){
      write("1",file=paste(args[5],".new_choose_blocks.txt",sep=""),append=T)
    }else{
      write("0",file=paste(args[5],".new_choose_blocks.txt",sep=""),append=T)
    }
  }
}

