library(stringr)
library(qiimer)
growth_import <-read.csv("GrowthParameters.csv",stringsAsFactors=FALSE)
AT_import <- read.csv("ColdgrowATFreq_CrossSect.csv", stringsAsFactors=FALSE, header=FALSE)
All_ATs <- AT_import[,1]
closest.allele<-c()
mytree<-read.tree("RAxML_bipartitionsBranchLabels_3.9.17.tree")
pw.dist<-cophenetic.phylo(mytree)
good.alleles<-str_pad(growth_import$rpoBAT,4,pad="0")
allele_samp <- All_ATs
for (a in 1:length(allele_samp)){
  newsamp<-str_pad(allele_samp[a],4,pad="0")
  if (newsamp%in%good.alleles){
    splitat.best<-substr(newsamp,regexpr("[^0]",newsamp),nchar(newsamp))
    closest.allele<-append(closest.allele,splitat.best)
  }
  else{
    at.query<-colnames(pw.dist)[which(grepl(paste("AT",newsamp,sep=""),colnames(pw.dist))==TRUE)]
    mindist<-Inf
    for (g in 1:length(good.alleles)){
      gquery<-colnames(pw.dist)[which(grepl(paste("AT",good.alleles[g],sep=""),colnames(pw.dist))==TRUE)]
      getdist<-as.numeric(dist_subset(pw.dist,c(gquery,at.query)))
      print(gquery)
      print(getdist)
      if (length(getdist)>0 && getdist<mindist){
        mindist<-getdist
        minname<-gquery
        ming<-good.alleles[g]
      }
    }
    splitat.best<-substr(ming,regexpr("[^0]",ming),nchar(ming))
    closest.allele<-append(closest.allele,splitat.best)
  }}
closest.allele

AT_import$closestAT <- closest.allele
