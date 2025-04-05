### 1 lod interval, hk method, sim.geno, ndraw=500, step=1.
setwd("/Users/limengxie/Desktop/JASHS_Revision/final_mqm_181/")
#install.packages("qtl")
#########################################################################################
#url <- "http://www.rqtl.org/download/old/qtl_1.38-4.tar.gz"
#install.packages(url,lib="/Users/limengxie/Library/R/3.3/library",type="source",repo=NULL)
#install.packages("qtl")
library(qtl)
library(rlecuyer)
library(snow)
#install.packages("qtl")
cross <- read.cross ("csv", "./","Mapping181.csv",estimate.map=F)
cross<-jittermap(cross)
## estimate genotyping error first
#loglik <- err <- c(0.005, 0.01,0.015, 0.02,0.025,0.03,0.035,0.04, 0.045,0.05) 
#for(i in seq(along=err)) {
 # cat(i, "of", length(err), "\n")
 # tempmap <- est.map(cross, error.prob=err[i])
  #loglik[i] <- sum(sapply(tempmap, attr, "loglik")) }
#lod <- (loglik - max(loglik))/log(10)
#plot(err, lod, xlab="Genotyping error rate", xlim=c(0.0, 0.05),
 #    ylab=expression(paste(log[10], " likelihood")))
###
cross<-calc.genoprob(cross, map.function="kosambi",error.prob = 0.025)
## normality test
norm<-{}
norm2<-{}
for(i in 1:(nphe(cross)-1)){
  x<-cross$pheno[i]
  x<-x[,1]
  y<-shapiro.test(x)
  ShapiroWilk.pvalue<-y$p.value
  Phenotype<-colnames(cross$pheno[i])
  norm<-cbind(Phenotype, ShapiroWilk.pvalue)
  norm2<-rbind(norm2, norm)
}
write.csv(file="Normality test p-values.csv", norm2)
## Batch effect
pdf(file="Batch effects.pdf", width=11, height=8.5)
for(i in 1:(nphe(cross)-1)){
  par(mfrow=c(1,2), las=1, cex=0.8)
  means<-apply(cross$pheno[i], 1, mean)
  plot(means)
  plot(sample(means), xlab="Random Index", ylab="means", main=colnames(cross$pheno[i]))
}
dev.off()
############

##Genotype graph
pdf(file="Graphical Genotype.pdf", width=11, height=8.5)
geno.image(cross)
dev.off()

##########################################################################################
pdf(file="Estimated Genetic Map.pdf", width=11, height=8.5)
par(cex=0.6)
plot.map(cross, show.marker.names=FALSE)
dev.off()


## permutation first
cross.sc1<-scanone(cross, pheno.col=1:(nphe(cross)-1), method="hk", use="all.obs", model="normal")
cross.sc1.perms<-scanone(cross, pheno.col=1:(nphe(cross)-1), method="hk", n.perm=1000, verbose=TRUE, model="normal", n.cluster=12)

##########################################################################################
sum<-summary(cross.sc1, threshold=0, perms=cross.sc1.perms, pvalues=TRUE, format="tabByCol", ci.function="lodint", drop=1, expandtomarkers=TRUE)
space<-" "
for(i in 1:(nphe(cross)-1)){
  x<-capture.output(sum[[i]])
  cat(colnames(cross$pheno[i]), file="10_Initial QTL hits by phenotype.txt", sep="\n", append=TRUE)
  cat(x, file="10_Summary of top hits by phenotype.txt", sep="\n", append=TRUE)
  cat(space, file="10_Summary of top hits by phenotype.txt", sep="\n", append=TRUE)
}
##########################################################################################

#(19) Export csv file with LOD scores at every marker. 

##########################################################################################
lods<-capture.output(cross.sc1)
cat(lods, file="11_LOD scores for every marker.txt", sep="\n")
##########################################################################################

#(20) Produce PDF file indicating genome scans with alpha=.05 threshold marked

##########################################################################################
z<-summary(cross.sc1.perms, alpha=.05)
pdf(file="QTL Plots.pdf", width=11, height=8.5)
for(i in 1:(nphe(cross)-1)){
  plot(cross.sc1, lodcolumn=i, lwd=1.5, gap=0, bandcol="gray70", incl.markers=TRUE, main=colnames(cross$pheno[i]), xlab=c("Threshold for alpha=.05 using 1000 permutations", z[i]))
  add.threshold(cross.sc1, perms=cross.sc1.perms, alpha=0.05, lodcolumn=i, gap=0)
}
dev.off()
##########################################################################################

#(21) Produce a csv file with the map positions for each marker.

#*NOTE*# This may produce a warning consistent with the number of chromosomes. 
#This is expected as the program is making the user aware that it is appending the column names to the file. 

##########################################################################################
newmap<-pull.map(cross)
for(i in 1:length(names(newmap))){
  snps<-names(newmap[[i]])
  gm<-c(snps, newmap[[i]])
  gm2<-matrix(gm, ncol=2)
  write.table(file="13_Genetic Map Positions.csv", sep=",", append=TRUE, gm2)
}
##########################################################################################

#(22)run sim.geno (necessary for some of the following functions)

##########################################################################################
cross2<-sim.geno(cross, n.draws=500, error.prob=0.025, step=1, map.function="kosambi")
##########################################################################################

######################################################################################
########																	##########
######## Steps #(23) - #(32) represent a loop that will scan for additional ##########
######## qtl, run model selection, and fit the qtl model for each phenotype ##########
######## in turn automatically. Run the entire code (Steps #(23) - #(32) )  ##########
######## simultaneously.                                                    ##########
########																	##########
######################################################################################


######################################
##Model Selection and QTL refinement##
######################################

#(23) Specifies the phenotype you wish to work with. If you are troubleshooting or otherwise want to run each step individually for a given phenotype rather than the whole things a loop, just replace k in ' phenotype<-k ' with the number of the phenotype you want to use, ignore the line ' for(k in 2:nphe(cross)){ '  and run each step individually. 

##########################################################################################
for(k in 1:(nphe(cross)-1)){
  
  phenotype<-k
  
  pheno<-colnames(cross$pheno[phenotype])
  ##########################################################################################
  
  #(24) Makes a qtl object containing the ALL the qtl from file="10_Initial QTL hits by phenotype" 
  #for the phenotype you specified.  
  
  ##########################################################################################
  chromo<-sum[[pheno]]
  chr<-{}
  pos<-{}
  for(i in 1:nrow(chromo)){
    chr1<-chromo[i,1]
    chr2<-as.numeric(as.character(chr1))
    chr<-c(chr, chr2)
  }
  for(i in 1:nrow(chromo)){
    pos1<-chromo[i,2]
    pos<-c(pos, pos1)
  }
  if(is.na(chr[1])==FALSE) {
    qtl<-makeqtl(cross, chr=chr, pos=pos, what="prob")
    createqtl<- paste("Q", 1:qtl$n.qtl, sep="")
    formula<-as.formula(paste("y ~ ", paste(createqtl, collapse= "+")))
    ##########################################################################################
    
    #(25) Scans for additional linked QTL conditioning on the QTL already detected. .
    
    #*NOTE*# You may receive warning messages about dropping individuals with missing phenotype data 
    #and/or that the column names in scanone input do not match those in perms input. These are both 
    #expected under some circumstances and do not effect the output of the code. Also you may recieve 
    #a warning that there is no chromosome number NA if no additional QTL are to be found.
    
    ##########################################################################################
    cross.aq<-addqtl(cross, pheno.col=phenotype, qtl=qtl, formula=formula, method="hk", model="normal")
    sub.perms<-subset(cross.sc1.perms, lodcolumn=phenotype)
    xx<-capture.output(summary(cross.aq, perms=sub.perms, alpha=.05, pvalues=TRUE, format="tabByCol", ci.function="lodint", drop=1, expandtomarkers=TRUE))
    xxy<-c(pheno,xx)
    cat(xxy, file="14_Additional QTL hits by phenotype.txt", sep="\n", append=TRUE)
    cat(space, file="14_Additional QTL hits by phenotype.txt", sep="\n", append=TRUE)
    ##########################################################################################
    
    #(26) Adds additional QTL (if any) in the file "14_Additional QTL hits by phenotype.txt" 
    #to the qtl object and update the formula object.
    
    ##########################################################################################
    sum.aq<-summary(cross.aq, perms=sub.perms, alpha=.05, pvalues=TRUE, format="tabByCol", ci.function="lodint", drop=0.95, expandtomarkers=TRUE)
    chr.aq<-{}
    pos.aq<-{}
    for(i in 1:nrow(sum.aq$lod)){
      chr.aq1<-sum.aq$lod[i,1]
      chr.aq2<-as.numeric(as.character(chr.aq1))
      chr.aq<-c(chr.aq, chr.aq1)
    }
    if(is.na(chr.aq)==TRUE) { chr.aq<-chr } else { chr.aq<-c(chr, chr.aq) }
    
    for(i in 1:nrow(sum.aq$lod)){
      pos1.aq<-sum.aq$lod[i,2]
      pos.aq<-c(pos.aq, pos1.aq)
    }
    if(is.na(pos.aq)==TRUE) { pos.aq<-pos } else { pos.aq<-c(pos,pos.aq) }
    
    qtl<-makeqtl(cross, chr=chr.aq, pos=pos.aq, what="prob")
    createqtl<- paste("Q", 1:qtl$n.qtl, sep="")
    formula<-as.formula(paste("y ~ ", paste(createqtl, collapse= "+")))
    ##########################################################################################
    
    #(27) Uses forward selection and backward elimination model selection to probe the model space 
    #for the best fit QTL model explaining your data.
    
    ##########################################################################################
    pen<-summary(sub.perms)
    cross.sw<-stepwiseqtl(cross, pheno.col=phenotype, qtl=qtl, formula=formula, method="hk", penalties=pen, model="normal", additive.only=TRUE)
    swQTL<-capture.output(print(cross.sw))
    swQTL2<-c(pheno, swQTL)
    cat(swQTL2, file="15_Additional QTL hits from stepwise analysis.txt", sep="\n", append=TRUE)
    cat(space, file="15_Additional QTL hits from stepwise analysis.txt", sep="\n", append=TRUE)
    ##########################################################################################
    
    #(28) Adds additional QTL (if any) to the qtl object found by stepwise model selection.
    
    ##########################################################################################
    sum.sw<-summary(cross.sw)
    pos.sw<-{}
    for(i in 1:nrow(sum.sw)){
      chr.sw1<-sum.sw[i,2]
      chr.sw2<-as.numeric(as.character(chr.sw1))
      chr.sw<-c(chr.sw, chr.sw2)
    }
    for(i in 1:nrow(sum.sw)){
      pos1.sw<-sum.sw[i,3]
      pos.sw<-c(pos.sw, pos1.sw)
    }
    
    qtl2<-makeqtl(cross, chr=chr.sw, pos=pos.sw, what="prob")
    createqtl<- paste("Q", 1:qtl2$n.qtl, sep="")
    formula<-as.formula(paste("y ~ ", paste(createqtl, collapse= "+")))
    rqtl<-refineqtl(cross, pheno.col=phenotype, qtl=qtl2, method="hk", model="normal")
    ##########################################################################################
    
    #(29) Writes a file containing the peak marker and the flanking markers representing the 1 LOD interval of each QTL.
    
    ##########################################################################################
    Q<-"Q"
    space<-" "
    for (i in 1:rqtl$n.qtl){
      interval<-capture.output(lodint(rqtl, qtl.index=i, drop=1, expandtomarkers=TRUE))
      q<-paste(Q, i, sep="")
      interval.new<-c(pheno, q, interval, space)
      cat(interval.new, file="16_Summary of Final QTL Intervals.txt", sep="\n", append=TRUE)
    }
    ##########################################################################################
    
    #(30) Writes a csv file containing the results of ANOVA for the full and reduced models, 
    #the % variance explained by each QTL and the estimated effect size 
    #(half the distance between the means for each genotype class in the case of RILs).
    
    #*NOTE*# You may receive a warning here about dropping individuals with missing phenotypes. 
    #This is expected if such a case exists and does not effect the output of the code.
    
    ##########################################################################################
    cross.ests<-fitqtl(cross, pheno.col=phenotype, qtl=rqtl, formula=formula, method="hk", dropone=TRUE, get.ests=TRUE, model="normal")
    ests<-capture.output(summary(cross.ests))
    ests<-c(pheno, ests)
    write(ests, file="17_ANOVA results and QTL effect estimates.txt", sep="\n", append=TRUE)
    ##########################################################################################
    
    #(31) Generates a pdf file of the effect plots and a text file with means and standard error for each genotype class. 
    #Generates a separate pdf file for every QTL.
    
    #*NOTE*# This will likely produce a warning indicating that column names are being appended to file.
    
    ##########################################################################################
    for(i in 1:length(chr.sw)) {
      b<-paste("Q",i, ".pdf",sep="")
      file<-paste("18_Marker effect plots", pheno, b)
      mar<-find.marker(cross2, chr=chr.sw[i], pos=pos.sw[i])
      pdf(file=file, width=11, height=8.5)
      plot.pxg(cross, marker=mar, pheno.col=phenotype)
      dev.off()
      phenoqtl<-paste(pheno, b)
      pheno.eff<-c(phenoqtl, space)
      means<-effectplot(cross2, pheno.col=phenotype, mname1=mar, draw=FALSE)
      cat(pheno.eff, file="19_means and SE.txt", sep="\n", append=TRUE)
      write.table(means$Means, file="19_means and SE.txt", sep=",", col.names="Means", row.names=TRUE, append=TRUE)
      cat(space, file="19_means and SE.txt", sep="\n", append=TRUE)
      write.table(means$SEs, file="19_means and SE.txt", sep=",", col.names="Standard Error", row.names=TRUE, append=TRUE)
      cat(space, file="19_means and SE.txt", sep="\n", append=TRUE)
    }
    ##########################################################################################
    
    #(32) Contingency statement if no QTL are detected in the initial genome scan for a given phenotype. You can ignore this step if you are not running the entire model selection section as a loop. 
    
    ##########################################################################################
  } else { null<-"	There were no LOD peaks above the threshold"
  cat(pheno, file="14_Additional QTL hits by phenotype.txt", sep="\n", append=TRUE)
  cat(null, file="14_Additional QTL hits by phenotype.txt", sep="\n", append=TRUE)
  cat(space, file="14_Additional QTL hits by phenotype.txt", sep="\n", append=TRUE)
  cat(pheno, file="15_Additional QTL hits from stepwise analysis.txt", sep="\n", append=TRUE)
  cat(null, file="15_Additional QTL hits from stepwise analysis.txt", sep="\n", append=TRUE)
  cat(space, file="15_Additional QTL hits from stepwise analysis.txt", sep="\n", append=TRUE)
  cat(pheno, file="16_Summary of Final QTL Intervals.txt", sep="\n", append=TRUE)
  cat(null, file="16_Summary of Final QTL Intervals.txt", sep="\n", append=TRUE)
  cat(space, file="16_Summary of Final QTL Intervals.txt", sep="\n", append=TRUE)
  cat(pheno, file="17_ANOVA results and QTL effect estimates.txt", sep="\n", append=TRUE)
  cat(null, file="17_ANOVA results and QTL effect estimates.txt", sep="\n", append=TRUE)
  cat(space, file="17_ANOVA results and QTL effect estimates.txt", sep="\n", append=TRUE)
  cat(pheno, file="19_means and SE.txt", sep="\n", append=TRUE)
  cat(null, file="19_means and SE.txt", sep="\n", append=TRUE)
  cat(space, file="19_means and SE.txt", sep="\n", append=TRUE) }
}
##########################################################################################

