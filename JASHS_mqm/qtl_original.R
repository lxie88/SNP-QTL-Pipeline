#As coded, this script is designed to implement the scanone() function of Rqtl with the EM algorithm. 
#See A Guide to Mapping with R/QTL by Karl Broman (ISBN: 978-0-387-92124-2) for a comprehensive tutorial and guidelines.
#Also a helpful forum moderated by Dr. Broman (as of the date of this publication) 
#can be found at http://groups.google.com/group/rqtl-disc/

#This script differs from v5 as the model selection and model fitting is fully automated and will run through all your phenotypes at once. It is highly recommneded that you check your phenotype and genotype data before running the autmated model selection. If the script throws an error and does not complete, it's best to investigate step by step using v5. 

####LOAD THE CROSS####
setwd("/users/limengxie/Desktop/JASHS_mqm/")
#(1) Load R/qtl

##########################################################################################
library(qtl)
library(rlecuyer)
library(snow)
##########################################################################################

#(2) Read in the genotype and phenotype data. Be sure the file structure is of class "bc"
#for a RIL or a backcross population. As in AA=A, and BB=H, AB=NA. Otherwise structure the 
#data as an intercross (AA=A, BB=B, and AB=H). Hets are treated as missing data in a selfed 
#RIL population. See ?read.cross() for more information. 

#READING THE DATA IN WITH THIS STATEMENT SIMULTANEOUSLY ESTIMATES THE GENETIC MAP POSITIONS 
#OF EACH OF THE MARKERS IF NOT ALREADY PRESENT IN YOUR DATA FILE. DEPENDING ON THE NUMBER OF 
#INDIVIDUALS AND THE NUMBER OF MARKERS, THIS FUNCTION COULD TAKE A WHILE. IF POSISTION INFORMATION 
#IN YOUR DATA FILE REPRESENTS PHYSICAL DISTANCES YOU MAY CONSIDER USING rescalemap() TO CONVERT 
#PHYSICAL POSITIONS TO CENTIMORGANS. BE SURE TO USE A SCALE APPROPRIATE FOR YOUR ORGANISM. 
#SCALE ARGUMENT SHOULD BE SET TO 1/(#bp	in a cM). SEE ?rescalemap() FOR MORE INFO.

#Note the "csvr" format. If your data is not formatted as a rotated comma separated file,
#then the format argument should be changed to reflect it's format. See ?read.cross() for more info.

##########################################################################################
cross<-read.cross(format ="csv", file="./Mapping181.csv", estimate.map=F)
##########################################################################################

#(3) Change population class from bc to RIL if you are working with a RIL population. 
#If you are working with a bc or an intercross and coded your data appropriately, you may skip this step. 

#R/QTL has strange support for RILs, it assumes you have an F2 if you code the genotypes correctly. 
#Hence why your file structure has to simulate a BC. In a BC there are only two genotype classes (AA and AB). 
#Likewise you only have two genotype classes in a RIL population (AA and BB). 
#So if in the data file you marked all of your BB individuals with an H (ie. made your RILs look like a BC) 
#this next step will tell R/QTL that those *H* genotypes are actually BB and not hets, and account 
#for the map expansion that occurs in a RIL population. Use convert2risib(cross) for RILs developed by sib-mating. 

##########################################################################################
#cross<-convert2riself(cross)
#cross<-jittermap(cross)
#cross<-calc.genoprob(cross, map.function="kosambi",error.prob = 0.025)
#newmap=est.map(cross, map.function="kosambi", n.cluster=4, tol=0.01,
 #              maxit=1000)
#cross= replace.map(cross, newmap)

##########################################################################################

#(4) Slightly adjust marker positions to to avoid identical positions for markers on different chromosomes.

##########################################################################################
cross<-jittermap(cross)
##########################################################################################

#(5) Calculate the underlying genotype probabilities using a kosambi recombination model.

##########################################################################################
cross<-calc.genoprob(cross, map.function="kosambi",error.prob = 0.025)
##########################################################################################


############################
####CHECK PHENOTYPE DATA####
############################


##########################################################################################
#BE CAREFUL! THE FOLLOWING CODE ASSUMES YOUR FIRST PHENOTYPE IS SIMPLY LINE IDs AND AS SUCH IT IS 
#NOT CONSIDERED. IF THIS IS NOT THE CASE IN YOUR DATA FILE, ADD A DUMMY PHENOTYPE AS THE FIRST PHENOTYPE IN YOUR DATA.
##########################################################################################

#(6) Generate Histogram plots to visualize each phenotypic distribution.

##########################################################################################
pdf(file="1_Phenotype Histograms.pdf", width=11, height=8.5)
for(i in 1:(nphe(cross)-1)) {
    plot.pheno(cross, pheno.col=i)
}
dev.off()	
##########################################################################################

#(7) Generate .csv file for Shapiro-wilk normality test. Low p-values indicate departures from normality.

##########################################################################################
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
write.csv(file="2_Normality test p-values.csv", norm2)
##########################################################################################

#(8) Generate pdf of plot to check for batch effects.

##########################################################################################
pdf(file="3_Batch effects.pdf", width=11, height=8.5)
for(i in 1:(nphe(cross)-1)){
    par(mfrow=c(1,2), las=1, cex=0.8)
    means<-apply(cross$pheno[i], 1, mean)
    plot(means)
    plot(sample(means), xlab="Random Index", ylab="means", main=colnames(cross$pheno[i]))
}
dev.off()
##########################################################################################


###########################
####CHECK GENOTYPE DATA####
###########################


#(9) Generate a graphical genotype.

##########################################################################################
pdf(file="3,5_Graphical Genotype.pdf", width=11, height=8.5)
geno.image(cross)
dev.off()
##########################################################################################

#(10) Generate a genetic map and save as a pdf. 
#Look for map expansion to identify markers whose map position may be in error.

##########################################################################################
pdf(file="4_Estimated Genetic Map.pdf", width=11, height=8.5)
par(cex=0.6)
plot.map(cross, show.marker.names=FALSE)
dev.off()
##########################################################################################

#(11) Generate a .csv file with p-values for a X^2 test for segregation distortion 
#i.e. departures from mendelian expectation.

##########################################################################################
sd<-geno.table(cross)
sd<-sd[ sd$P.value < 1e-5, ]
write.csv(file="5_Chi-square for segregation distortion.csv", sd)
##########################################################################################

#(12) Generate a pdf of a histogram to compare genotypes for each pair of individuals 
#and identify pairs that have unusually similiar genotypes. i.e. sample mix up.

##########################################################################################
pdf(file="6_Histogram of genotype comparisons.pdf")
genotype.comparisons<-comparegeno(cross)
hist(genotype.comparisons, breaks=200, xlab="Proportion of identical genotypes")
rug(genotype.comparisons)
dev.off()
##########################################################################################

#(13) Generate csv file to identify which individuals are unusually similar and their proportion of similiarity

##########################################################################################
#cg.high<-which(genotype.comparisons>0.95, arr.ind=TRUE)
#proportion<-{}
#for(i in 1:nrow(cg.high)){
 #   x<-genotype.comparisons[cg.high[i,1], cg.high[i,2]]
  #  proportion<-c(proportion, x)
#}
#Genotype1<-cg.high[,1]
#Genotype2<-cg.high[,2]
#cg.high<-cbind(Genotype1, Genotype2, proportion)
#write.csv(file="7_Unusually similiar genotypes.csv", cg.high)
##########################################################################################

#(14) Identify markers where >15% of the population is NA and export as csv. 
#To change this threshold change the value of x<-x*0.15 to the desired tolerance.

##########################################################################################
x<-(nind(cross))
x<-x*0.15
missing.ind<-nmissing(cross, what="mar")
missing.ind<-missing.ind[missing.ind > x]
nmi<-"Number of Individuals without marker data"
missing.ind<-c(nmi, missing.ind)
write.csv(file="8_Poorly typed markers.csv", missing.ind)
##########################################################################################

#(15) Identify Individuals where >15% of the markers is NA and export as csv. 
#To change this threshold change the value of x<-x*0.15 to the desired tolerance.

##########################################################################################
x<-(totmar(cross))
x<-x*0.15
number.missing.mar<-nmissing(cross, what="ind")
Line.ID<-1:nind(cross)
missing.mar<-cbind(Line.ID, number.missing.mar)
missing.mar<-missing.mar[missing.mar[,2]>x,]
write.csv(file="9_Poorly typed individuals.csv", missing.mar)
##########################################################################################

############################################
####INITIAL GENOME SCAN FOR ADDITIVE QTL####
############################################

##########################################################################################
#YOU MAY CONSIDER SPLITTING NORMAL FROM NON-NORMAL AND/OR BINARY PHENOTYPES AND RUNNING SEPARATELY.

#FOR BINARY TRAITS CHANGE THE model="normal" ARGUMENTS TO model="binary" FOR STEPS #(16), #(17), AND #(25) ONLY. 

#FOR CONSIDERATION OF BINARY TRAITS THE CODE STILL REQUIRES LINE IDs AS YOUR FIRST PHENOTYPE, BUT WILL GIVE AN ERROR 
#AS YOUR LINE IDs ARE NOT BINARY. FOR CONSIDERATION OF BINARY TRAITS, INCLUDE A DUMMY PHENOTYPE CONSISTING ENTIRELY 
#OF ZEROS RATHER THAN LINE IDs AS THE FIRST PHENOTYPE. 

#FOR CONSIDERATION OF NON-PARAMETRIC TRAITS CHANGE model="normal" TO model="np" ONLY FOR STEPS #(16), AND #(17).

##########################################################################################

#(16) Perform single QTL Scan for all phenotypes using Haley-Knott Regression. 
#This may produce a warning if there are individuals with missing phenotype data. 
#If this occurs it does not change the output of the code. 

#*NOTE*# THIS STEP MAY GIVE A NON-CONVERGENCE ERROR FOR BINARY TRAITS.

##########################################################################################
cross.sc1<-scanone(cross, pheno.col=1:(nphe(cross)-1), method="em", use="all.obs", model="normal")
##########################################################################################

#(17) Calculate LOD significance thresholds based on 1000 permutations<--THIS COULD TAKE A WHILE. 
#If you have a multi-core processor, you may consider invoking the SNOW package ' library(snow) ' 
#and specifying the argument n.cluster with the number of parallel computations you'd like to do 
#(usually the number of chromosomes in the genome). 

#*NOTE*# THIS STEP MAY GIVE A NON-CONVERGENCE ERROR FOR BINARY TRAITS.

##########################################################################################
cross.sc1.perms<-scanone(cross, pheno.col=(nphe(cross)-1), method="em", n.perm=1000, verbose=TRUE, model="normal", n.cluster=4)
##########################################################################################

#(18) Produce a text file indicating the most signficant marker on each chromosme is extracted (they may or may not cross the significance threshold determined by the permutation test.), 
#the drop 1.5 LOD intervals, LODs, and pvalues for QTL detected.

##########################################################################################
sum<-summary(cross.sc1, threshold=0, perms=cross.sc1.perms, pvalues=TRUE, format="tabByCol", ci.function="lodint", drop=1.5, expandtomarkers=TRUE)
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
#pdf(file="QTL Plots.pdf", width=11, height=8.5)
#for(i in 1:(nphe(cross)-1)){
 #   plot(cross.sc1, lodcolumn=i, lwd=1.5, gap=0, bandcol="gray70", incl.markers=TRUE, main=colnames(cross$pheno[i]), xlab=c("Threshold for alpha=.05 using 1000 permutations", z[i]))
  #  add.threshold(cross.sc1, perms=cross.sc1.perms, alpha=0.05, lodcolumn=i, gap=0)
#}
#dev.off()
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
        xx<-capture.output(summary(cross.aq, perms=sub.perms, alpha=.05, pvalues=TRUE, format="tabByCol", ci.function="lodint", drop=1.5, expandtomarkers=TRUE))
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
        chr.sw<-{}
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
        
        #(29) Writes a file containing the peak marker and the flanking markers representing the 1.5 LOD interval of each QTL.
        
        ##########################################################################################
        Q<-"Q"
        space<-" "
        for (i in 1:rqtl$n.qtl){
            interval<-capture.output(lodint(rqtl, qtl.index=i, drop=1.5, expandtomarkers=TRUE))
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
