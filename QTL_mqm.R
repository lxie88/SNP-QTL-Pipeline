###############################################################################
# 1-LOD QTL Mapping Script (HK method)
#
# This script performs QTL detection, permutation testing, and model refinement 
# for multiple phenotypes in a single cross object using the R/qtl package.
# 
# Key steps include:
#   - Reading data and setting up the cross object
#   - Checking genotype quality and performing error-prob calculations
#   - Running normality checks and diagnosing batch effects
#   - Performing initial genome scans (scanone) with permutations
#   - Extracting QTL peaks and refining QTL positions
#   - Using stepwise model selection (stepwiseqtl) to find the best QTL model
#   - Extracting 1-LOD support intervals for each QTL
#   - Plotting and exporting final results, such as LOD profiles and effect plots
#   - Generating genotype class means and standard errors for each QTL
#
# By default, this script:
#   - Uses the "hk" (Haley-Knott regression) method
#   - Simulates genotypes via sim.geno with 500 draws
#   - Calculates QTL intervals by dropping 1 LOD from the peak
#   - Analyzes every phenotype in cross$pheno except the last column (e.g., if it's an ID)
#
# Author: [Your Name or Lab]
# Date: [Date or Version]
###############################################################################

### Housekeeping: set working directory and load required packages
setwd("/users/limengxie/Desktop/JASHS_mqm/")
# install.packages("qtl")  # uncomment if needed

library(qtl)
library(rlecuyer)
library(snow)

###############################################################################
# Step 1: Read in the cross data
# 
# Here, we read in a csv file with genotypic/phenotypic data. We then apply 
# jittermap() to slightly perturb marker positions if some share the exact same 
# genetic position.
###############################################################################
cross <- read.cross("csv", 
                    "/users/limengxie/Desktop/JASHS_mqm/", 
                    "Mapping181_11.csv",
                    estimate.map = FALSE)
cross <- jittermap(cross)

###############################################################################
# Optional (commented out): Estimate an optimal genotyping error probability
# 
# This code chunk systematically tries a range of error.prob values to see which 
# one yields the maximum log likelihood. The user can then pick the rate that 
# best fits the data. We plot the LOD difference relative to the max (peak).
###############################################################################
# loglik <- err <- c(0.005, 0.01,0.015, 0.02,0.025,0.03,0.035,0.04, 0.045,0.05)
# for(i in seq(along=err)) {
#   cat(i, "of", length(err), "\n")
#   thkpmap <- est.map(cross, error.prob=err[i])
#   loglik[i] <- sum(sapply(thkpmap, attr, "loglik"))
# }
# lod <- (loglik - max(loglik))/log(10)
# plot(err, lod, xlab="Genotyping error rate",
#      xlim=c(0.0, 0.05),
#      ylab=expression(paste(log[10], " likelihood")))

###############################################################################
# Step 2: Set error probability and calculate genotype probabilities
# 
# We choose an error.prob (e.g., 0.025) and recalculate genotype probabilities 
# for each marker location, using Kosambi mapping function.
###############################################################################
cross <- calc.genoprob(cross, 
                       map.function = "kosambi",
                       error.prob   = 0.025)

###############################################################################
# (Optional) Normality checks (Shapiro-Wilk) on each phenotype 
# 
# This chunk is commented out but can be used to gauge whether data are normal.
###############################################################################
# norm <- {}
# norm2 <- {}
# for(i in 1:(nphe(cross)-1)) {
#   x <- cross$pheno[i]
#   x <- x[,1]
#   y <- shapiro.test(x)
#   ShapiroWilk.pvalue <- y$p.value
#   Phenotype <- colnames(cross$pheno[i])
#   norm <- cbind(Phenotype, ShapiroWilk.pvalue)
#   norm2 <- rbind(norm2, norm)
# }
# write.csv(file="Normality test p-values.csv", norm2)

###############################################################################
# (Optional) Batch effects 
# 
# Quick diagnostic: plots the mean of each individual’s phenotype. This chunk 
# can help identify outliers or batch effects across individuals.
###############################################################################
pdf(file="Batch effects.pdf", width=11, height=8.5)
for(i in 1:(nphe(cross)-1)) {
  par(mfrow=c(1,2), las=1, cex=0.8)
  means <- apply(cross$pheno[i], 1, mean)
  plot(means)
  plot(sample(means), xlab="Random Index", ylab="means", 
       main=colnames(cross$pheno[i]))
}
dev.off()

###############################################################################
# Step 3: Graphical Genotype representation
# 
# This provides a heatmap-like view of genotypes across the entire cross. 
# Large stretches of missing data or suspected errors can be visually identified.
###############################################################################
pdf(file="Graphical Genotype.pdf", width=11, height=8.5)
geno.image(cross)
dev.off()

###############################################################################
# Step 4: Save a plot of the genetic map
###############################################################################
pdf(file="Estimated Genetic Map.pdf", width=11, height=8.5)
par(cex=0.6)
plot.map(cross, show.marker.names=FALSE)
dev.off()

###############################################################################
# Step 5: Initial QTL scans (scanone) + Permutations
# 
# We scan for QTL across the genome using the “hk” method and then run 1000 
# permutations to assess significance thresholds for each phenotype.
###############################################################################
cross.sc1 <- scanone(cross, 
                     pheno.col = 1:(nphe(cross)-1),
                     method    = "hk",
                     use       = "all.obs",
                     model     = "normal")

cross.sc1.perms <- scanone(cross,
                           pheno.col = 1:(nphe(cross)-1),
                           method    = "hk",
                           n.perm    = 1000,
                           verbose   = TRUE,
                           model     = "normal",
                           n.cluster = 4)

###############################################################################
# Step 6: Summarize initial QTL hits
# 
# We create a summary (with threshold=0 to list all peaks) that includes p-values.
# The 1-LOD support interval is obtained by setting drop=1 in the summary.
###############################################################################
sum <- summary(cross.sc1,
               threshold      = 0,
               perms          = cross.sc1.perms,
               pvalues        = TRUE,
               format         = "tabByCol",
               ci.function    = "lodint",
               drop           = 1,
               expandtomarkers= TRUE)

space <- " "

for(i in 1:(nphe(cross)-1)){
  x <- capture.output(sum[[i]])
  cat(colnames(cross$pheno[i]), 
      file="10_Initial QTL hits by phenotype.txt", 
      sep="\n", 
      append=TRUE)
  cat(x, 
      file="10_Summary of top hits by phenotype.txt", 
      sep="\n", 
      append=TRUE)
  cat(space, 
      file="10_Summary of top hits by phenotype.txt", 
      sep="\n", 
      append=TRUE)
}

###############################################################################
# Step 7: Export LOD scores at every marker to a text file
###############################################################################
lods <- capture.output(cross.sc1)
cat(lods, file="11_LOD scores for every marker.txt", sep="\n")

###############################################################################
# Step 8: Plot QTL profiles with alpha=0.05 threshold lines
###############################################################################
z <- summary(cross.sc1.perms, alpha = .05)  # LOD threshold at alpha=.05

pdf(file="QTL Plots.pdf", width=11, height=8.5)
for(i in 1:(nphe(cross)-1)) {
  plot(cross.sc1, lodcolumn = i, lwd = 1.5, gap = 0, 
       bandcol = "gray70", incl.markers = TRUE,
       main     = colnames(cross$pheno[i]),
       xlab     = c("Threshold for alpha=.05 using 1000 permutations", z[i]))
  add.threshold(cross.sc1, 
                perms     = cross.sc1.perms, 
                alpha     = 0.05, 
                lodcolumn = i, 
                gap       = 0)
}
dev.off()

###############################################################################
# Step 9: Export a CSV file listing the map positions for each marker
# 
# Because we iterate over each chromosome’s marker map, we might see warnings 
# about column names. That’s normal if the table is appended in a loop.
###############################################################################
newmap <- pull.map(cross)
for(i in 1:length(names(newmap))) {
  snps <- names(newmap[[i]])
  gm   <- c(snps, newmap[[i]])
  gm2  <- matrix(gm, ncol=2)
  write.table(file="13_Genetic Map Positions.csv", 
              sep=",", 
              append=TRUE, 
              gm2)
}

###############################################################################
# Step 10: Run sim.geno (needed for effect plotting)
# 
# We simulate genotypes with 500 draws (ndraw=500), error prob of 0.025, 
# step=1, and Kosambi mapping. This is required for some advanced plotting 
# features (e.g., effectplot).
###############################################################################
cross2 <- sim.geno(cross, 
                   n.draws     = 500, 
                   error.prob  = 0.025, 
                   step        = 1, 
                   map.function= "kosambi")

###############################################################################
# Steps 11–20: Model refinement in a loop over phenotypes
# 
# The code below (in one big for-loop) automatically:
#   - Identifies initial QTL peaks from the summary
#   - Performs addqtl() to look for additional QTL
#   - Uses stepwiseqtl() for forward/backward model selection
#   - Refines QTL positions (refineqtl)
#   - Exports 1-LOD intervals and final ANOVA summaries
#   - Plots marker effects and writes genotype class means/SE
#
# The “if(is.na(chr[1])==FALSE)” block ensures we only proceed if we have 
# detected at least one QTL peak above threshold. Otherwise, we note that 
# there were no significant QTL for that phenotype.
###############################################################################
for(k in 1:(nphe(cross)-1)){
  
  # (11) Identify which phenotype is being analyzed in this iteration
  phenotype <- k
  pheno     <- colnames(cross$pheno[phenotype])
  
  # (12) Make a QTL object from the initial summary ("10_Initial QTL hits...") 
  # for the chosen phenotype
  chromo <- sum[[pheno]]
  chr    <- {}
  pos    <- {}
  
  for(i in 1:nrow(chromo)){
    chr1 <- chromo[i,1]
    chr2 <- as.numeric(as.character(chr1))
    chr  <- c(chr, chr2)
  }
  for(i in 1:nrow(chromo)){
    pos1 <- chromo[i,2]
    pos  <- c(pos, pos1)
  }
  
  # Only proceed if the initial summary found QTL
  if(is.na(chr[1]) == FALSE) {
    qtl <- makeqtl(cross, chr=chr, pos=pos, what="prob")
    createqtl <- paste("Q", 1:qtl$n.qtl, sep="")
    formula   <- as.formula(paste("y ~ ", paste(createqtl, collapse="+")))
    
    # (13) addqtl(): Search for additional linked QTL, conditioning on the 
    # QTL already found
    cross.aq   <- addqtl(cross, 
                         pheno.col = phenotype,
                         qtl       = qtl, 
                         formula   = formula, 
                         method    = "hk", 
                         model     = "normal")
    sub.perms  <- subset(cross.sc1.perms, lodcolumn=phenotype)
    xx         <- capture.output(summary(cross.aq, 
                                         perms=sub.perms, 
                                         alpha=.05, 
                                         pvalues=TRUE, 
                                         format="tabByCol", 
                                         ci.function="lodint", 
                                         drop=1, 
                                         expandtomarkers=TRUE))
    xxy        <- c(pheno, xx)
    cat(xxy, file="14_Additional QTL hits by phenotype.txt", sep="\n", append=TRUE)
    cat(space, file="14_Additional QTL hits by phenotype.txt", sep="\n", append=TRUE)
    
    # (14) Incorporate any newly found QTL into our QTL object
    sum.aq  <- summary(cross.aq,
                       perms          = sub.perms,
                       alpha          = .05,
                       pvalues        = TRUE,
                       format         = "tabByCol",
                       ci.function    = "lodint",
                       drop           = 0.95,
                       expandtomarkers= TRUE)
    chr.aq  <- {}
    pos.aq  <- {}
    for(i in 1:nrow(sum.aq$lod)){
      chr.aq1 <- sum.aq$lod[i,1]
      chr.aq2 <- as.numeric(as.character(chr.aq1))
      chr.aq  <- c(chr.aq, chr.aq1)
    }
    # If no new QTL, keep old ones
    if(is.na(chr.aq)==TRUE) { chr.aq <- chr } else { chr.aq <- c(chr, chr.aq) }
    
    for(i in 1:nrow(sum.aq$lod)){
      pos1.aq <- sum.aq$lod[i,2]
      pos.aq  <- c(pos.aq, pos1.aq)
    }
    if(is.na(pos.aq)==TRUE) { pos.aq <- pos } else { pos.aq <- c(pos,pos.aq) }
    
    # (15) Update qtl object and build formula
    qtl        <- makeqtl(cross, chr=chr.aq, pos=pos.aq, what="prob")
    createqtl  <- paste("Q", 1:qtl$n.qtl, sep="")
    formula    <- as.formula(paste("y ~ ", paste(createqtl, collapse="+")))
    
    # (16) Stepwise model selection 
    # This chunk uses user-defined permutations to set the penalty threshold
    pen     <- summary(sub.perms)
    cross.sw<- stepwiseqtl(cross, 
                           pheno.col       = phenotype, 
                           qtl             = qtl, 
                           formula         = formula, 
                           method          = "hk", 
                           penalties       = pen, 
                           model           = "normal", 
                           additive.only   = TRUE)
    swQTL   <- capture.output(print(cross.sw))
    swQTL2  <- c(pheno, swQTL)
    cat(swQTL2, file="15_Additional QTL hits from stepwise analysis.txt", sep="\n", append=TRUE)
    cat(space, file="15_Additional QTL hits from stepwise analysis.txt", sep="\n", append=TRUE)
    
    # (17) Refine QTL positions from stepwise results
    sum.sw   <- summary(cross.sw)
    chr.sw   <- {}
    pos.sw   <- {}
    for(i in 1:nrow(sum.sw)){
      chr.sw1 <- sum.sw[i,2]
      chr.sw2 <- as.numeric(as.character(chr.sw1))
      chr.sw  <- c(chr.sw, chr.sw2)
    }
    for(i in 1:nrow(sum.sw)){
      pos1.sw <- sum.sw[i,3]
      pos.sw  <- c(pos.sw, pos1.sw)
    }
    qtl2    <- makeqtl(cross, chr=chr.sw, pos=pos.sw, what="prob")
    createqtl<- paste("Q", 1:qtl2$n.qtl, sep="")
    formula <- as.formula(paste("y ~ ", paste(createqtl, collapse="+")))
    rqtl    <- refineqtl(cross, 
                         pheno.col=phenotype, 
                         qtl     =qtl2, 
                         method  ="hk",
                         model   ="normal")
    
    # (18) Save the 1-LOD interval markers around each refined QTL
    Q <- "Q"
    for(i in 1:rqtl$n.qtl){
      interval    <- capture.output(lodint(rqtl, 
                                          qtl.index       = i, 
                                          drop            = 1, 
                                          expandtomarkers = TRUE))
      q           <- paste(Q, i, sep="")
      interval.new<- c(pheno, q, interval, space)
      cat(interval.new, 
          file="16_Summary of Final QTL Intervals.txt", 
          sep="\n", 
          append=TRUE)
    }
    
    # (19) Fit the final QTL model and export an ANOVA summary 
    cross.ests <- fitqtl(cross, 
                         pheno.col=phenotype, 
                         qtl     =rqtl, 
                         formula =formula, 
                         method  ="hk", 
                         dropone =TRUE, 
                         get.ests=TRUE, 
                         model   ="normal")
    ests       <- capture.output(summary(cross.ests))
    ests       <- c(pheno, ests)
    write(ests, file="17_ANOVA results and QTL effect estimates.txt", 
          sep="\n", 
          append=TRUE)
    
    # (20) Generate effect plots + means & SE for each genotype class
    # One PDF file per QTL, plus a text file with genotype means and SE
    for(i in 1:length(chr.sw)) {
      b    <- paste("Q", i, ".pdf", sep="")
      file <- paste("18_Marker effect plots", pheno, b)
      
      mar  <- find.marker(cross2, chr=chr.sw[i], pos=pos.sw[i])
      pdf(file=file, width=11, height=8.5)
      plotPXG(cross, marker=mar, pheno.col=phenotype)
      dev.off()
      
      phenoqtl  <- paste(pheno, b)
      pheno.eff <- c(phenoqtl, space)
      means     <- effectplot(cross2, 
                              pheno.col=phenotype, 
                              mname1   =mar,
                              draw     =FALSE)
      cat(pheno.eff, 
          file="19_means and SE.txt", 
          sep="\n", 
          append=TRUE)
      write.table(means$Means, 
                  file     ="19_means and SE.txt", 
                  sep      =",", 
                  col.names="Means", 
                  row.names=TRUE, 
                  append   =TRUE)
      cat(space, 
          file="19_means and SE.txt", 
          sep="\n", 
          append=TRUE)
      write.table(means$SEs, 
                  file     ="19_means and SE.txt", 
                  sep      =",", 
                  col.names="Standard Error", 
                  row.names=TRUE, 
                  append   =TRUE)
      cat(space, 
          file="19_means and SE.txt", 
          sep="\n", 
          append=TRUE)
    }
    
  } else {
    # If no QTL found above threshold for this phenotype
    null <- "\tThere were no LOD peaks above the threshold"
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
    cat(space, file="19_means and SE.txt", sep="\n", append=TRUE)
  }
} # End for-loop over phenotypes
###############################################################################
