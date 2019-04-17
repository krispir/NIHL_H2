
#----------------------------------------------------------------------------------------------------
# Load necessary libraries
#----------------------------------------------------------------------------------------------------

library("xcms")
library("CAMERA")

#----------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------
# CAMERA FIX FOR IGRAPH/XCMS GROUPS COLLISION (from https://support.bioconductor.org/p/69414/)
# THIS IS NECESSARY TO RUN IN SOME VERSIONS OF R FOR CAMERA TO WORK...
#----------------------------------------------------------------------------------------------------

imports = parent.env(getNamespace("CAMERA"))
unlockBinding("groups", imports)
imports[["groups"]] = xcms::groups
lockBinding("groups", imports)

#----------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------
# Generate an xcmsSet from CDF files
#----------------------------------------------------------------------------------------------------

# Point to folder containing CDF files
setwd("E:/Forskningsprojekt/Cisplatin2/1701_Cisplatin2_Perilymph/CDF/POS/study")

# Generate xcms object from CDF files
xset <- xcmsSet(method="centWave", # Peak picking algorithm ("centWave" is recommended for UPLC-HRMS data)
                peakwidth=c(5, 50), # Guideline peakwidth domain for the centWave algorithm. Usually set very wide.
                ppm=20, # ppm cutoff for ROI generation
                noise=0, # IC cutoff for noise elimination. A lower value means longer computation times.
                snthresh=10, # S/N threshold for peak picking.
                mzdiff=0.02, # m/z window used for peak picking. (I am unsure of the differences to ppm, maybe complementary)
                prefilter=c(3, 100), # 
                mzCenterFun="wMean", 
                integrate=1, 
                fitgauss=TRUE, # Generate fit data 
                verbose.columns=FALSE, 
                nSlaves=2)
# Group features across samples
xset.g1 <- group(xset, 
                 method="density", 
                 bw=12.4, 
                 mzwid=0.02, 
                 minfrac=1, 
                 minsamp=1, 
                 max=50)
# Perform retention time correction
xset.c <- retcor(xset.g1, 
                 method="obiwarp", 
                 plottype="none", 
                 distFunc="cor_opt", 
                 profStep=1, 
                 center=3, 
                 response=1, 
                 gapInit=0.4, 
                 gapExtend=2.7, 
                 factorDiag=2, 
                 factorGap=1, 
                 localAlignment=0)
# Group features across samples
xset.g2 <- group(xset.c, 
                 method="density", 
                 bw=12.4, 
                 mzwid=0.02, 
                 minfrac=1, 
                 minsamp=1, 
                 max=50)

# Perform fillPeaks on missing integrals
xset.f <- fillPeaks(xset.g2, nSlaves=2)

#----------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------
# Generate CAMERA data from xset
#----------------------------------------------------------------------------------------------------

# CAMERA POS data with Stansrup rules
file<-system.file("rules/rules_jan_pos.csv", package="CAMERA")
rules<-read.csv(file)

# Perform CAMERA on the dataset
xset.f.an <- xsAnnotate(xset.f, 
                        polarity="positive")
xset.f.anF <- groupFWHM(xset.f.an, 
                        perfwhm=0.6)
xset.f.anI <- findIsotopes(xset.f.anF, 
                           mzabs=0.01, 
                           intval="into", 
                           minfrac=0.2)
xset.f.anC <- groupCorr(xset.f.anI, 
                        cor_eic_th=0.8, 
                        cor_exp_th=0.9, 
                        calcIso=TRUE, 
                        calcCiS=TRUE, 
                        calcCaS=TRUE)
xset.f.anA <- findAdducts(xset.f.anC, 
                          ppm=5, 
                          mzabs=0.01, 
                          polarity="positive", 
                          rules=rules)
						  
#----------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------
# Save data structures
# REMEMBER TO CHANGE PATH BEFORE RUNNING THE SCRIPT!
#----------------------------------------------------------------------------------------------------
						  
# Point to output location

# !!! OBS!: CHANGE THIS PATH BEFORE RUNNING THE SCRIPT! !!!

setwd(paste("C:/Users/kripi733/Box Sync/Forskning/001 Projekt/15-BK5309 Cisplatin 2/", 
            "003 Results/Manus/Bullerdatan/Data och figurer/000 COMPLETE DATA TREATMENT PATHS/", 
			"01 Data pretreatment", sep=""))

# Save data
save(xset, file="xset.Rda")
save(xset.f, file="xset.f.Rda")
save(xset.f.anA, file="xset.f.anA.Rda")

#----------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------
# Create a peak list with feature names and camera data
#----------------------------------------------------------------------------------------------------

xset.f.anA.pT <- data.frame(cbind(name=groupnames(xset.f),
                                  getPeaklist(xset.f.anA,intval="into")),
                            row.names=NULL,check.names=FALSE)

# Get a complete sample list from the xset@phenoData
samples.exclqc <- as.character(dimnames(xset.f@phenoData)[[1]][xset.f@phenoData!="QC"])
samples.all <- as.character(dimnames(xset.f@phenoData)[[1]])

# Use the below commented code to manually verify that the ordering of the samples is the same in 
# samples.all and columns of peakTable (TRUE=same sample)

# for (i in 1:length(samples.all))
# {
#       if (grepl(samples.all[i],dimnames(xset.f.anA.pT)[[2]][c(16:110)][i]))
#       {
#             print(TRUE)
#       } else
#       {
#             cat(paste(samples.all[i],":",
#                       dimnames(xset.f.anA.pT)[[2]][c(16:(length(dimnames(xset.f.anA.pT)[[2]])-3))][i],
#                       "\n",sep=""))
#       }
# }

# Loop through and print the replacements as they are made for review
# This allows the manual verification that all replacements are correct.
for (i in 1:length(samples.all)) {
     print(paste(dimnames(xset.f.anA.pT)[[2]][c(16:(length(dimnames(xset.f.anA.pT)[[2]])-3))][i],
			" -> ",samples.all[i],sep=""))
     dimnames(xset.f.anA.pT)[[2]][c(16:(length(dimnames(xset.f.anA.pT)[[2]])-3))][i] <- samples.all[i]
}

#----------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------
# MFC normalize the dataset to an average of the QC samples
#----------------------------------------------------------------------------------------------------

# General formula for median fold change normalization to a representative sample
# t(t(alldata) / apply(alldata/apply(qcs,1,mean),2,median))
# where alldata = (m,n) and repsamp and repsamp is array of length m

# Make a list of QC sample IDs
# MODIFY THESE TO FIT YOUR DATASET! THIS LIST SHOULD BE THE LIST OF SAMPLE IDS FOR YOUR QC SAMPLES
samples.qc <- c("QC1","QC2","QC3","QC4","QC5","QC6","QC7","QC8",
                "QC9","QC10","QC11","QC12","QC13","QC14","QC15")

# Normalize according to mfc formula
# NOTE: Procedure has been validated by manual calculation of mfc normalized data in excel and 
# comparison to procedure output.
xset.mfc <- xset.f.anA.pT
xset.mfc[,samples.all] <- t(t(xset.mfc[,samples.all]) / 
							apply(xset.mfc[,samples.all]/apply(xset.mfc[,samples.qc],1,mean),2,median))

# Calculate %RSD across QC injections
xset.mfc <- data.frame(cbind(xset.mfc,
                             rsdqc=apply(xset.mfc[,samples.qc],1,sd) / 
								apply(xset.mfc[,samples.qc],1,mean)),check.names=FALSE)

#----------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------
# Filter features
#----------------------------------------------------------------------------------------------------

# Cutoff values for filtering of features
# CHANGE TO FIT YOUR DATA!
rtmin <- 60
rtmax <- 960
maxrsd <- 0.3

# Filter features with 60<rt<960 and rsdqc>0.3
xset.mfc.filt <- xset.mfc[xset.mfc[,"rt"]>rtmin & xset.mfc[,"rt"]<rtmax & xset.mfc[,"rsdqc"]<=maxrsd,]
xset.mfc.filt[,1] <- factor(xset.mfc.filt[,1],levels=xset.mfc.filt[,1])

#----------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------
# Output normalized datasets
#----------------------------------------------------------------------------------------------------

write.table(xset.mfc, file="xset.mfc.csv", sep=";", dec=".", row.names=TRUE, col.names=TRUE, na="NA")
write.table(xset.mfc.filt, file="xset.mfc.filt.csv", sep=";", 
			dec=".", row.names=TRUE, col.names=TRUE, na="NA")

#----------------------------------------------------------------------------------------------------