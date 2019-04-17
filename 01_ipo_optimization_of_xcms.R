#source("https://bioconductor.org/biocLite.R")
#biocLite("xcms")
#biocLite("CAMERA")
biocLite("plyr")

library("xcms")
library("CAMERA")
library("IPO")

## CAMERA FIX FOR IGRAPH/XCMS groups COLLISION (from https://support.bioconductor.org/p/69414/)
imports = parent.env(getNamespace("CAMERA"))
unlockBinding("groups", imports)
imports[["groups"]] = xcms::groups
lockBinding("groups", imports)


########################################
########################################
### OPTIMIZATION OF POSITIVE QC DATA ###
########################################
########################################

# Start positive
start.pos <- Sys.time()

# Point to CDF files
setwd("G:/Forskningsprojekt/Cisplatin2/1701_Cisplatin2_Perilymph/CDF/POS/optQC")

#?list.files
mzdatafiles <- list.files(recursive = TRUE, full.names = TRUE)
#mzdatafiles

peakpickingParameters <- getDefaultXcmsSetStartingParams('centWave')
#setting levels for min_peakwidth to 10 and 20 (hence 15 is the center point)
peakpickingParameters$min_peakwidth <- c(10,20) 
peakpickingParameters$max_peakwidth <- c(26,42)
#setting only one value for ppm therefore this parameter is not optimized
peakpickingParameters$ppm <- 20
resultPeakpickingPOS <- 
  optimizeXcmsSet(files = mzdatafiles, 
                  params = peakpickingParameters, 
                  nSlaves = 4, 
                  subdir = 'rsmDirectory')
optimizedXcmsSetObjectPOS <- resultPeakpickingPOS$best_settings$xset
save(optimizedXcmsSetObjectPOS, file="optimizedXcmsSetObjectPOS.Rda")

retcorGroupParameters <- getDefaultRetGroupStartingParams()
retcorGroupParameters$profStep <- 1
resultRetcorGroupPOS <-
  optimizeRetGroup(xset = optimizedXcmsSetObjectPOS, 
                   params = retcorGroupParameters, 
                   nSlaves = 4, 
                   subdir = "rsmDirectory")

writeRScript(resultPeakpickingPOS$best_settings$parameters, 
             resultRetcorGroupPOS$best_settings, 
             nSlaves=4)

# End positive
end.pos <- Sys.time()
time.pos <- end.pos - start.pos
time.pos

########################################
########################################
### OPTIMIZATION OF NEGATIVE QC DATA ###
########################################
########################################

# Start positive
start.neg <- Sys.time()

# Point to CDF files
setwd("E:/Forskningsprojekt/Cisplatin2/1701_Cisplatin2_Perilymph/CDF/NEG/study")

mzdatafiles <- list.files(recursive = TRUE, full.names = TRUE)

peakpickingParameters <- getDefaultXcmsSetStartingParams('centWave')
#setting levels for min_peakwidth to 10 and 20 (hence 15 is the center point)
peakpickingParameters$min_peakwidth <- c(10,20) 
peakpickingParameters$max_peakwidth <- c(26,42)
#setting only one value for ppm therefore this parameter is not optimized
peakpickingParameters$ppm <- 20
resultPeakpickingPOS <- 
     optimizeXcmsSet(files = mzdatafiles[110:111], 
                     params = peakpickingParameters, 
                     nSlaves = 4, 
                     subdir = 'rsmDirectory')
optimizedXcmsSetObjectPOS <- resultPeakpickingPOS$best_settings$xset
save(optimizedXcmsSetObjectPOS, file="optimizedXcmsSetObjectPOS.Rda")

retcorGroupParameters <- getDefaultRetGroupStartingParams()
retcorGroupParameters$profStep <- 1
resultRetcorGroupPOS <-
     optimizeRetGroup(xset = optimizedXcmsSetObjectPOS, 
                      params = retcorGroupParameters, 
                      nSlaves = 4, 
                      subdir = "rsmDirectory")

writeRScript(resultPeakpickingPOS$best_settings$parameters, 
             resultRetcorGroupPOS$best_settings, 
             nSlaves=4)

# End positive
end.neg <- Sys.time()
time.neg <- end.neg - start.neg
time.neg

# Save workspace
save.image()
