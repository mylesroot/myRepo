# Import a BigWig file as an RleList object
library (rtracklayer)
library(genomation)
library(soGGi)

bw = ("~/Documents/Shaun Data/CRAC7_4th_repeat_NNNGCGCAGC_L5Ac_collapsed.bs1.forward.bw")

bigWig = import(bw)

covs = coverage(bigWig)

bigWig


# Window the exon1 / Intron1 regions

# For each window, get the counts for that region

# Windowing can take multiple bigWig files as input

# Take outputs from Soggi and then join them before the plotting.


# CREATE BETTER FILTERING FUNCTION AND ALLOW SOGGI REGION PLOT TO SELECT GENES/INVERSE FOR SEQLENGTHS

bigSoggi = regionPlot(bamFile = bw, testRanges = myBed, nOfWindows = 5, style = "percentOfRegion", format="bigwig", method = "bin")

soggi = regionPlot(bamFile = bw, testRanges = bedPlus, nOfWindows = 5, style = "percentOfRegion", format="bigwig", method = "bin",
                   distanceOutRegionStart = 0, distanceOutRegionEnd = 0, distanceAround = 100)












