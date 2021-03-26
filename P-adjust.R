
sigPairsDF.T0 <- read.csv("sigPairsDF.T0.csv")
sigPairsDF.T21 <- read.csv("sigPairsDF.T21.csv")

sigPairsDF.T0$FDR <- p.adjust(sigPairsDF.T0$X6,method="BH")
sigPairsDF.T21$FDR <- p.adjust(sigPairsDF.T21$X6,method="BH")

write.csv(sigPairsDF.T0,"sigPairsDF.T0.FDR.csv")
write.csv(sigPairsDF.T21,"sigPairsDF.T21.FDR.csv")
