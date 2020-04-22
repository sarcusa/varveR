# #fast single ensemble model
# transtable <- data.frame(codes = 1:6,overcount = c(0.05,0.10,0.30,0,-2,-3), undercount =  c(0.05,0.10,0.30,.50,-2,-3))
#
# tt <- transtable
# seq1 <- col172_sequence
# seq2 <- col173_sequence
#
#
#
# seq1 <- translateCodesToProbabilities(seq1,translationTable = tt,baselineProb = .05)
# seq2 <- translateCodesToProbabilities(seq2,translationTable = tt,baselineProb = .05)
#
# #generate thickness estimate
# sseq1 <- simulateOverAndUndercounting(seq1)
# g <- which(!is.na(sseq1$newThicks))
# sseq1$ensThick <- sseq1$newThicks[g]
# sseq1$tiePoints <- sseq1$newTiePoint[g]
#
# sseq2 <- simulateOverAndUndercounting(seq2)
# g <- which(!is.na(sseq2$newThicks))
# sseq2$ensThick <- sseq2$newThicks[g]
# sseq2$tiePoints <- sseq2$newTiePoint[g]
#
# #
# ensList <- list(sseq1,sseq2)
# #get all the tie points
# allTiePoints <- c("topOfCore",na.omit(unique(c(unique(sseq1$newTiePoint)),na.omit(unique(sseq2$newTiePoint)))))
#
#
# #simulate the varves
# modeledVarves <- fastVarveModel(ensList,allMarkerLayers = allTiePoints)
#
