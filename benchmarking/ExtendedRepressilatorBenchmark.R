library(sRACIPE)
library(microbenchmark)

Extended_Repressilator <- data.frame(Source = c("1", "1", "1", "2", "2", "2", "3", "3",
                                                "4", "4", "4", "5", "5", "5", "6", "6",
                                                "7", "7", "7", "8", "8", "8", "9", "9"),
                                     Target = c("2", "3", "4", "1", "3", "5", "1", "6",
                                                "5", "6", "7", "4", "6", "8", "4", "9",
                                                "8", "9", "1", "7", "9", "2", "7", "3"),
                                     type = c(1, 1, 2, 1, 1, 2, 1, 2, 1, 1, 2,
                                              1, 1, 2, 1, 2, 1, 1, 2, 1, 1, 2,
                                              1, 2))

iterCounts = c(seq(from=1, to=20, by=1), seq(from=30, to = 250, by=10))
outputDir <- file.path(getwd(), "extendedRepressilator")
setwd(outputDir)

#integrateStepSize and numStepsConverge are chosen as default parameters in cellcycle

#iterResultsExtended1 <- lapply(iterCounts, function(n) {
#  microbenchmark(sracipeSimulate(circuit = Extended_Repressilator, nIC = 3, numConvergenceIter = n,
#                                 integrateStepSize = 0.01, numStepsConverge = 1000, numModels = 100), times = 15L)
#
#})

iterResultsExtended1 <- lapply(iterCounts, function(n) {
  i <- 0

  microbenchmark(
    {
      SEName <- paste0(n, "extendedRepressilator", ((i)%%15)+1)
      i <<- i + 1
      obj = sracipeSimulate(circuit = Extended_Repressilator, nIC = 3, numConvergenceIter = n,
                            integrateStepSize = 0.01, numStepsConverge = 1000, numModels = 100)
      save(obj, file = paste0(SEName, ".RDA"))
    },
    times = 15L
  )

})
save(iterResultsExtended1, file = "iterResultsExtended1.RDA")
