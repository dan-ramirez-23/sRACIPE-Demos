library(sRACIPE)
library(microbenchmark)

data("CoupledToggleSwitchSA")

iterCounts = c(seq(from=1, to=20, by=1), seq(from=30, to = 250, by=10))
outputDir <- file.path(getwd(), "coupledToggleSwitch")
setwd(outputDir)

#iterResultsCoupled <- lapply(iterCounts, function(n) {
#  microbenchmark(sracipeSimulate(circuit = CoupledToggleSwitchSA, nIC = 3, numConvergenceIter = n,
#                                 integrateStepSize = 0.01, numStepsConverge = 1000, numModels = 100), times = 15L)
#})

i <- 0

iterResultsCoupled <- lapply(iterCounts, function(n) {

  microbenchmark(
    {
      SEName <- paste0(n, "CoupledToggleSwitch", ((i)%%15)+1)
      i <<- i + 1
      obj = sracipeSimulate(circuit = CoupledToggleSwitchSA, nIC = 3, numConvergenceIter = n,
                                     integrateStepSize = 0.01, numStepsConverge = 1000, numModels = 100)
      save(obj, file = paste0(SEName, ".RDA"))
    },
    times = 15L
  )

})


save(iterResultsCoupled, file = "iterResultsCoupled.RDA")
