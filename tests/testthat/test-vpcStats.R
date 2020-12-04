context("Virtual Predictive Check statistics")

test_that(".addDoseRelTime returns the same dataset with an additional column", {
  data <- data.frame(id = c(rep(1, 6), rep(2, 6)), time = c(seq_len(6), seq_len(6)))
  mockery::stub(.addDoseRelTime, "getDoseInformation",
                data.frame(id = c(1, 1, 2, 2), time = c(1, 4, 1, 4)))
  mockery::stub(.addDoseRelTime, ".getSubjocc", c("id"))
  datawithdose <- .addDoseRelTime(data)
  expect_equal(nrow(datawithdose), nrow(data))
  expect_equal(ncol(datawithdose), ncol(data) + 1)
  expect_equal(subset(datawithdose, select = c("id", "time")), data)
  expect_equal(datawithdose$time_rel, rep(c(0, 1, 2, 0, 1, 2), 2))
})

test_that(".addDoseRelTime when number of dose different per subjocc", {
  data <- data.frame(id = c(rep(1, 6), rep(2, 6)), time = c(seq_len(6), seq_len(6)))
  mockery::stub(.addDoseRelTime, "getDoseInformation",
                data.frame(id = c(1, 2, 2), time = c(1, 1, 4)))
  mockery::stub(.addDoseRelTime, ".getSubjocc", c("id"))
  datawithdose <- .addDoseRelTime(data)
  expect_equal(nrow(datawithdose), nrow(data))
  expect_equal(ncol(datawithdose), ncol(data) + 1)
  expect_equal(subset(datawithdose, select = c("id", "time")), data)
  expect_equal(datawithdose$time_rel, c(0, 1, 2, 3, 4, 5, 0, 1, 2, 0, 1, 2))
})

test_that(".getObservationData returns mlx.getObservationData in case of discrete data", {
  mockery::stub(.getObservationData, "mlx.getObservationInformation", c(a = 2, b = 3))
  expect_equal(.getObservationData("a", "discrete"), 2)
  expect_equal(.getObservationData("b", "discrete"), 3)
})

test_that(".getObservationData returns mlx.getObservationData in case of event data", {
  mockery::stub(.getObservationData, "mlx.getObservationInformation", c(a = 2, b = 3))
  expect_equal(.getObservationData("a", "event"), 2)
  expect_equal(.getObservationData("b", "event"), 3)
})

test_that(".getObservationData returns vpc chart data observation in case of continuous data", {
  tmpdir <- tempdir()
  directory <- paste0(tmpdir, "/ChartsData/VisualPredictiveCheck/")
  unlink(directory, recursive = TRUE)
  dir.create(directory, recursive = TRUE)
  obsName <- "a"
  data <- data.frame(id = c(rep(1, 6), rep(2, 6)), time = c(seq_len(6), seq_len(6)))
  write.csv(data, file = paste0(directory, obsName, "_observations.txt"), row.names = FALSE)
  mockery::stub(.getObservationData, "mlx.getProjectSettings", list(directory = tmpdir))
  expect_equal(.getObservationData(obsName, "continuous"), data)
})

test_that(".getObservationData returns observation with id column", {
  tmpdir <- tempdir()
  directory <- paste0(tmpdir, "/ChartsData/VisualPredictiveCheck/")
  unlink(directory, recursive = TRUE)
  dir.create(directory, recursive = TRUE)
  obsName <- "a"
  data <- data.frame(ID = c(rep(1, 6), rep(2, 6)), time = c(seq_len(6), seq_len(6)))
  write.csv(data, file = paste0(directory, obsName, "_observations.txt"), row.names = FALSE)
  mockery::stub(.getObservationData, "mlx.getProjectSettings", list(directory = tmpdir))
  expect_true("id" %in% names(.getObservationData(obsName, "continuous")))
})

test_that(".getObservationData returns downloaded vpc chart data observation in case of continuous data", {
  tmpdir <- tempdir()
  directory <- paste0(tmpdir, "/ChartsData/VisualPredictiveCheck/")
  unlink(directory, recursive = TRUE)
  obsName <- "a"
  data <- data.frame(id = c(rep(1, 6), rep(2, 6)), time = c(seq_len(6), seq_len(6)))
  mockery::stub(.getObservationData, "mlx.getProjectSettings", list(directory = tmpdir))
  mockery::stub(.getObservationData, "mlx.getScenario", list())
  mockery::stub(.getObservationData, "mlx.setScenario", NULL)
  mockery::stub(
    .getObservationData, "mlx.computeChartsData",
    function(exportVPCSimulations) {
      if(!exportVPCSimulations) {
        NULL
      } else {
        tmpdir <- tempdir()
        dir.create(directory, recursive = TRUE)
        write.csv(data, file = paste0(directory, obsName, "_observations.txt"), row.names = FALSE)
      }
    }
  )
  expect_equal(.getObservationData(obsName, "continuous"), data)
})

test_that(".getSimulationData returns vpc chart data simulations", {
  tmpdir <- tempdir()
  directory <- paste0(tmpdir, "/ChartsData/VisualPredictiveCheck/")
  unlink(directory, recursive = TRUE)
  dir.create(directory, recursive = TRUE)
  obsName <- "a"
  data <- data.frame(id = c(rep(1, 6), rep(2, 6)), time = c(seq_len(6), seq_len(6)))
  write.csv(data, file = paste0(directory, obsName, "_simulations.txt"), row.names = FALSE)
  mockery::stub(.getSimulationData, "mlx.getProjectSettings", list(directory = tmpdir))
  expect_equal(.getSimulationData(obsName), data)
})

test_that(".getSimulationData returns simulations with id column", {
  tmpdir <- tempdir()
  directory <- paste0(tmpdir, "/ChartsData/VisualPredictiveCheck/")
  unlink(directory, recursive = TRUE)
  dir.create(directory, recursive = TRUE)
  obsName <- "a"
  data <- data.frame(ID = c(rep(1, 6), rep(2, 6)), time = c(seq_len(6), seq_len(6)))
  write.csv(data, file = paste0(directory, obsName, "_simulations.txt"), row.names = FALSE)
  mockery::stub(.getSimulationData, "mlx.getProjectSettings", list(directory = tmpdir))
  expect_true("id" %in% names(.getSimulationData(obsName)))
})

test_that(".getSimulationData returns downloaded vpc chart data simulations", {
  tmpdir <- tempdir()
  directory <- paste0(tmpdir, "/ChartsData/VisualPredictiveCheck/")
  unlink(directory, recursive = TRUE)
  obsName <- "a"
  data <- data.frame(id = c(rep(1, 6), rep(2, 6)), time = c(seq_len(6), seq_len(6)))
  mockery::stub(.getSimulationData, "mlx.getProjectSettings", list(directory = tmpdir))
  mockery::stub(.getSimulationData, "mlx.getScenario", list())
  mockery::stub(.getSimulationData, "mlx.setScenario", NULL)
  mockery::stub(
    .getSimulationData, "mlx.computeChartsData",
    function(exportVPCSimulations) {
      if(!exportVPCSimulations) {
        NULL
      } else {
        tmpdir <- tempdir()
        dir.create(directory, recursive = TRUE)
        write.csv(data, file = paste0(directory, obsName, "_simulations.txt"), row.names = FALSE)
      }
    }
  )
  expect_equal(.getSimulationData(obsName), data)
})
