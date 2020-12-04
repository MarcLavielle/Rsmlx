context("Virtual Predictive Check Binning")

test_that(".addBinsIndex returns error when innapropriate columnName", {
  data <- data.frame(id = c(rep(1, 6), rep(2, 6)), time = c(seq_len(6), seq_len(6)))
  expect_error(.addBinsIndex(data, columnName = "another_name"))
})

test_that(".addBinsIndex returns a input dataset when binsValue is NULL", {
  data <- data.frame(id = c(rep(1, 6), rep(2, 6)), time = c(seq_len(6), seq_len(6)))
  expect_equal(.addBinsIndex(data, columnName = "time", binsValue = NULL), data)
})

test_that(".addBinsIndex returns the same dataset with additional columns", {
  data <- data.frame(id = c(rep(1, 6), rep(2, 6)), time = c(seq_len(6), seq_len(6)))
  binsValue <- data.frame(bins = c(1, 2), bins_start = c(0, 3), bins_stop = c(3, 6))
  binName <- "toto"
  datawithbins <- .addBinsIndex(data, columnName = "time", binsValue = binsValue, binName = binName)
  expect_equal(nrow(datawithbins), nrow(data))
  expect_true(ncol(datawithbins) >= ncol(data))
  expect_equal(subset(datawithbins, select = c("id", "time")), data)
  expect_equal(datawithbins$toto, rep(c(1, 1, 1, 2, 2, 2), 2))
})

test_that(".getBins returns NULL when no settings and no category", {
  data <- data.frame(id = c(rep(1, 6), rep(2, 6)), time = c(seq_len(6), seq_len(6)))
  expect_null(.getBins(data, bycat = FALSE, settings = NULL))
})

test_that(".getBins create categories from data when categorical binning with no categories specified", {
  data <- c(1, 2, 3, 1, 1, 2, 5, 5, 3)
  categories <- c(1, 2, 3, 5)
  expectedBins <- data.frame(
    bins_start = categories, bins_stop = categories,
    bins_middle = categories, bins = seq_len(length(categories))
  )
  expect_equal(.getBins(data, bycat = TRUE), expectedBins)
})

test_that(".getBins returns bins based of categories in arg when categorical binning", {
  data <- c(1, 2, 3, 1, 1, 2, 5, 5, 3)
  categories <- c(5, 8, 12)
  expectedBins <- data.frame(
    bins_start = categories, bins_stop = categories,
    bins_middle = categories, bins = seq_len(length(categories))
  )
  expect_equal(.getBins(data, bycat = TRUE, categories = categories), expectedBins)
})

test_that(".getBins returns fixed bins when used fixed bins", {
  data <- c(1, 2, 3, 1, 1, 2, 5, 5, 3, 6)
  settings <- list(useFixedBins = TRUE, fixedBins = c(2, 5))
  expectedBins <- data.frame(
    bins_start = c(1, 2, 5), bins_stop = c(2, 5, 6),
    bins_middle = c(mean(c(1, 1, 1, 2, 2)), mean(c(3, 5, 5, 3)), mean(c(6))), bins = seq_len(3)
  )
  expect_equal(.getBins(data, settings = settings), expectedBins)
})

test_that(".getBins returns result of computeBins function when no fixed bins", {
  data <- c(1, 2, 3, 1, 1, 2, 5, 5, 3, 6)
  settings <- list(
    useFixedBins = FALSE, criteria = "toto", useFixedNbBins = "tata",
    fixedNbBins = "titi", binRange = "tonton", nbBinData = "tutu"
  )
  mockery::stub(.getBins, "mlx.computeBins", list(values = c(1, 7, 40), middles = c(10, 74)))
  expectedBins <- data.frame(
    bins_start = c(1, 7), bins_stop = c(7, 40),
    bins_middle = c(10, 74), bins = seq_len(2)
  )
  expect_equal(.getBins(data, settings = settings), expectedBins)
})

test_that("getBinsSettings returns error when inapropriate arguments", {
  expect_error(getBinsSettings(is.fixedBins = 5))
  expect_error(getBinsSettings(is.fixedBins = "toto"))
  expect_error(getBinsSettings(fixedBins = c(-1, 10, 12)))
  expect_error(getBinsSettings(fixedBins = c("a", 10, 12)))
  expect_error(getBinsSettings(criteria = "wrongcriteria"))
  expect_error(getBinsSettings(criteria = 7))
  expect_error(getBinsSettings(is.fixedNbBins = 5))
  expect_error(getBinsSettings(is.fixedNbBins = "toto"))
  expect_error(getBinsSettings(nbBins = "a"))
  expect_error(getBinsSettings(nbBins = TRUE))
  expect_error(getBinsSettings(nbBins = -5))
  expect_error(getBinsSettings(nbBins = 6.5))
  expect_error(getBinsSettings(binRange = c(-1, 10, 12)))
  expect_error(getBinsSettings(binRange = c("a", 10, 12)))
  expect_error(getBinsSettings(binRange = c(FALSE)))
  expect_error(getBinsSettings(binRange = c(6.53, 10, 12)))
  expect_error(getBinsSettings(nbBinData = c(-1, 10, 12)))
  expect_error(getBinsSettings(nbBinData = c("a", 10, 12)))
  expect_error(getBinsSettings(nbBinData = c(FALSE)))
  expect_error(getBinsSettings(nbBinData = c(6.53, 10, 12)))
})

test_that("getBinsSettings returns a binSettingsClass object which is a list", {
  expect_equal(class(getBinsSettings()), "binSettingsClass")
  expect_true(is.list(getBinsSettings()))
})

test_that("getBinsSettings always returns the same setting names", {
  expect_equal(names(getBinsSettings()), names(getBinsSettings(is.fixedBins = FALSE)))
})

test_that("getBinsSettings output take into account argument values", {
  expect_equal(getBinsSettings(fixedBins = c(5, 6, 10))$fixedBins, c(5, 6, 10))
  expect_equal(getBinsSettings(binRange = c(120, 140))$binRange, c(120, 140))

})

test_that("computeVpcBins returns error when inapropriate arguments", {
  mockery::stub(
    computeVpcBins,
    '.getBins',
    function(data, settings = NULL, bycat = FALSE, categories = c()) {
      if (bycat) {
        return(data.frame(type = "categorical"))
      } else {
        return(data.frame(type = "continuous"))
      }
    }
  )
  data <- c(1, 2, 3, 1, 1, 2, 5, 5, 3)
  expect_error(computeVpcBins())
  expect_error(computeVpcBins(data = c("a", 2, 1, 2)))
  expect_error(computeVpcBins(data = FALSE))
  expect_error(computeVpcBins(data = data, split = c(1, 2, 1, 2)))
  expect_warning(computeVpcBins(data = data, type = "toto"))
  expect_error(computeVpcBins(data = data, binsSettings = list(a = 1, b = 2)))
  expect_warning(computeVpcBins(data = data, is.binsName = "toto"))
  expect_warning(computeVpcBins(data = data, is.binsName = 3))
  expect_warning(computeVpcBins(data = data, dataName = 3))
  expect_warning(computeVpcBins(data = data, dataName = FALSE))
})

test_that("computeVpcBins returns continuous bins by split when type = continuous", {
  data <- c(1, 2, 3, 1, 1, 2, 5, 5, 3)
  split <- c(1, 2, 1, 2, 1, 2, 2, 2, 1)
  mockery::stub(
    computeVpcBins,
    ".getBins",
    function(data, settings = NULL, bycat = FALSE, categories = c()) {
      if (bycat) {
        return(data.frame(type = "categorical"))
      } else {
        return(data.frame(type = "continuous"))
      }
    }
  )
  expectedBins <-  data.frame(type = c("continuous", "continuous"), split = c(1, 2))
  expect_equal(computeVpcBins(data = data, split = split), expectedBins)

})
