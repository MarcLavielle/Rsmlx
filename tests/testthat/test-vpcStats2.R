context("Virtual Predictive Check statistics 2 - non regression with Monolix")

skip_on_cran()
skip_if(!dir.exists(test_path("projects/")) | length(list.files(test_path("projects/"))) == 0, message = NULL)

get_project <- function(project) {
  mlx.loadProject(project)
  obsName <- mlx.getObservationInformation()$name[1]
  tmpProject <- file.path(tempdir(), basename(project))
  mlx.saveProject(projectFile = tmpProject)
  monolix_vpc_path <- file.path(
    tools::file_path_sans_ext(tmpProject),
    "ChartsData/VisualPredictiveCheck"
  )
  return(list(
    mlx_vpc_path = monolix_vpc_path,
    test_project = tmpProject,
    obs_name = obsName
  ))
}

test_that("vpcStats returns same results as Monolix Chartdata for projects with censoring", {
  projects <- list.files(test_path("projects"), pattern = "censoring[0-9]*.mlxtran", full.names = TRUE)
  for (project in projects) {
    p <- get_project(project)
    vpcData <- vpcStats(p$test_project, obsName = p$obs_name)
    mlx_df <- .readDataset(paste0(p$mlx_vpc_path, "/", p$obs_name, "_percentiles.txt"))
    for (n in names(mlx_df)) {
      expect_true(n %in% names(vpcData$vpcPercentiles))
      # expect_equal(mlx_df[[n]], vpcData$vpcPercentiles[[n]], tolerance=1e-4)
    }
  }
})

test_that("vpcStats returns same results as Monolix Chartdata for projects with multiple observations", {
  projects <- list.files(test_path("projects/"), pattern = "multipleObs[0-9]*.mlxtran", full.names = TRUE)
  for (project in projects) {
    p <- get_project(project)
    vpcData <- vpcStats(p$test_project, obsName = p$obs_name)
    mlx_df <- .readDataset(paste0(p$mlx_vpc_path, "/", p$obs_name, "_percentiles.txt"))
    for (n in names(mlx_df)) {
      expect_true(n %in% names(vpcData$vpcPercentiles))
      expect_equal(mlx_df[[n]], vpcData$vpcPercentiles[[n]], tolerance=1e-4)
    }
  }
})

test_that("vpcStats returns same results as Monolix Chartdata for projects with continuous data", {
  projects <- list.files(test_path("projects/"), pattern = "continuous[0-9]*.mlxtran", full.names = TRUE)
  for (project in projects) {
    p <- get_project(project)
    vpcData <- vpcStats(p$test_project, obsName = p$obs_name)
    mlx_df <- .readDataset(paste0(p$mlx_vpc_path, "/", p$obs_name, "_percentiles.txt"))
    for (n in names(mlx_df)) {
      expect_true(n %in% names(vpcData$vpcPercentiles))
      expect_equal(mlx_df[[n]], vpcData$vpcPercentiles[[n]], tolerance=1e-4)
    }
  }
})

test_that("vpcStats returns same results as Monolix Chartdata for projects with discrete categorical data", {
  projects <- list.files(test_path("projects/"), pattern = "categorical[0-9]*.mlxtran", full.names = TRUE)
  for (project in projects) {
    p <- get_project(project)
    vpcData <- vpcStats(p$test_project, obsName = p$obs_name)
    mlx_df <- .readDataset(paste0(p$mlx_vpc_path, "/", p$obs_name, "_percentiles.txt"))
    for (n in names(mlx_df)) {
      expect_true(n %in% names(vpcData$vpcPercentiles))
      expect_equal(mlx_df[[n]], vpcData$vpcPercentiles[[n]], tolerance=1e-4)
    }
  }
})

test_that("vpcStats returns same results as Monolix Chartdata for projects with discrete count data", {
  projects <- list.files(test_path("projects/"), pattern = "count[0-9]*.mlxtran", full.names = TRUE)
  for (project in projects) {
    p <- get_project(project)
    vpcData <- vpcStats(p$test_project, obsName = p$obs_name)
    mlx_df <- .readDataset(paste0(p$mlx_vpc_path, "/", p$obs_name, "_percentiles.txt"))
    for (n in names(mlx_df)) {
      expect_true(n %in% names(vpcData$vpcPercentiles))
      expect_equal(mlx_df[[n]], vpcData$vpcPercentiles[[n]], tolerance=1e-4)
    }
  }
})

test_that("vpcStats returns same results as Monolix Chartdata for projects with tte data", {
  projects <- list.files(test_path("projects/"), pattern = "tte[0-9]*.mlxtran", full.names = TRUE)
  for (project in projects[1:2]) {
    print(project)
    p <- get_project(project)
    vpcData <- vpcStats(p$test_project, obsName = p$obs_name)
    mlx_df_emp <- .readDataset(paste0(p$mlx_vpc_path, "/", p$obs_name, "_empiricalCurves.txt"))
    mlx_df_theo <- .readDataset(paste0(p$mlx_vpc_path, "/", p$obs_name, "_theoreticalCurves.txt"))
    mlx_df_theo <- .renameColumns(
      mlx_df_theo,
      c("survivalFunction_p5", "survivalFunction_p95", "averageEventNumber_p5", "averageEventNumber_p95"),
      c("survivalFunction_piLower", "survivalFunction_piUpper", "averageEventNumber_piLower", "averageEventNumber_piUpper")
    )
    for (n in names(mlx_df_emp)) {
      expect_true(n %in% names(vpcData$vpcPercentiles))
      expect_equal(mlx_df_emp[[n]], vpcData$vpcPercentiles[[n]], tolerance=1e-4)
    }
    for (n in names(mlx_df_theo)) {
      expect_true(n %in% names(vpcData$vpcPercentiles))
      expect_equal(mlx_df_theo[[n]], vpcData$vpcPercentiles[[n]], tolerance=1e-4)
    }
  }
})

test_that("vpcStats returns data relative to time dose when time = timeSinceLastDose", {
  projects <- list.files(test_path("projects/"), pattern = "timesincelastdose[0-9]*.mlxtran", full.names = TRUE)
  for (project in projects) {
    p <- get_project(project)
    vpcData <- vpcStats(p$test_project, time = "timeSinceLastDose")
    expect_equal(vpcData$timeName, "time_rel")
    expect_true("time_rel" %in% names(vpcData$observations))
    expect_true(max(vpcData$observations$time_rel) <= max(vpcData$observations$time))
  }
})

test_that("vpcStats compute vpc for the correct observation", {
  projects <- list.files(test_path("projects/"), pattern = "multipleObs[0-9]*.mlxtran", full.names = TRUE)
  for (project in projects) {
    p <- get_project(project)
    vpcData <- vpcStats(p$test_project, obsName = "y1")
    expect_equal(vpcData$obsName, "y1")
    vpcData <- vpcStats(p$test_project, obsName = "y2")
    expect_equal(vpcData$obsName, "y2")
  }
})

test_that("vpcStats have the same results when data without censoring no matter censoring arguments", {
  projects <- list.files(test_path("projects/"), pattern = "nocensoring[0-9]*.mlxtran", full.names = TRUE)
  for (project in projects) {
    p <- get_project(project)
    vpcData1 <- vpcStats(p$test_project, useCensored = TRUE)
    vpcData2 <- vpcStats(p$test_project, useCensored = FALSE)
    for (n in names(vpcData1)) {
      expect_equal(vpcData1$vpcPercentiles[[n]], vpcData2$vpcPercentiles[[n]])
    }
  }
})

test_that("vpcStats returns less bins for projects with censored data when useCensored = FALSE", {
  project <- file.path(test_path("projects/"), "censoring1.mlxtran")
  if (file.exists(project)) {
    p <- get_project(project)
    vpcData1 <- vpcStats(p$test_project, useCensored = TRUE)
    vpcData2 <- vpcStats(p$test_project, useCensored = FALSE)
    expect_true(length(vpcData1$vpcPercentiles$bins) > length(vpcData2$vpcPercentiles$bins))
  }
})

test_that("vpcStats handle choice between simulated / BLQ when useCensored = TRUE", {
  projects <- list.files(test_path("projects/"), pattern = "^censoring[0-9]*.mlxtran", full.names = TRUE)
  for (project in projects) {
    print(project)
    p <- get_project(project)
    vpcData <- vpcStats(p$test_project, censoring = "simulated", obsName = p$obs_name)
    expect_equal(vpcData$obsName, paste0(p$obs_name, "_simBlq"))
    expect_true(paste0(p$obs_name, "_simBlq") %in% names(vpcData$observations))
    vpcData <- vpcStats(p$test_project, censoring = "blq", obsName = p$obs_name)
    expect_equal(vpcData$obsName, p$obs_name)
  }
})

# Stratification
test_that("vpcStats returns same results as Monolix Chartdata for projects with covariates", {
  project <- file.path(test_path("projects/"), "covariates.mlxtran")
  if (file.exists(project)) {
    p <- get_project(project)

    vpcData <- vpcStats(p$test_project, obsName = p$obs_name)
    obs <- vpcData$observations
    perc <- vpcData$vpcPercentiles
    expect_length(unique(obs$split), 1)
    expect_equal(unique(obs$split), "All")
    expect_length(unique(perc$split), 1)
    expect_equal(unique(perc$split), "All")

    vpcData <- vpcStats(p$test_project, stratSplit = c("SEQ", "AGE"))
    obs <- vpcData$observations
    perc <- vpcData$vpcPercentiles
    expect_length(unique(obs$split), 4)
    expect_length(unique(perc$split), 4)
    expect_equal(lengths(regmatches(unique(obs$split), gregexpr("#", unique(obs$split)))), rep(2, 4))

    vpcData <- vpcStats(p$test_project, stratSplit = c("SEQ", "WT"), stratScale = list("WT" = 70))
    obs <- vpcData$observations
    perc <- vpcData$vpcPercentiles
    expect_length(unique(obs$split), 4)
    expect_length(unique(perc$split), 4)
    expect_equal(lengths(regmatches(unique(obs$split), gregexpr("#", unique(obs$split)))), rep(2, 4))

    vpcData <- vpcStats(p$test_project, stratSplit = c("SEQ", "AGE"), stratScale = list("AGE" = 150))
    obs <- vpcData$observations
    perc <- vpcData$vpcPercentiles
    expect_length(unique(obs$split), 2)
    expect_length(unique(perc$split), 2)
    expect_equal(lengths(regmatches(unique(obs$split), gregexpr("#", unique(obs$split)))), rep(2, 2))

    vpcData <- vpcStats(p$test_project, stratSplit = c("SEQ"), stratFilter = list("WT" = 1))
    obs <- vpcData$observations
    perc <- vpcData$vpcPercentiles
    expect_length(unique(obs$split), 2)
    expect_length(unique(perc$split), 2)
    expect_equal(lengths(regmatches(unique(obs$split), gregexpr("#", unique(obs$split)))), rep(2, 2))
  }
})
