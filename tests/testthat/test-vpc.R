context("Virtual Predictive Check Stats and Plot")

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

test_that("vpc returns a dataframe when plot = FALSE", {
  project <- file.path(test_path("projects"), "continuous1.mlxtran")
  if (file.exists(project)) {
    p <- get_project(project)
    res <- vpc(p$test_project, plot = FALSE)
    expect_type(res, "list")
    expect_length(res, 4)
    expect_true(all(names(res) %in% c("vpcPercentiles", "observations", "obsName", "timeName")))
    expect_s3_class(res$observations, "data.frame")
    expect_s3_class(res$vpcPercentiles, "data.frame")
  }
})

test_that("vpc returns a ggplot object when plot = TRUE", {
  project <- file.path(test_path("projects"), "continuous1.mlxtran")
  if (file.exists(project)) {
    p <- get_project(project)
    res <- vpc(p$test_project, plot = TRUE)
    expect_s3_class(res, "ggplot")
  }
})

test_that("vpc returns a ggplot object with one subplot for projects with tte data when meanNumberEventsCurve not displayed", {
  project <- file.path(test_path("projects"), "tte_cut.mlxtran")
  if (file.exists(project)) {
    p <- get_project(project)
    res <- vpc(p$test_project, plot = TRUE, eventDisplay = c("survivalCurve", "empiricalCurve", "predictionInterval"))
    expect_s3_class(res, "ggplot")
  }
})

test_that("vpc returns a ggplot object with one subplot for projects with tte data when survivalCurve not displayed", {
  project <- file.path(test_path("projects"), "tte_cut.mlxtran")
  if (file.exists(project)) {
    p <- get_project(project)
    res <- vpc(p$test_project, plot = TRUE, eventDisplay = c("meanNumberEventsCurve", "empiricalCurve", "predictionInterval"))
    expect_s3_class(res, "ggplot")
  }
})

test_that("vpc returns a ggplot object with two subplots for projects with tte data when survivalCurve and meanNumberEventsCurve are displayed", {
  project <- file.path(test_path("projects"), "tte_cut.mlxtran")
  if (file.exists(project)) {
    p <- get_project(project)
    res <- vpc(p$test_project, plot = TRUE, eventDisplay = c("survivalCurve", "meanNumberEventsCurve", "empiricalCurve", "predictionInterval"))
    expect_s3_class(res, "gtable")
  }
})
