context("Virtual Predictive Check Stats and Plot")

skip_on_cran()
skip_if_not_installed("lixoftConnectors")

initRsmlx(warnings = FALSE, info = FALSE)
demo_path <- file.path(path.expand("~"), "lixoft", "monolix", paste0("monolix", mlx.getLixoftConnectorsState()$version), "demos")
skip_if(!dir.exists(demo_path), message = NULL)

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
  joint_demos <- file.path(demo_path, "4.joint_models/4.1.continuous_PKPD")
  projects <- list.files(path = joint_demos , pattern = '[.]mlxtran$', full.names = T, include.dirs = F, recursive = T)
  for (project in projects) {
    p <- get_project(project)
    res <- vpc(p$test_project, obsName = p$obs_name, plot = FALSE)
    expect_type(res, "list")
    expect_length(res, 4)
    expect_true(all(names(res) %in% c("vpcPercentiles", "observations", "obsName", "timeName")))
    expect_s3_class(res$observations, "data.frame")
    expect_s3_class(res$vpcPercentiles, "data.frame")
  }
})

test_that("vpc returns a ggplot object when plot = TRUE", {
  joint_demos <- file.path(demo_path, "4.joint_models/4.1.continuous_PKPD")
  projects <- list.files(path = joint_demos , pattern = '[.]mlxtran$', full.names = T, include.dirs = F, recursive = T)
  for (project in projects) {
    p <- get_project(project)
    res <- vpc(p$test_project, obsName = p$obs_name, plot = TRUE)
    expect_s3_class(res, "ggplot")
  }
})

test_that("vpc returns a ggplot object with one subplot for projects with tte data when meanNumberEventsCurve not displayed", {
  tte_demos <- file.path(demo_path, "3.models_for_noncontinuous_outcomes/3.3.time_to_event_data_model")
  projects <- list.files(path = tte_demos , pattern = '[.]mlxtran$', full.names = T, include.dirs = F, recursive = T)
  for (project in projects[1:2]) {
    p <- get_project(project)
    res <- vpc(p$test_project, obsName = p$obs_name, plot = TRUE, eventDisplay = c("survivalCurve", "empiricalCurve", "predictionInterval"))
    expect_s3_class(res, "ggplot")
  }
})

test_that("vpc returns a ggplot object with one subplot for projects with tte data when survivalCurve not displayed", {
  tte_demos <- file.path(demo_path, "3.models_for_noncontinuous_outcomes/3.3.time_to_event_data_model")
  projects <- list.files(path = tte_demos , pattern = '[.]mlxtran$', full.names = T, include.dirs = F, recursive = T)
  for (project in projects[1:2]) {
    p <- get_project(project)
    res <- vpc(p$test_project, obsName = p$obs_name, plot = TRUE, eventDisplay = c("meanNumberEventsCurve", "empiricalCurve", "predictionInterval"))
    expect_s3_class(res, "ggplot")
  }
})

test_that("vpc returns a ggplot object with two subplots for projects with tte data when survivalCurve and meanNumberEventsCurve are displayed", {
  tte_demos <- file.path(demo_path, "3.models_for_noncontinuous_outcomes/3.3.time_to_event_data_model")
  projects <- list.files(path = tte_demos , pattern = '[.]mlxtran$', full.names = T, include.dirs = F, recursive = T)
  for (project in projects[1:2]) {
    p <- get_project(project)
    res <- vpc(p$test_project, obsName = p$obs_name, plot = TRUE, eventDisplay = c("survivalCurve", "meanNumberEventsCurve", "empiricalCurve", "predictionInterval"))
    expect_s3_class(res, "gtable")
  }
})
