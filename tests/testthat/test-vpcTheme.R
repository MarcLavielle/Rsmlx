context("Virtual Predictive Check Theme")

test_that("createVpcTheme returns an error when inapropriate arguments", {
  expect_error(createVpcTheme(update = list(obs_color = 5)))
  expect_error(createVpcTheme(update = list(obs_shape = 12.5)))
  expect_error(createVpcTheme(update = list(emp_perc_legend = -1)))
  expect_error(createVpcTheme(update = list(theo_pi_median_color = 4)))
})

test_that("createVpcTheme returns a warning when wrong element", {
  expect_warning(createVpcTheme(update = list(toto = 5)))
})

test_that("createVpcTheme returns a vpc_theme object which is a list", {
  expect_s3_class(createVpcTheme(), "vpc_theme")
  expect_type(createVpcTheme(), "list")
})

test_that("createVpcTheme replace the value in list by the input argument", {
  expect_equal(
    createVpcTheme(update = list(obs_color = rgb(80, 130, 180, maxColorValue = 255)))$obs_color,
    rgb(80, 130, 180, maxColorValue = 255)
  )
  expect_equal(createVpcTheme(update = list(emp_perc_size = 5))$emp_perc_size, 5)
  expect_equal(createVpcTheme(update = list(theo_perc_legend = "toto"))$theo_perc_legend, "toto")
  expect_equal(createVpcTheme(update = list(outlier_areas_legend = NA))$outlier_areas_legend, NA)
  expect_equal(createVpcTheme(update = list(theo_pi_median_fill = NA))$theo_pi_median_fill, NA)
})
