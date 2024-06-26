% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/buildmlx.R
\name{buildmlx}
\alias{buildmlx}
\title{Automatic statistical model building}
\usage{
buildmlx(
  project = NULL,
  final.project = NULL,
  model = "all",
  prior = NULL,
  weight = NULL,
  coef.w1 = 0.5,
  paramToUse = "all",
  covToTest = "all",
  covToTransform = "none",
  center.covariate = FALSE,
  criterion = "BICc",
  linearization = FALSE,
  ll = T,
  test = T,
  direction = NULL,
  steps = 1000,
  n.full = 10,
  max.iter = 20,
  explor.iter = 2,
  fError.min = 0.001,
  seq.cov = FALSE,
  seq.cov.iter = 0,
  seq.corr = TRUE,
  p.max = 0.1,
  p.min = c(0.075, 0.05, 0.1),
  print = TRUE,
  nb.model = 1
)
}
\arguments{
\item{project}{a string: the initial Monolix project}

\item{final.project}{a string: the final Monolix project (default adds "_built" to the original project)}

\item{model}{components of the model to optimize c("residualError", "covariate", "correlation"), (default="all")}

\item{prior}{list of prior probabilities for each component of the model (default=NULL)}

\item{weight}{list of penalty weights for each component of the model (default=NULL)}

\item{coef.w1}{multiplicative weight coefficient used for the first iteration only (default=0.5)}

\item{paramToUse}{list of parameters possibly function of covariates (default="all")}

\item{covToTest}{components of the covariate model that can be modified   (default="all")}

\item{covToTransform}{list of (continuous) covariates to be log-transformed (default="none")}

\item{center.covariate}{TRUE/FALSE center the covariates of the final model (default=FALSE)}

\item{criterion}{penalization criterion to optimize c("AIC", "BIC", "BICc", gamma) (default=BICc)}

\item{linearization}{TRUE/FALSE whether the computation of the likelihood is based on a linearization of the model (default=FALSE)}

\item{ll}{TRUE/FALSE  compute the observe likelihood and the criterion to optimize at each iteration}

\item{test}{TRUE/FALSE  perform additional statistical tests for building the model (default=TRUE)}

\item{direction}{method for covariate search c("full", "both", "backward", "forward"), (default="full" or "both")}

\item{steps}{maximum number of iteration for stepAIC (default=1000)}

\item{n.full}{maximum number of covariates for an exhaustive comparison of all possible covariate models (default=10)}

\item{max.iter}{maximum number of iterations (default=20)}

\item{explor.iter}{number of iterations during the exploratory phase (default=2)}

\item{fError.min}{minimum fraction of residual variance for combined error model (default = 1e-3)}

\item{seq.cov}{TRUE/FALSE whether the covariate model is built before the correlation model}

\item{seq.cov.iter}{number of iterations before building the correlation model (only when seq.cov=F, default=0)}

\item{seq.corr}{TRUE/FALSE whether the correlation model is built iteratively (default=TRUE)}

\item{p.max}{maximum p-value used for removing non significant relationships between covariates and individual parameters (default=0.1)}

\item{p.min}{vector of 3 minimum p-values used for testing the components of a new model (default=c(0.075, 0.05, 0.1))}

\item{print}{TRUE/FALSE display the results (default=TRUE)}

\item{nb.model}{number of models to display at each iteration (default=1)}
}
\value{
a new Monolix project with a new statistical model.
}
\description{
Automatic statistical model building is available directly in the lixoftConnectors package using the function
\code{runModelBuilding}.
}
\details{
buildmlx uses SAMBA (Stochastic Approximation for Model Building Algorithm), an iterative procedure to accelerate and optimize the process of model building by identifying 
at each step how best to improve some of the model components. This method allows to find 
the optimal statistical model which minimizes some information criterion in very few steps.

Penalization criterion can be either a custom penalization of the form gamma*(number of parameters),
AIC (gamma=2) or BIC (gamma=log(N)).

Several strategies can be used for building the covariate model at each iteration of the algorithm: 
\code{direction="full"} means that all the possible models are compared (default when the number of covariates
is less than 10). Othrwise, \code{direction} is the mode of stepwise search of \code{stepAIC {MASS}}, 
can be one of "both", "backward", or "forward", with a default of "both" when there are at least 10 covariates.
See https://monolix.lixoft.com/rsmlx/ for more details.
}
\examples{
\dontrun{
# RsmlxDemo1.mlxtran is a Monolix project for modelling the pharmacokinetics (PK) of warfarin 
# using a PK model with parameters ka, V, Cl.

# By default, buildmlx will compute the best statistical model in term of BIC, i.e , 
# the best covariate model, the best correlation model for the three random effects and the best 
# residual error model in terms of BIC. 
# In this example, three covariates (wt, age, sex) are available with the data and will be used 
# for building the covariate model for the three PK parameters:
r1 <- buildmlx(project="RsmlxDemo1.mlxtran")
  
# Here, the covariate model will be built for V and Cl only and log-transformation of all 
# continuous covariates will also be considered:
r2 <- buildmlx(project="RsmlxDemo1.mlxtran", paramToUse=c("V", "Cl"), covToTransform="all") 

# Only the covariate model will be  built, using AIC instead of BIC:
r3 <- buildmlx(project="RsmlxDemo1.mlxtran", model="covariate", criterion="AIC") 
}

# See http://monolix.lixoft.com/rsmlx/buildmlx/ for detailed examples of use of buildmlx
# Download the demo examples here: http://monolix.lixoft.com/rsmlx/installation

}
\seealso{
\code{getModelBuildingSettings} settings for model building with lixoftConnectors \cr
\code{runModelBuilding} run model building with lixoftConnectors \cr
\code{getModelBuildingResults} results for model building with lixoftConnectors
}
