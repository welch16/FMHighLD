% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/init.R
\name{init_iteration}
\alias{init_iteration}
\title{Performs the initial iteration of the FM-HighLD model}
\usage{
init_iteration(formula, data, multi_trait = FALSE)
}
\arguments{
\item{formula}{formula used for fitting the underlying linear model used by
FM-HighLD}

\item{data}{`data.frame` used for the FM-HighLD model}

\item{multi_trait}{a logical indicator determining if the model is for
multi_trait or single-trait fine-mapping}
}
\value{
a list with the models, and probability matrix for the EM-algorithm
}
\description{
Performs the initial iteration of the FM-HighLD model
}
