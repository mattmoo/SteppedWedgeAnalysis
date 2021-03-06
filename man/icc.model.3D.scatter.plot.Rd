% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.helpers.R
\name{icc.model.3D.scatter.plot}
\alias{icc.model.3D.scatter.plot}
\title{Creates a 3D scatter plot of ICC vs two other variables, ICC on y-axis.}
\usage{
icc.model.3D.scatter.plot(
  icc.model,
  vars.to.plot = NULL,
  x.vals = seq(0, 10, by = 0.1),
  y.vals = seq(0, 10, by = 0.1)
)
}
\arguments{
\item{icc.model}{A model for predicting ICC from given variables.}

\item{vars.to.plot}{Vector of two string names of the variable that will be
on the x- and y-axes. Will autodetect if NULL, but that will fail if there
are more than three potential variables to plot ICC against. (Default:
NULL)}

\item{x.vals}{x values to generate values for (Default: 1:10 by 0.1)}

\item{y.vals}{y values to generate values for (Default: 1:10 by 0.1)}
}
\value{
A plotly 3D scatter plot of predicted values from a model.
}
\description{
Creates a 3D scatter plot of ICC vs two other variables, ICC on y-axis.
}
