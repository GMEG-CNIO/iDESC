% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zp_prediction.R
\name{zp_prediction}
\alias{zp_prediction}
\title{zp_prediction: Predict dropout rate from LOESS curve}
\usage{
zp_prediction(mat, norm_factor, span, loess_control)
}
\arguments{
\item{mat}{Count matrix}

\item{norm_factor}{normalization factors from function: "normalization_factors"}

\item{span}{smoothing parameter}

\item{loess_control}{Optional. If set to "control", applies \code{loess.control(surface = "direct")}
to prevent memory issues with large datasets. Use when you encounter warnings like
"k-d tree limited by memory".}
}
\value{
A vector of LOESS predicted dropout rate
}
\description{
zp_prediction: Predict dropout rate from LOESS curve
}
