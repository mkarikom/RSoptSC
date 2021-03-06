% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PlottingFunctions.R
\name{ColorHue}
\alias{ColorHue}
\title{Get a vector of n equally spaced rgb colors}
\usage{
ColorHue(n, starthue = 15, endhue = 360, luminance = 65, chroma = 100)
}
\arguments{
\item{n}{integer number of hex codes to return}

\item{starthue}{real hue argument for the grDevices::hcl() generates value 1}

\item{endhue}{real hue argument}

\item{luminance}{the luminance of the hcl value}

\item{chroma}{the chroma of the hcl value}
}
\description{
Get a vector of n equally spaced rgb colors.
}
