#' msaeHB : Multivariate Small Area Estimation using Hierarchical Bayesian Method
#'
#' Implements area level of multivariate small area estimation using hierarchical Bayesian (HB) method under Normal and T distribution. The 'rjags' package is employed to obtain parameter estimates. For the reference, see Rao and Molina (2015) <doi:10.1002/9781118735855>.
#'
#' @section Author(s):
#' Azka Ubaidillah \email{azka@@stis.ac.id} and Novia Permatasari  \email{novia.permatasari@@bps.go.id}
#'
#' \strong{Maintainer}: Novia Permatasari  \email{novia.permatasari@@bps.go.id}
#'
#'
#' @section Functions:
#' \describe{
#'   \item{\code{mHBNormal}}{Estimate multivariate small area estimation under normal distribution}
#'   \item{\code{mHBT}}{Estimate multivariate small area estimation under normal distribution}
#'}
#' @section Reference:
#' \itemize{
#'   \item{Rao, J.N.K & Molina. (2015). Small Area Estimation 2nd Edition. New York: John Wiley and Sons, Inc. <doi:10.1002/9781118735855>.}
#' }
#'
#'
#' @docType package
#' @name msaeHB
#'
#' @import rjags
#' @import coda
#' @import stats
#' @import grDevices
#' @import graphics

NULL
