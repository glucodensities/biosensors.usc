#' biosensors.usc Package
#'
#' Biosensor data have the potential to improve disease control and detection. However,
#' the analysis of these data under free-living conditions is not feasible with current
#' statistical techniques. To address this challenge, we introduce a new functional
#' representation of biosensor data, termed the glucodensity, together with a data
#' analysis framework based on distances between them. The new data analysis procedure
#' is illustrated through an application in diabetes with continuous-time glucose
#' monitoring (CGM) data. In this domain, we show marked improvement with respect to
#' state-of-the-art analysis methods. In particular, our findings demonstrate that (i)
#' the glucodensity possesses an extraordinary clinical sensitivity to capture the
#' typical biomarkers used in the standard clinical practice in diabetes; (ii) previous
#' biomarkers cannot accurately predict glucodensity, so that the latter is a richer
#' source of information and; (iii) the glucodensity is a natural generalization of the
#' time in range metric, this being the gold standard in the handling of CGM data.
#' Furthermore, the new method overcomes many of the drawbacks of time in range metrics
#' and provides more indepth insight into assessing glucose metabolism.
#'
#' @docType package
#'
#' @author Juan C. Vidal \email{juan.vidal@usc.es}
#' @author Marcos Matabuena \email{marcos.matabuena@usc.es}
#'
#' @name biosensors.usc
#' @useDynLib biosensors.usc, .registration=TRUE
NULL
# > NULL
