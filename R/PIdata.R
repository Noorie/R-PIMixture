#' PIdata includes two datasets, PIdata1 and PIdata2, which are simulation data.
#'
#' PIdata1 and PIdata2 are from simple random samples and stratified simple random samples, respectively. The two variables of "samp.weight" and "strata" should be included on the dataset of stratified random samples. The variable of "strata.frac" is needed for the design-based variance calculation with the option "sample.design=2".
#'
#' @docType data
#'
#' @usage data(PIdata)
#'
#'
#' @keywords datasets
#'
#' @format Objects of class data.frame, and missing values are coded as -999.
#'
#' @examples
#' data(PIdata)
#' head(PIdata1)
#' head(PIdata2)
"PIdata1"
