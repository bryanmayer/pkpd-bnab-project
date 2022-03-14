
#' Calculate neutralization given experimental titer
#' @param titer positive value, ratio of concentration of IC50
#' @return scalar or numerical vector
titer2neut = function(titer) {1 - 1/(1+titer)}

#' Calculate experimental titer given neutralization
#' @param neut neutralization between 0 and 1
#' @return scalar or numerical vector
neut2titer = function(neut) {neut/(1 - neut)}

#' Calculate IIP given neutralization
#' @param neut neutralization between 0 and 1
#' @return scalar or numerical vector
neut2iip = function(neut) pmin(-log10(1 - neut), -log10(.Machine$double.eps))


#' Calculate BH neutralization given two titers
#'
#' @param ID50_1 titer bnab1
#' @param ID50_2 titer bnab1
#'
#' @return scalar or numerical vector
#'
titer2BH = function(ID50_1, ID50_2) {
  1 - (1 - titer2neut(ID50_1)) * (1 - titer2neut(ID50_2))
}

#' Calculate BH experimental titer given two titers
#'
#' @param ID50a titer bnab1
#' @param ID50b titer bnab2
#' @param NaN_as_NA flag to ignore simulation settings that created boundary conditions for NaNs
#'
#' @return scalar or numerical vector
#'
calc_BHID50 = function(ID50a, ID50b, NaN_as_NA = F){
  zeros = which(ID50a < sqrt(.Machine$double.eps) & ID50b < sqrt(.Machine$double.eps))
  ID50a = pmax(ID50a, sqrt(.Machine$double.eps))
  ID50b = pmax(ID50b, sqrt(.Machine$double.eps))

  x = (-(ID50a + ID50b) + sqrt((ID50a + ID50b)^2 + 4*ID50b*ID50a))/(2*ID50a*ID50b)
  if(any(is.nan(x))) {
    if(!NaN_as_NA) stop(paste("BH ID50 nAn, input:", ID50a, ID50b))
    x[is.nan(x)] <- NA_real_
  }
  out = 1/x

  # this handles some cluster issues with high/float values
  infinites = which(is.infinite(out))

  if(length(infinites > 0)) out[infinites] = pmax(ID50a, ID50b)[infinites]
  if(length(zeros > 0)) out[zeros] = 0

  out

}

#' A convoluted way to do a labeller function, would remove if code was refactored
#' @param interaction_var a variable vector that looks like levels below
#'
interaction_labeller = function(interaction_var){
  factor(interaction_var,
         labels = c("Min. Neut.", "Max. Neut.", "Additivity", "Bliss-Hill"),
         levels = c("minNeut", "maxNeut", "additivity", "BH")
  )
}

