#' Selected Dataset of HAM10000
#'
#' A dataset containing 1000 pigmented skin lesion 64*64 (RGB) images from
#' HAM10000 as predictors and corresponding diagnostic category as responses.
#'
#' @format A list with the first item as predictors, the second as responses.
#' \describe{
#'     \item{hamimgs}{predictors, a matrix with each row representing an image}
#'     \item{hamlab}{responses, a matrix with each row a dummy variable
#'     representing categories.}
#'     ...
#' }
#' @source \url{https://www.kaggle.com/kmader/skin-cancer-mnist-ham10000}
#'
"HAM"
