#' @export
RunMH_my <- function (Target, init, B, concentration, h, dat, pars) 
{
  zz = proc.time()
  p <- length(init)
  Y <- array(0, c(B, p))
  moveCount <- rep(0, p)
  init <- init/sum(init)
  a <- concentration
  ycurrent <- SALTSampler::Logit(init)
  for (i in 1:B) {
    for (j in sample(p)) {
      ycand <- SALTSampler::PropStep(ycurrent, j, h[j])
      if (any(is.na(ycand) | is.infinite(ycand) | is.nan(ycand))) {
        move <- FALSE
      }
      else {
        move <- (log(runif(1)) < attr(ycand, "dbt") + 
                   Target(ycand, ycurrent, a, dat, pars))
      }
      if (!is.na(move) & move) {
        ycurrent <- ycand
        moveCount[j] <- moveCount[j] + 1
      }
    }
    Y[i, ] <- ycurrent
  }
  runTime <- proc.time() - zz
  S <- matrix(exp(SALTSampler::LogPq(Y)$logp), nrow(Y))
  return(list(Y = Y, S = S, runTime = runTime, moveCount = moveCount, 
              p = p, init = init, B = B, concentration = concentration, 
              h = h, dat = dat))
}
