#' compdag
#'
#' @param x data of the first community
#' @param y data of the second scommunity
#' @param paramx the fixed mu for the first community
#' @param paramy the fixed mu for the second community
#'
#' @return cd=1 indicates that the community-level causal direction is x to y;
#'         cd=-1 indicates that the community-level causal direction is y to x;
#'         cd = 0 indicate that there is no causal relationship between x and y.
#'         mEst is the estimated microbe-level causal effects.
#'         signalEst is the estimated community-level causal effect.
#'         x2y is the looic of the model x to y;
#'         y2x is the looic of the model y to x;
#'         x0y is the looic of the model that there is no casual relationship.
#' @export


compdag <- function(x,y,paramx,paramy){
  px = nrow(x)
  py = nrow(y)
  n = ncol(x)
  #fit baseline
  datax <- list(n=n,p1=px,y1=x)
  fitx <- rstan::sampling(stanmodels$compDAG0,datax,iter=5000,chains=1,seed=1,thin = 20)
  paramsx <- rstan::extract(fitx)
  datay <- list(n=n,p1=py,y1=y)
  fity <- rstan::sampling(stanmodels$compDAG0,datay,iter=5000,chains=1,seed=1,thin = 20)
  paramsy <- rstan::extract(fity)

  #fit regression
  dataxy <- list(n=n,p1=px,p2=py,y1=x,y2=y,mu2=paramx)
  fitxy <- rstan::sampling(stanmodels$compDAG1,dataxy,iter=5000,chains=1,seed=1,thin = 20)
  paramsxy = rstan::extract(fitxy)
  datayx <- list(n=n,p1=py,p2=px,y1=y,y2=x,mu2=paramy)
  fityx <- rstan::sampling(stanmodels$compDAG1,datayx,iter=5000,chains=1,seed=1,thin = 20)
  paramsyx = rstan::extract(fityx)

  #calculate likelihood
  loglx = paramsx$log_lik
  logly = paramsy$log_lik
  logxy = paramsxy$log_lik
  logyx = paramsyx$log_lik

  logx2y = loglx+logxy
  logy2x = logly+logyx
  logx0y = loglx+logly

  x2yloo = rstan::loo(logx2y)$looic
  y2xloo = rstan::loo(logy2x)$looic
  x0yloo = rstan::loo(logx0y)$looic

  # cd=1:causal direction is from x to y, cd=-1:causal direction is from y to x, cd=0:there is no causal relationship between y and x.
  if (x2yloo==min(c(x2yloo,y2xloo,x0yloo)))
    cd = 1
  if (y2xloo==min(c(x2yloo,y2xloo,x0yloo)))
    cd = -1
  if (x0yloo==min(c(x2yloo,y2xloo,x0yloo)))
    cd = 0

  # get estimations
  mEstxy <- apply(paramsxy$mt1,c(2,3),mean)
  signalEstxy <- mean(paramsxy$pi1)
  # mu1Estxy <- mean(paramsxy$a)
  # calculate LOOIC for two directions


  return(list(cd = cd,mEst = mEstxy,signalEst = signalEstxy,x2y=x2yloo,y2x=y2xloo,x0y = x0yloo))
}

