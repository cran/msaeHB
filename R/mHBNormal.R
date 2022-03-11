#' @title Multivariate Small Area Estimation using Hierarchical Bayesian under Normal Distribution
#' @description This function implements small area estimation using hierarchical bayesian to variable of interest that assumed to be a multivariate normal distribution.
#' @param formula an object of class list of formula, describe the model to be fitted
#' @param vardir vector containing name of sampling variances of direct estimators  in the following order : \code{var1, var2, . , var(k) , cov12, . cov1k, cov23, . , cov(k-1)(k)}
#' @param iter.update number of updates with default \code{3}
#' @param iter.mcmc number of total iterations per chain with default \code{10000}
#' @param thin thinning rate, must be a positive integer with default \code{2}
#' @param burn.in number of iterations to discard at the beginning with default \code{2000}
#' @param data dataframe containing the variables named in \code{formula} and \code{vardir}
#' @return The function returns a list with the following objects:
#' \describe{
#'    \item{Est}{A vector with the values of Small Area mean Estimates using Hierarchical bayesian method }
#'    \item{coefficient}{A dataframe with the estimated model coefficient}
#'    \item{plot}{Trace, Density, Autocorrelation Function Plot of MCMC samples}
#' }
#' @examples
#'   ## Load dataset
#'   data(datasaeNorm)
#'   ## Using parameter 'data'
#'   Fo <- list(f1=Y1~X1+X2,
#'              f2=Y2~X1+X2)
#'   vardir <- c("v1", "v2", "v12")
#'   m1 <- mHBNormal(formula=Fo, vardir=vardir,
#'   iter.update = 1, iter.mcmc = 1000,
#'   thin = 2, burn.in = 200, data=datasaeNorm)
#'
#' @export mHBNormal
mHBNormal <- function(formula, vardir, iter.update=3,
                      iter.mcmc=10000, thin = 2,
                      burn.in =2000, data){


  result <- list(Est = NA, refVar = NA, coefficient = NA,
                 plot=NA)

  Y.var = sapply(formula, function(x){all.vars(x)[1]})
  Y= data[,Y.var]
  m=nrow(Y)
  P=ncol(data.frame(Y))

  if(P==1){
    stop("Use Normal() on saeHB.")
  }
  if (!any(is.na(Y)))  {

    X=lapply(formula, function(x){model.matrix(x,data)})
    S.E = df2R(colMeans(data[,vardir]),P)

    nvar <- sapply(X, ncol)
    nvar_tot = sum(nvar)

    tau.ua=tau.ub=tau.u=rep(1,P)
    m.u = rep(0,P)
    mu.b = rep(0,nvar_tot)
    t.b = rep(1,nvar_tot)

    R=diag(1,P)
    k.a=v.a=P
    k.b=v.b=1

    Iter=iter.update

    for (i in 1:iter.update) {
      lt <- list()
      for(i in c(1:P)){
        if(i == 1){lt[["X"]] = X[[i]]} else {lt[["X"]] = cbind(lt[["X"]], X[[i]])}
        lt[["mu.b"]] = mu.b
        lt[["t.b"]] = t.b
        lt[["nvar"]] = nvar
      }
      dat <- list("m" = m, "P"=P, "Y" = Y, "R"=R, "V"=S.E/(v.a/v.b),
                  "k.a"=k.a,"k.b"=k.b,"v.a"=v.a,"v.b"=v.b,
                  "m.u"=m.u)
      dat <- append(dat, lt)
      inits <- list(b = mu.b)

      cat("model {
					for (i in 1:m) {
							Y[i,1:P] ~ dmnorm(MU[i,1:P],tau.y[1:P, 1:P])
							MU[i,1:P]<- sum(b*X[i,])+U[i,1:P]
							U[i,1:P]~dmnorm(m.u,tau.u)
					}

					for(q in 1:sum(nvar)){
					  b[q] ~ dnorm(mu.b[q],t.b[q])
					}

					tau.u[1:P,1:P] ~ dwish(R[1:P,1:P],k)
					tau.y[1:P,1:P] ~ dwish(V[1:P,1:P],v)
					k ~ dgamma(k.a,k.b)
					v ~ dgamma(v.a,v.b)

			}", file="saeHBnormal.txt")

      jags.m <- jags.model( file = "saeHBnormal.txt", data=dat, inits=inits, n.chains=1, n.adapt=500 )
      file.remove("saeHBnormal.txt")
      params <- c("MU","b", "tau.u","k","v")
      samps <- coda.samples( jags.m, params, n.iter=iter.mcmc, thin=thin)
      samps1 <- window(samps, start=burn.in+1, end=iter.mcmc)

      result_samps=summary(samps1)
      beta=result_samps$statistics[(m*P+1):(m*P+sum(nvar)),1:2]
      for (i in 1:nvar_tot){
        mu.b[i]  = beta[i,1]
        t.b[i] = 1/(beta[i,2]^2)
      }

      k.a = result_samps$statistics[m*P+nvar_tot+1,1]^2/result_samps$statistics[m*P+nvar_tot+1,2]^2
      k.b = result_samps$statistics[m*P+nvar_tot+1,1]/result_samps$statistics[m*P+nvar_tot+1,2]^2
      v.a = result_samps$statistics[nrow(result_samps$statistics),1]^2/result_samps$statistics[nrow(result_samps$statistics),2]^2
      v.b = result_samps$statistics[nrow(result_samps$statistics),1]/result_samps$statistics[nrow(result_samps$statistics),2]^2
      Tau.u=matrix(c(result_samps$statistics[(m*P+nvar_tot+2):(nrow(result_samps$statistics)-1),1]),P,P)
      R=result_samps$statistics[m*P+nvar_tot+1,1]*solve(Tau.u)

    }

    result_samps=summary(samps)
    b.varnames <- list()
    result_mcmc <- samps[,(m*P+1):(m*P+sum(nvar))]
    beta=result_samps$statistics[(m*P+1):(m*P+sum(nvar)),1:2]
    mu=result_samps$statistics[1:(m*P),1:2]

    Estimation=data.frame(mu)

    Quantiles <- as.data.frame(result_samps$quantiles[1:(m*P+nvar_tot),])
    q_mu <- Quantiles[1:(m*P),]
    q_beta <- (Quantiles[(m*P+1):(m*P+nvar_tot),])
    beta <- cbind(beta,q_beta)
    Estimation <- data.frame(Estimation,q_mu)
    colnames(Estimation) <- c("MEAN","SD","2.5%","25%","50%","75%","97.5%")

  } else {

    stop("Data contains NA values")

  }

  result$Est         = Estimation
  result$coefficient = beta
  result$plot        = list(graphics.off() ,par(mar=c(2,2,2,2)),autocorr.plot(result_mcmc,col="brown2",lwd=2),plot(result_mcmc,col="brown2",lwd=2))
  return(result)

}
