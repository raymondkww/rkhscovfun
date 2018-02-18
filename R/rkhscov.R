################################
##### Sobolev RK function ######
################################

k1 <- function(t){
  return(t-.5)
}
k2 <- function(t){
  return( (k1(t)^2-1/12)/2 )
}
k4 <- function(t){
  return( (k1(t)^4-k1(t)^2/2+7/240)/24 )
}
K.sob <- function(s,t){
  ans <- 1 + k1(s)*k1(t) + k2(s)*k2(t) - k4(abs(s-t))
  return(ans)
}

#############################
#### svec transformation ####
#############################
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

svec <- function(X){
  # symmetric X of size n times n
  ind <- as.vector(row(X)> col(X))
  ind2 <- as.vector(row(X)>= col(X))
  X[ind] <- X[ind] * sqrt(2)
  return(X[ind2])
}

svec.inv <- function(x){
  nn <- length(x)
  n <- (-1 + sqrt(1+4 *2*nn))/2
  X <- matrix(nr=n, nc=n)
  ind <- as.vector(row(X)> col(X))
  ind2 <- as.vector(row(X)>= col(X))
  ind3 <- as.vector(row(X)< col(X))
  X[ind2] <- x
  X[ind] <- X[ind]/sqrt(2)
  X[ind3] <- t(X)[ind3]
  return(X)
}

smat <- function(X){
  # symmetric X of size n^2 times n^2
  n <- sqrt(dim(X)[1])
  if (!is.wholenumber(n)) stop("n is not an integer")
  temp <- matrix(1:(n*n), nr=n, nc=n)
  index <- temp[as.vector(row(temp)> col(temp))]
  index2 <- t(temp)[as.vector(row(temp)> col(temp))]
  ind2 <- as.vector(row(temp)>= col(temp))
  Y <- X
  Y[,index]  <- (Y[,index] + Y[,index2])/sqrt(2)
  Y[index,]  <- (Y[index,] + Y[index2,])/sqrt(2)
  return(Y[ind2, ind2])
}


#####################
#### preparation ####
#####################

# using pivoted cholesky (warning is expected from chol)
getM.chol <- function(tt, tol=1e-6){
  K <- getK(tt)
  temp <- chol(K, pivot=T, tol=tol)
  r <- attr(temp, "rank")
  oo <- order(attr(temp, "pivot"))
  return(t(temp[1:r, oo]))
}

# sampling a subset of time points
getM.sam <- function(tt, nsam){
  K <- getK(tt)
  ind <- sample(length(tt), nsam)
  u <- svd(K[,ind])$u
  a <- chol(t(u) %*% K %*% u)
  return(u %*% t(a))
}

prep2 <- function(time, x, subject, Mmethod="eig", tol=1e-6, nsam=50){

  if (any(time<0)|| any(time>1)){
    stop("time has to be between 0 and 1!")
  }

  # prepare for fitting
  Xs <- list()
  tts <- list()
  ms <- NULL
  i <- 0
  for (zz in unique(subject)){
    if (sum(subject==zz)>1){
      i <- i+1
      tts[[i]] <- time[subject == zz]
      Xs[[i]] <- as.double(x[subject == zz])
      ms <- c(ms, length(tts[[i]]))
    }
  }

  # getting Mi
  n <- length(Xs)
  # note: some element of time may be removed, so use unlist(tts) instead of time
  if (Mmethod=="eig"){
    M <- getM(unlist(tts))
  } else if (Mmethod=="chol"){
    M <- getM.chol(unlist(tts), tol=tol)
  } else if (Mmethod=="sam"){
    M <- getM.sam(unlist(tts), nsam=nsam)
  }

  ii <- c(0, cumsum(ms))
  Ms <- list()
  for (i in (1:n)){
    Ms[[i]] <- M[(ii[i]+1):ii[i+1],]
  }

  obj <- Qlossprep_cpp(Xs, Ms, 1:n)
  return(list(Xs=Xs, Ms=Ms, ms=ms, M=M, tts=tts, R=obj$R, Qv=obj$Qv, c=obj$c))
}

#################
#### rkhscov ####
#################
# TO-DO: checking on the control parameters
rkhscov.alg.control <- function(Mmethod="eig", preptol=1e-6, nsam=50, L=1,
                                eta=2, alpha=0.9, maxit=10000, traceit=FALSE,
                                max.rank=-1, tol=1e-12, variant=2, cond=2){
  return(list(Mmethod=Mmethod, preptol=preptol, nsam=nsam, L=L, eta=eta, alpha=alpha, maxit=maxit, traceit=traceit, max.rank=max.rank, tol=tol, variant=variant, cond=cond))
}


# penalty: lam * (gam * sum_j w_j |eig_j| + (1-gam)/2 * sum_j eig_j^2) if pos=false

rkhscov <- function(time, x, subject, lam, gam=1, weight=NULL, centered=FALSE, B0v=NULL, pos=TRUE, control=list())
{
  if ((min(time) < 0) || (max(time) > 1)){
    stop("Need to rescale time to [0,1].")
  }

  control <- do.call(rkhscov.alg.control, control)

  # centering
  if (!centered){
    fit <- gss::ssanova(x~time, alpha=1)
    x <- x - gss:::fitted.ssanova(fit)
  }
  pout <- prep2(time, x, subject, Mmethod=control$Mmethod, tol=control$preptol, nsam=control$nsam)
  r <- ncol(pout$M)
  if (is.null(weight)){
    weight <- rep(1, r)
  } else {
    if (length(weight)<r) weight <- c(rep(1e6*weight[1],r-length(weight)), weight)
    else stop("length of weight does not match!")
    if (any(sort(weight, decreasing=T)!=weight))
      warning("weight should be decreasing to match the order of the eigenvalues!")
  }
  # weight may not work with accelerated gradient descent

  # run accelarated proximal gradient method
  if (is.null(B0v)) B0v <- svec(diag(rep(1e-8,r)))
  res <- rkhscov_pg_cpp(RR=pout$R, RQv=pout$Qv, c=pout$c, lam=lam, gam=gam, B0v=B0v,
                        weight=weight, L=control$L, eta=control$eta,
                        alpha=control$alpha, maxit=control$maxit,
                        traceit=as.logical(control$traceit), tol=control$tol,
                        max_rank=control$max.rank, pos=as.logical(pos),
                        variant=control$variant, cond=control$cond)
  return(list(e=res$e, obj=res$obj, dBs=res$dBs, conv=res$conv, r=r,
              k=length(res$e$values), Ms=pout$Ms, Minv=MASS::ginv(pout$M), t=pout$tts,
              call=match.call(), pos=pos, weight=res$weight, lam=lam, gam=gam))
}

###################
#### k-fold CV ####
###################

gen.groups <- function(n, nfold){
  leave.out <- trunc(n/nfold)
  o <- sample(1:n)
  groups <- vector("list", nfold)
  for (j in (1:(nfold-1))){
    jj <- (1+(j-1)*leave.out)
    groups[[j]]<-(o[jj:(jj+leave.out-1)]) 
  }
  groups[[nfold]] <- o[(1+(nfold-1)*leave.out):n]
  return(groups=groups)
}

# TO-DO: checking on the control parameters
rkhscov.cv.control <- function(lams=NULL, gams=1, lamu=NULL, lam2.rtol=1e-1,
                               lam.min.ratio=1e-8, nlam=20, ngam=20, fine=TRUE, nfine=20){
  return(list(lams=lams, gams=gams, lamu=lamu, lam2.rtol=lam2.rtol,
              lam.min.ratio=lam.min.ratio, nlam=nlam, ngam=ngam, fine=fine,
              nfine=nfine))
}

rkhscov.cv <- function(time, x, subject, nfold=5, weight=NULL, centered=FALSE,
                       ncpu=nfold, pos=TRUE, control.alg=list(),
                       control.cv=list()){
  if ((min(time) < 0) || (max(time) > 1)){
    stop("Need to rescale time to [0,1].")
  }

  control.alg <- do.call(rkhscov.alg.control, control.alg)
  control.cv <- do.call(rkhscov.cv.control, control.cv)

  # centering
  if (!centered){
    fit <- gss::ssanova(x~time, alpha=1)
    x <- x - gss:::fitted.ssanova(fit)
  }
  pout0 <- prep2(time, x, subject, Mmethod=control.alg$Mmethod, tol=control.alg$preptol, nsam=control.alg$nsam)
  r <- ncol(pout0$M)
  n <- length(pout0$Ms)

  # initialize weight
  if (is.null(weight)){
    weight <- rep(1, r)
  } else {
    if (length(weight)<r) weight <- c(rep(1e6*weight[1],r-length(weight)), weight)
    else stop("length of weight does not match!")
    if (any(sort(weight, decreasing=T)!=weight))
      warning("weight should be decreasing to match the order of the eigenvalues!")
  }
  # weight may not work with accelerated gradient descent

  # preparation for CV
  groups <- gen.groups(n, nfold)

  # construction of lam sequence
  if (is.null(control.cv$gams)){
    gams <- seq(1e-3, 1, len=control.cv$ngam)
  } else {
    gams <- control.cv$gams
  }

  if (is.null(control.cv$lams)){
    if (is.null(control.cv$lamu)){
      if (any(gams<1e-20)){
        if (length(gams)==1)
          lamu <- start_lambda(pout0$R, pout0$Qv, pout0$c, weight, 0, pos,
                               control.cv$lam2.rtol) # when we supply gams as c(0)
        else
          stop("gams includes 0, but not of length 1")
      } else {
        lamu <- start_lambda(pout0$R, pout0$Qv, pout0$c, weight, 1, pos,
                             control.cv$lam2.rtol) # When gam>0, it gives the same value
      }
    } else {
        lamu <- control.cv$lamu
    }
    laml <- lamu * control.cv$lam.min.ratio
    lams <- exp(seq(log(laml), log(lamu), len=control.cv$nlam)) # on log scale
  } else {
    lams <- control.cv$lams
  }


  if (ncpu>1) {
    cl <- parallel::makeCluster(min(ncpu,parallel::detectCores()))
    doParallel::registerDoParallel(cl)
  }

  B0v <- svec(diag(rep(1e-8, r)))
  cvs.grid <- matrix(nr=length(lams), nc=length(gams))
  for (k in (1:length(gams))){
    gam <- gams[k]

    cvjs <- foreach::foreach(j=1:nfold, .verbose=T, .packages=c("rkhscovfun")) %dopar% {
      rkhscov_pg_cvj_cpp(pout0$Xs, pout0$Ms,  setdiff(1:n, groups[[j]]),
                         groups[[j]], lams, gam, B0v, weight, L=control.alg$L,
                         eta=control.alg$eta, alpha=control.alg$alpha, maxit=control.alg$maxit,
                         traceit=as.logical(FALSE), tol=control.alg$tol,
                         max_rank=control.alg$max.rank, pos=as.logical(pos),
                         variant=control.alg$variant, cond=control.alg$cond)
    }

    #cvs.se <- apply(sapply(cvjs, function(x){x[,1]}), 1, function(x){sd(x)})/sqrt(nfold)
    cvs <- apply(matrix(sapply(cvjs, function(x){x[,2]}),nc=nfold), 1, sum)/n
    cvs.grid[,k] <- cvs
  }
  cv.best.ind <- arrayInd(which.min(cvs.grid), dim(cvs.grid))

  # fine tuning on lambda
  if ((control.cv$fine && (cv.best.ind[1]>1)) && (cv.best.ind[1]<length(lams))){
    B0v <- svec(diag(rep(1e-8, r)))
    lams2 <- seq(lams[cv.best.ind[1]-1], lams[cv.best.ind[1]+1], len=control.cv$nfine)
    gam <- gams[cv.best.ind[2]]
    cvjs2 <- foreach::foreach(j=1:nfold, .verbose=T, .packages=c("rkhscovfun")) %dopar% {
      rkhscov_pg_cvj_cpp(pout0$Xs, pout0$Ms, setdiff(1:n, groups[[j]]), groups[[j]],
                         lams2, gam, B0v, weight, L=control.alg$L, eta=control.alg$eta,
                         alpha=control.alg$alpha, maxit=control.alg$maxit,
                         traceit=as.logical(FALSE), tol=control.alg$tol,
                         max_rank=control.alg$max.rank, pos=as.logical(pos),
                         variant=control.alg$variant, cond=control.alg$cond)
    }
    # combine
    lams <- c(lams, lams2); oo <- order(lams); lams <- lams[oo]
    cvs2 <- apply(matrix(sapply(cvjs2, function(x){x[,2]}),nc=nfold), 1, sum)/n
    cvs <- c(cvs.grid[,cv.best.ind[2]], cvs2); cvs <- cvs[oo]

    cv.best.ind[1] <- which.min(cvs)
  }

  if (ncpu>1) parallel::stopCluster(cl)

  B0v <- svec(diag(rep(1e-8, r)))
  res <- rkhscov_pg_cpp(RR=pout0$R, RQv=pout0$Qv, c=pout0$c, lam=lams[cv.best.ind[1]],
                        gam=gams[cv.best.ind[2]], B0v=B0v, weight=weight,
                        L=control.alg$L, eta=control.alg$eta, alpha=control.alg$alpha,
                        maxit=control.alg$maxit*2, traceit=as.logical(control.alg$traceit),
                        tol=control.alg$tol, max_rank=control.alg$max.rank,
                        pos=as.logical(pos), variant=control.alg$variant,
                        cond=control.alg$cond)
  sobj <- list(e=res$e, obj=res$obj, dBs=res$dBs, conv=res$conv,
               r=r, k=length(res$e$values), Ms=pout0$Ms, Minv=MASS::ginv(pout0$M), t=pout0$tts,
               weight=res$weight, lam=lams[cv.best.ind[1]], gam=gams[cv.best.ind[2]])

  # flag information only for lambda
  if (cv.best.ind[1]==1) flag <- 1
  else if (cv.best.ind[1]==length(lams)) flag <- 2
  else flag <- 0

  return(list(sobj=sobj, cvs.grid=cvs.grid, cvs=cvs, cv.best.ind=cv.best.ind,
              lams=lams, gams=gams,
              call=match.call(), pos=pos, weight=weight, flag=flag))
}

#####################
#### predictions ####
#####################

fitted.rkhscov <- function(sobj){
  out <- list()
  p <- ncol(sobj$e$vectors)
  B <- sobj$e$vectors %*% diag(sobj$e$values, nr=p, nc=p) %*% t(sobj$e$vectors)
  for (i in (1:length(sobj$Ms)))
    out[[i]] <- sobj$Ms[[i]] %*% B %*% t(sobj$Ms[[i]])
  return(out)
}


compute.A <- function(sobj){
  p <- ncol(sobj$e$vectors)
  B <- sobj$e$vectors %*% diag(sobj$e$values, nr=p, nc=p) %*% t(sobj$e$vectors)
  A <- t(sobj$Minv) %*% B %*% sobj$Minv
  return(A)
}


predict.rkhscov <- function(newtt, sobj){
  A <- compute.A(sobj)
  K <- outer(newtt, unlist(sobj$t), K.sob)
  return(K %*% A %*% t(K))
}

###################################
#### functions related to FPCA ####
###################################

# using average
# more stable
Qmat3 <- function(ss, o=10000, fix.pos=T){
  N <- length(ss)
  tt <- seq(0, 1, len=o)
  Z <- outer(ss, tt, K.sob)
  out <- Z %*% t(Z) / o

  if (fix.pos){
    ee <- eigen(out, symmetric=T)
    ee$values[ee$values<0] <- 0
    out <- ee$vectors %*% diag(ee$values) %*% t(ee$vectors)
  }
  return(out)
}


fpca.rkhscov <- function(sobj, Q=NULL){
  if (is.null(Q)){
    Q <- Qmat3(unlist(sobj$t), fix.pos=T)
  }
  R <- sobj$Minv %*% Q %*% t(sobj$Minv)
  R <- (R + t(R)) / 2
  #Rroot <- matsq_cpp(R, tol=1e-20) # matsq_cpp produce R =Rroot %*% t(Rroot), not Rroot %*% Rroot
  ee <- eigen(R); ee$values[ee$values < 0] <- 0
  Rroot <- ee$vectors %*% diag(sqrt(ee$values)) %*% t(ee$vectors)
  p <- ncol(sobj$e$vectors)
  B <- sobj$e$vectors %*% diag(sobj$e$values, nr=p, nc=p) %*% t(sobj$e$vectors)
  e <- eigen(Rroot %*% B %*% Rroot, symmetric=T)
  U <- t(sobj$Minv) %*% MASS::ginv(Rroot) %*% e$vectors

  rank <- length(sobj$e$values)
  return(list(t=sobj$t, U=U[,1:rank], values=e$values[1:rank]))
}

compute.fpca <- function(tt, fpca.obj){
  Z <- outer(unlist(fpca.obj$t), tt, K.sob)
  return(t(fpca.obj$U) %*% Z)
}



