# mean function
mu <- function(tt){
  return( 3*sin(3*pi*(tt+0.5)+2*tt^3) )
}

# eigenfunctions
phis <- function(tt, k){
  if (k%%2){
    return(sqrt(2)*cos(2*ceiling(k/2)*pi*tt))
  } else {
    return(sqrt(2)*sin(2*ceiling(k/2)*pi*tt))
  }
}

# square root of eigenvalues
slams <- function(k){
  return(1/(k+1))
}

generate.demo.data <- function(n, m, sigma, nk){
  # simulate observations
  times <- NULL
  y <- NULL
  for (i in (1:n)){
    tt <- runif(m)
    times <- c(times, tt)
    yy <- mu(tt)
    for (k in (1:nk)){
      yy <- yy + rnorm(1) * slams(k) * phis(tt, k) 
    }
    yy <- yy + rnorm(m) * sigma
    y <- c(y, yy)
  }
  ids <- rep(1:n, each=m)
  return(list(times=times, ids=ids, y=y))
}
