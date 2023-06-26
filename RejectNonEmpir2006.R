library(spatstat)
library(ggplot2)
library(spatstat)
library(calculus)
library(cubature)
library(stats)
#Fonctions utiles
# Saliency map
x.max <- 10; x.min <- 0; y.max <- 10; y.min <- 0

calc.dist <- function(x,y,X,Y){
  return(sqrt((x-X)^2+(y-Y)^2))
}

alpha1 <- function(x,y){
  res <- ((x-5)^2+(y-5)^2)
  return(res/(50))
}

alpha <- function(x,y){
  res <- (x^2+y^2)
  return(res/200)
}

#Kernel BetaPrime
BP <- function(x,y,X,Y,c,d){
  return(dbetapr(x=calc.dist(x,y,X,Y),shape1 = c,shape2 = d)*1/(2*pi))
}

# Couvrement par des boules
ballcover <- function(x,y,mat,r,l1,beta){
  n <- length(mat)/3 #nombre de points
  l <- min(l1,n)
  PP <- as.ppp(mat[,c(1,2)],W=c(0,10,0,10))
  W <- discs(centres = PP,radii = r)
  #res <- beta^(n-l+1-(1:(max((n-l),1))))*inside.owin(x=mat[1:(max((n-l),1)),1],y=mat[1:(max((n-l),1)),2],disc(radius = r,centre = c(x,y)))*(mat[1:(max((n-l),1)),3]==1)*(l1==l)
  res <- beta^(l-(1:l))*inside.owin(x=mat[(n-l+1):n,1],y=mat[(n-l+1):n,2],disc(radius = r,centre = c(x,y)))*(mat[(n-l+1):n,3]==1)
  #return(sum(res)/sum(beta^((max((n-l),1))+1-(1:(max((n-l),1))))))
  return(sum(res)/sum(beta^(l-(1:l))*inside.owin(x=mat[(n-l+1):n,1],y=mat[(n-l+1):n,2],disc(radius = r,centre = c(x,y)))^(inside.owin(x=mat[(n-l+1):n,1],y=mat[(n-l+1):n,2],disc(radius = r,centre = c(x,y)))!=0)))
}

#Max fonction rééquilibrage
prob.acc <- function(a,b){
  max <- optimize(f=function(x) x^(a)*(1-x)^(b),lower = 0.01,upper = 0.99,maximum = TRUE)$objective
  return(max)
}

# Fonction à échantillonner
target_density <- function(x,y,X,Y,c,d){
  alpha(x,y)*BP(x,y,X,Y,c,d)
}

# Fonction de proposition
proposal_density <- function(x, y,X,Y) {
  dnorm(x, mean = X, sd = 4) * dnorm(y, mean = Y, sd = 4)
}

a.r <- function(X,Y,r,mat,l,beta,a,b,c,d){
  # Borne sup pour target density
  upper_bound <- optim(par = c(X*0.9,Y*0.9),fn=function(M) target_density(M[1],M[2],X,Y,c,d)/proposal_density(M[1],M[2],X,Y),method="L-BFGS-B",lower = c(x.min,y.min),upper = c(x.max,y.max),control=list(fnscale = -1))$value
  found <- FALSE
  #print(upper_bound)
  max.pi <- (a/(a+b))^a*(b/(a+b))^b
  # Algo de rejet
  while(found == FALSE) {
    
    # Générer des candidats
    x_in_window <- FALSE
    y_in_window <- FALSE
    while (x_in_window == FALSE) {
      candidate_x <- rnorm(1, mean = X, sd = 4)
      if(candidate_x<10 & candidate_x>0){
        x_in_window <- TRUE
      }
    }
    while (y_in_window == FALSE) {
      candidate_y <- rnorm(1, mean = Y, sd = 4)
      if(candidate_y<10 & candidate_y>0){
        y_in_window <- TRUE
      }
    }
    
    # Proba d'acceptation
    #acceptance_prob <- target_density(candidate_x, candidate_y,X,Y)/upper_bound/proposal_density(candidate_x, candidate_y,X,Y)
    acceptance_prob <- target_density(candidate_x, candidate_y,X,Y,c,d)/upper_bound/proposal_density(candidate_x, candidate_y,X,Y)
    #max.pi <- prob.acc(a,b)#a,ballcover(candidate_x,candidate_y,mat,r,l,beta))
    prop <- alpha(candidate_x,candidate_y)^a*(1-alpha(candidate_x,candidate_y))^b/max.pi
    #prop <- dbeta(ballcover(candidate_x,candidate_y,mat,r,l,beta),shape1 = a,shape2 = b)
    if(runif(1) < acceptance_prob & runif(1)< prop){
      sample <- c(candidate_x, candidate_y)
      #print(prop)
      found <- TRUE
    }
  }
  return(sample)
}

#ESTIMATION
S1 <- function(mat,c,d){
  n <- nrow(mat)
  X1 <- mat[1:(n-1),1:2]
  X2 <- mat[2:n,1:2]
  Xdist <- calc.dist(X1[1],X1[2],X2[1],X2[2])
  res <- log(1/(2*pi)*(Xdist^(c-2)*(1+Xdist^(-c-d)))/beta(c,d))
  return(sum(res))
}

S2 <- function(mat,a,b){
  n <- nrow(mat)
  X <- mat[2:n,1:2]
  Xalpha <- alpha(X[,1],X[,2])
  res <- a*log(Xalpha)+b*log(1-Xalpha)
  return(sum(res)-(n-1)*a*log(a/(a+b))+b*log(b/(a+b)))
}

integrand <- function(mat,x,y,a,b,c,d){
  Xi <- tail(mat,n=1)[1]; Yi <- tail(mat,n=1)[2]
  dis <- calc.dist(x,y,Xi,Yi)
  max.pi <- (a/(a+b))^a*(b/(a+b))^b
  res <- alpha(x,y)*1/(2*pi)*dis^(c-2)*(1+dis)^(-c-d)/beta(c,d)*alpha(x,y)^a*(1-alpha(x,y))^b/max.pi
  return(res)
}

S3 <- function(mat,a,b,c,d){
  s <- 0
  n <- nrow(mat)-1
  for (i in 2:(n-1)) {
    s <- s+ log(calculus::integral(f=function(u,v) integrand(mat[1:i,],u,v,a,b,c,d),bounds=list(u=c(0,10),v=c(0,10)),relTol = 0.1,absTol = 1e-2)$value)
  }
  return(s)
}

loglik <- function(mat,a,b,c,d){
  return(S1(mat,c,d)+S2(mat,a,b)-S3(mat,a,b,c,d))
}







