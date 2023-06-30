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

max <- 1.6
alpha <- function(x, y) {
  return((1.5*exp(-0.1*((x-3)^2+(y-3)^2))+exp(-0.1*((x-9)^2+(y-9)^2)))/max)
}

alpha1 <- function(x,y){
  res <- (x^2+y^2)
  return(res/200)
}

#Kernel BetaPrime
BP <- function(x,y,X,Y,c,d){
  return(dbetapr(x=calc.dist(x,y,X,Y),shape1 = c,shape2 = d)*1/(2*pi))
}

# Couvrement p_2 par des boules
ballcover <- function(x,y,mat,r,l1,beta){
  n <- nrow(mat) #nombre de points
  l <- min(l1,n)
  PP <- as.ppp(mat[,c(1,2)],W=c(0,10,0,10))
  W <- discs(centres = PP,radii = r)
  ind1 <- which((inside.owin(x=mat[(n-l+1):n,1],y=mat[(n-l+1):n,2],disc(radius = r,centre = c(3,3))) == 1) & ((mat[(n-l+1):n,3]==1)))
  #res <- beta^(n-l+1-(1:(max((n-l),1))))*inside.owin(x=mat[1:(max((n-l),1)),1],y=mat[1:(max((n-l),1)),2],disc(radius = r,centre = c(x,y)))*(mat[1:(max((n-l),1)),3]==1)*(l1==l)
  res <- ((inside.owin(x=mat[(n-l+1):n,1],y=mat[(n-l+1):n,2],disc(radius = r,centre = c(x,y))) == 1) & ((mat[(n-l+1):n,3]==1)))
  #return(sum(res)/sum(beta^((max((n-l),1))+1-(1:(max((n-l),1))))))
  return(sum(res)/(sum(res) + sum(beta^(l-(1:l))*((inside.owin(x=mat[(n-l+1):n,1],y=mat[(n-l+1):n,2],disc(radius = r,centre = c(x,y))) == 1) & ((mat[(n-l+1):n,3]==0)))))^(sum(inside.owin(x=mat[(n-l+1):n,1],y=mat[(n-l+1):n,2],disc(radius = r,centre = c(x,y))))!=0))
}

# Couvrement par des boules
ballcover1 <- function(x,y,mat,r,l1,beta){
  n <- length(mat)/3 #nombre de points
  l <- min(l1,n)
  PP <- as.ppp(mat[,c(1,2)],W=c(0,10,0,10))
  W <- discs(centres = PP,radii = r)
  #res <- beta^(n-l+1-(1:(max((n-l),1))))*inside.owin(x=mat[1:(max((n-l),1)),1],y=mat[1:(max((n-l),1)),2],disc(radius = r,centre = c(x,y)))*(mat[1:(max((n-l),1)),3]==1)*(l1==l)
  res <- beta^(l-(1:l))*inside.owin(x=mat[(n-l+1):n,1],y=mat[(n-l+1):n,2],disc(radius = r,centre = c(x,y)))*(mat[(n-l+1):n,3]==1)
  #return(sum(res)/sum(beta^((max((n-l),1))+1-(1:(max((n-l),1))))))
  return(sum(res)/sum(beta^(l-(1:l))*inside.owin(x=mat[(n-l+1):n,1],y=mat[(n-l+1):n,2],disc(radius = r,centre = c(x,y)))^(inside.owin(x=mat[(n-l+1):n,1],y=mat[(n-l+1):n,2],disc(radius = r,centre = c(x,y)))!=0)))
}

lamb <- function(x,y,mat,r,tau){
  N <- sum(inside.owin(x=mat[,1],y=mat[,2],disc(radius = r,centre = c(x,y))))
  res <- (1-exp(-tau*N))/(1+exp(-tau*N))
  return(res)
}

#Max fonction rééquilibrage
prob.acc <- function(a,b){
  max <- optimize(f=function(x) x^(a)*(1-x)^(b),lower = 0,upper = 1,maximum = TRUE)$objective
  return(max)
}

# Fonction à échantillonner
target_density <- function(x,y,X,Y,c,d){
  alpha(x,y)*BP(x,y,X,Y,c,d)
}

# Fonction de proposition
proposal_density <- function(x, y,X,Y) {
  dnorm(x, mean = X, sd = 2.5) * dnorm(y, mean = Y, sd = 2.5)
}

a.r <- function(X,Y,r,mat,l,beta,a,b,c,d,tau){
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
      candidate_x <- rnorm(1, mean = X, sd = 2.5)
      if(candidate_x<10 & candidate_x>0){
        x_in_window <- TRUE
      }
    }
    while (y_in_window == FALSE) {
      candidate_y <- rnorm(1, mean = Y, sd = 2.5)
      if(candidate_y<10 & candidate_y>0){
        y_in_window <- TRUE
      }
    }
    
    # Proba d'acceptation
    #acceptance_prob <- target_density(candidate_x, candidate_y,X,Y)/upper_bound/proposal_density(candidate_x, candidate_y,X,Y)
    acceptance_prob <- target_density(candidate_x, candidate_y,X,Y,c,d)/upper_bound/proposal_density(candidate_x, candidate_y,X,Y)
    #max.pi <- prob.acc(a,b)#a,ballcover(candidate_x,candidate_y,mat,r,l,beta))
    conv <- lamb(candidate_x,candidate_y,mat,r,tau)
    prop <- conv*ballcover(candidate_x,candidate_y,mat,r,l,beta)^a*(1-ballcover(candidate_x,candidate_y,mat,r,l,beta))^b/max.pi+1-conv
    #prop <- dbeta(ballcover(candidate_x,candidate_y,mat,r,l,beta),shape1 = a,shape2 = b)
    if(runif(1) < acceptance_prob & runif(1)< prop){
      sample <- c(candidate_x, candidate_y)
      found <- TRUE
      prob <- NA
      n <- nrow(mat)
      l1 <- min(l,n)
      if (sum(inside.owin(x=mat[(n-l1+1):n,1],y=mat[(n-l1+1):n,2],disc(radius = r,centre = c(candidate_x,candidate_y))))>4){
        prob <- ballcover(candidate_x,candidate_y,mat,r,l,beta)
      }
      #print(prop)
    }
  }
  return(list(point=sample,prob=prob))
}