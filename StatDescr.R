prop.infec <- function(mat,R){
  n <- length(mat)/3#*is.marked+length(mat)/2*(1-is.marked)
  res <- 0
  time <- inside.owin(x=mat[1:(n-1),1],y=mat[1:(n-1),2],disc(radius = R,centre = c(mat[n,1],mat[n,2])))
  if(length(time) >0){
    res <- sum(mat[time,3])/length(time)
  }
  return(res)
}

first.contact <- function(mat){
  return(tail(nndist(mat[,1],mat[,2]),n=1))
}

proper.zone <- function(r,mat,W){
  size <- length(mat)/3
  # B A RENDRE EMPTY
  B <- disc(radius=r,centre = c(mat[1,1],mat[1,2]))
  for (i in 1:(size-1)) {
    if(mat[i,3]==1){
      B <- union.owin(B,disc(radius=r,centre = c(mat[i,1],mat[i,2])))
    }
  }
  return( area.owin(intersect.owin(setminus.owin(disc(radius=r,centre = c(mat[size,1],mat[size,2])),B),W))/area.owin(intersect.owin(disc(radius=r,centre = c(mat[size,1],mat[size,2])),W)))
}

ball.union.cover <- function(r,mat,W){
  size <- length(mat)/3
  # B A RENDRE EMPTY
  B <- disc(radius=r,centre = c(mat[1,1],mat[1,2]))
  for (i in 1:size) {
    if(mat[i,3]==1){
      B <- union.owin(B,disc(radius=r,centre = c(mat[i,1],mat[i,2])))
    }
  }
  return(area.owin(intersect.owin(B,W))/area.owin(W))
}

scanpath.length <- function(mat){
  ispos <- FALSE
  dist <- 0
  k <- 0
  size <- length(mat)/3
  while ((ispos == FALSE) & (size-k-1>0)) {
    dist <- dist + sqrt((mat[size-k,1]-mat[size-k-1,1])^2+(mat[size-k,2]-mat[size-k-1,2])^2)
    ispos <- mat[size-k-1,3]
    k <- k+1
  }
  return(dist)
}

cumul.recurr <- function(r,mat){
  # INSIDE OWIN
  n <- length(mat)/3
  N <- 0
  for (i in 2:n) {
    N <- N+sum(inside.owin(mat[1:(i-1),1],mat[1:(i-1),2],disc(radius = r,centre = c(mat[i,1],mat[i,2]))))
  }
  return(N)
}