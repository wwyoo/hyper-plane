#Author: William Weimin Yoo

cross <- function(x, y)
#brief: computes cross product of two vectors,
#       this is not the crossproduct used in R
#       but the vector product a x b 

#input: two numeric vectors, each of length three

#output: a numeric vector that is orthogonal to the
#        plane spanned by the input vectors

{
  z1 <- x[2L] * y[3L] - x[3L] * y[2L]
  z2 <- x[3L] * y[1L] - x[1L] * y[3L]
  z3 <- x[1L] * y[2L] - x[2L] * y[1L]
  return(c(z1, z2, z3))
} 

collision <- function(a, b)
#brief: check whether two intervals intersect, with 
#       the intervals formed by projecting convex hull onto axis

#input: two numeric vectors, where each vector contains axis coordinates
#        of the points spanning the convex hull

#output: returns TRUE if intervals intersect and FALSE otherwise

{
  interb <- c(min(a), max(a))  #base interval
  interc <- c(min(b), max(b))  #comparision interval
  #do they overlap?
  if((interb[1L] <= interc[2L]) & (interb[2L] >= interc[1L]))  
  {
    collision <- TRUE  #if yes, then collision occurs
  }
  else{
    collision <- FALSE
  }
  return(collision)
}

nlinclass <- function(X)
#brief: computes total number of possible binary classifications
#       given a set of points in R^3

#input: an n x 3 numric matrix where each row is a point in R^3

#output: an integer for the total number of possible 
#         binary classifications

#note: output is 2^n - number of non linearly separable groups.
#      Two groups of points are linearly separable iff the 
#      convex hulls spanned by the respective groups do not collide.
#      Collision is detected by projecting the 3D convex hulls into
#      1D intervals, and collision occurs if the intervals intersect 
#      in all dimensions

{
  n <- dim(X)[1L]  #total number of data points

  #2 points and less always be shattered
  if(n <= 2L)
  { 
    return(2 ^ n)
  }

  # 3 points will be shattered unless there are collinear
  if(n == 3L)
  {  
    # 3 points are collinear if the convex hull/triangular
    #spanned by those points have area zero
    d <- cross(X[2L, ] - X[1L, ], X[1L, ] - X[3L,])
    d <- sqrt(crossprod(d))
      if(d == 0)
      {
        return(2 ^ n - 2)
      }
      else{
        return(2 ^ n)
      }
  }

  if(n >= 4L)
  { 
    #project convex hull onto the x-y, x-z, and y-z planes
    xyhull <- chull(X[, -3L])
    xzhull <- chull(X[, -2L])
    yzhull <- chull(X[, -1L])

    #single point case
    #if you are in all the convex hull in all projections
    #then you cannot be separated by planes
    cm1 <- setdiff(1:n, xyhull)
    cm2 <- setdiff(1:n, xzhull)
    cm3 <- setdiff(1:n, yzhull)
    onepoint <- intersect(intersect(cm1, cm2), cm3)

    #need to evaluate up to |_n/2_| only because of symmetry
    upper <- floor(n / 2L)  
    count <- rep(0,  upper - 1L)  #count collision events
    for(i in 2L:upper)
    {
      all.combin <- combn(n, i)  #generate all n choose i combinations
      tol.combin <- dim(all.combin)[2L]  #total of possible combinations
      collision.tol <- rep(0, tol.combin)  #store collision events

      for(j in 1L:tol.combin)
      { 
        #x-y plane
        #project convex hull of the jth point combination onto x-y plane
        xyhull <- chull(X[all.combin[, j], -3L]) 

        #project convex hull of the rest onto x-y plane
        xyhullrest <- chull(X[-all.combin[, j], -3L]) 

          #project jth point hull onto x axis
          Pxyhull.x <- X[all.combin[, j], -3L][xyhull, -2L]  

          #project the complement hull onto x axis
          Pxyhullrest.x <- X[-all.combin[, j], -3L][xyhullrest, -2L]

          #did they collide in 1D x axis?
          collisionxy.x <- collision(a = Pxyhull.x, b = Pxyhullrest.x)

          #project jth point hull onto y axis
          Pxyhull.y <- X[all.combin[, j], -3L][xyhull, -1L]

          #project the complement hull onto y axis
          Pxyhullrest.y <- X[-all.combin[, j], -3L][xyhullrest, -1L]

          #did they collide in 1D y axis?
          collisionxy.y <- collision(a = Pxyhull.y, b = Pxyhullrest.y)

        #did they collide in the x-y plane?
        collision.xy <- collisionxy.x & collisionxy.y
        
        #the comments below will be similar to the one above 
        #suppress for brevity
        #x-z plane
        xzhull <- chull(X[all.combin[, j], -2L])
        xzhullrest <- chull(X[-all.combin[, j], -2L])

          Pxzhull.x <- X[all.combin[, j], -2L][xzhull, -3L]
          Pxzhullrest.x <- X[-all.combin[, j], -2L][xyhullrest, -3L]
          collisionxz.x <- collision(a = Pxzhull.x, b = Pxzhullrest.x)

          Pxzhull.z <- X[all.combin[, j], -2L][xzhull, -1L]
          Pxzhullrest.z <- X[-all.combin[, j], -2L][xzhullrest, -1L]
          collisionxz.z <- collision(a = Pxzhull.z, b = Pxzhullrest.z)

        collision.xz <- collisionxz.x & collisionxz.z

        #y-z plane
        yzhull <- chull(X[all.combin[, j], -1L])
        yzhullrest <- chull(X[-all.combin[, j], -1L])

          Pyzhull.y <- X[all.combin[, j], -1L][yzhull, -3L]
          Pyzhullrest.y <- X[-all.combin[, j], -1L][yzhullrest, -3L]
          collisionyz.y <- collision(a = Pyzhull.y, b = Pyzhullrest.y)

          Pyzhull.z <- X[all.combin[, j], -1L][yzhull, -2L]
          Pyzhullrest.z <- X[-all.combin[, j], -1L][yzhullrest, -2L]
          collisionyz.z <- collision(a = Pyzhull.z, b = Pyzhullrest.z)

        collision.yz <- collisionyz.y & collisionyz.z

      #collision in all three planes for the jth point combination?
      collision.tol[j] <- as.numeric(collision.xy & collision.xz & collision.yz)
      }
      count[i - 1L] <- sum(collision.tol)  #total counts of non linearly separable groups
    }
    #factors of two due to duality
    #for even number of points, there is only one middle value
    #for odd number of points, there are two replicated middle values
    if(n %% 2 == 0){
      return(2 ^ n - 2 * sum(count[-upper + 1L]) - 2 * length(onepoint) - count[upper - 1L])
    }
    else{
      return(2 ^ n - 2 * sum(count) - 2 * length(onepoint)) 
    }
  }
}

set.seed(1000)
X <- matrix(rnorm(30, 0,1), 10, 3)

#The following code is used to plot the collision graphs as shown in Figure 1
#Figure 1(a)
par(mar = c(5, 4, 0, 0) + 0.1)
plot(X[, -3], xlab = "x", ylab = "y")  #project points onto x-y plane
identify(X[,-3])  #label the points
combin <- c(1, 2)  #points 1-2
xyhull <- chull(X[combin, -3L]) #convex hull of the points 1-2
xyhullrest <- chull(X[-combin, -3L]) #convex hull of the rest
Pxyhull <- X[combin, -3L][xyhull, ]  #project 1-2 hull onto x-y plane
Pxyhullrest <- X[-combin, -3L][xyhullrest, ]  #project the hull of the rest onto x-y plane
lines(rbind(Pxyhull, Pxyhull[1, ]))  #plot 1-2 convex hull
lines(rbind(Pxyhullrest, Pxyhullrest[1, ]))  #plot the convex hull of the rest

#Figure 1(b)
x11()
par(mar = c(5, 4, 0, 0) + 0.1)
plot(X[, -2], xlab = "x", ylab = "z")
identify(X[,-2])
xzhull <- chull(X[combin, -2L])
xzhullrest <- chull(X[-combin, -2L])
Pxzhull <- X[combin, -2L][xzhull, ]
Pxzhullrest <- X[-combin, -2L][xzhullrest, ]
lines(rbind(Pxzhull, Pxzhull[1, ]))
lines(rbind(Pxzhullrest, Pxzhullrest[1, ]))

#Figure 1(c)
x11()
par(mar = c(5, 4, 0, 0) + 0.1)
plot(X[, -1], xlab = "y", ylab = "z")
identify(X[,-1])
yzhull <- chull(X[combin, -1L])
yzhullrest <- chull(X[-combin, -1L])
Pyzhull <- X[combin, -1L][yzhull, ]
Pyzhullrest <- X[-combin, -1L][yzhullrest, ]
lines(rbind(Pyzhull, Pyzhull[1, ]))
lines(rbind(Pyzhullrest, Pyzhullrest[1, ]))
