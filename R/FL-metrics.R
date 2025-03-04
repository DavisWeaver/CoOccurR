#### All GLCM functions are applied to normalised co-occurrence matrices
# These are square matrices of size n_levels x n_levels
# where n_levels is the number of discretized fitness levels
# To generate a random co-occurrence matrix in a system with 3 levels
# x = matrix(sample.int(9, size = 9, replace = TRUE), ncol = 3)
# x = x / sum(x)

#### P x-y ####
# This is the marginal probability of specific differences between x and y
p_xsuby.FitLandDF <- function(x, axis, ...) {
  nlevels=dim(x)[1]
  #p_k=c()
  p_k = rep(0, nlevels)
  for (i in 1:nlevels){
    for (j in 1:nlevels){
      k = abs(i-j)
      p_k[k+1] = p_k[k+1] + x[i,j]
    }
  }
  return(p_k)
}


#### P x+y ####
# This is the marginal probability of specific sums of x and y
p_xplusy.FitLandDF <- function(x, axis, ...) {
    nlevels=dim(x)[1]
    #p_k=c()
    p_k = rep(0, 2*nlevels)
    for (i in 1:nlevels){
      for (j in 1:nlevels){
        k = i+j
        p_k[k] = p_k[k] + x[i,j]
      }
    }
    return(p_k)
  }

##### CORRELATION #####
# COMMENTED OUT BECAUSE UNABLE TO FIND WEIGHTED VAR FUNCTION
# Correlation is the correlation between neighboring fitness values
#correlation.FitLandDF <- function(x, ...) {
#  nlevels=dim(x)[1]
#  mu_x = weighted.mean(1:nlevels, rowSums(x))
#  mu_y = weighted.mean(1:nlevels, colSums(x))
#  sd_x = sqrt(weighted.var(1:nlevels, rowSums(x)))
#  sd_y = sqrt(weighted.var(1:nlevels, colSums(x)))
#  corr =  0
#  for (i in 1:nlevels){
#    for (j in 1:nlevels){
#      corr = corr + ((i-mu_x)*(j-mu_y)) * x[i,j]
#    }
#  }
#  return((corr/(sd_x*sd_y))[1,1])
#}

#### SUM AVERAGE #####
# Sum average is the mean value of the sum in marginal distribution px+y
sum_avg.FitLandDF <- function(x, ...) {
  p_xplusy = p_xplusy.FitLandDF(x)[-1]
  nlevels=dim(x)[1]
  i = 2:(2*nlevels)
  sum_avg = i * p_xplusy
  return(sum(sum_avg))
}

#### SUM ENTROPY #####
# Sum entropy is the entropy of marginal distribution px+y
sum_entropy.FitLandDF <- function(x, ...) {
  xplusy = p_xplusy.FitLandDF(x)[-1] #Starts at index 2
  sum_entropy = - xplusy * log (xplusy)
  return(sum(sum_entropy, na.rm=TRUE))
}

#### DIFFERENCE ENTROPY #####
# Difference entropy is the entropy of marginal distribution px-y
diff_entropy.FitLandDF <- function(x, ...) {
  xsuby = p_xsuby.FitLandDF(x)
  diff_entropy = - xsuby * log (xsuby)
  return(sum(diff_entropy))
}

#### DISSIMILARITY #####
# Dissimilarity is the weighted sum of gray-level differences
diss.FitLandDF <- function(x, ...) {
  diss = 0
  nlevels=dim(x)[1]
  for (i in 1:nlevels){
    for (j in 1:nlevels){
      diss = diss + abs(i-j)*x[i,j]
    }
  }
  return(diss)
}
