# generic to get co-occurrence matrix from objects of various classes
#' Calculate Co-Occurrence Matrix
#'
#' Calculate generalized co-occurrence matrix from a variety of objects,
#' currently including fitness landscapes stored as a `FitLandDF` instance from
#' the `fitscape` package.
#'
#' @name comat
#' @param x object upon which co-occurrence matrix will be calculated
#' @param discrete function that discretizes object
#' @param neighbor function that returns TRUE if two numeric vectors are within
#'   acceptable distance of one another or a single-element `character` vector
#'   that describes how to identify acceptable neighbors/offsets
#' @param normalize function that normalizes the co-occurrence matrix
#' @param ... additional arguments
#' @return co-occurrence matrix
#' @export
#' @examples
#' # create fitness landscape as instance of FitLandDF object
#' a <- round(runif(64))
#' a <- array(a, dim = rep(4, 3))
#' my_landscape <- fitscape::FitLandDF(a)
#'
#' # calculate co-occurrence matrix using:
#' #   Manhattan distance of 1
#' #   discretization into 2 equal-sized buckets
#' #   normalization: multiply all elements so that sum of matrix equals unity
#' comat <- get_comatrix(my_landscape,
#'                       discrete = equal_discrete(2),
#'                       neighbor = manhattan(1))
#'
#' # print co-occurrence matrix
#' print(comat)
#'
get_comatrix <- function(x, ...) {
  UseMethod("get_comatrix")
}

# default not defined
#' @rdname comat
#' @export
get_comatrix.default <- function(x, ...) {
  stop("default behavior for get_comatrix not yet defined")
}

# for igraph objects
# 1. discretize node values to `nlevels` levels
# 2. initialize co-occurrence matrix with all -1 (zero if actually processed)
# 3. count co-occurrences
# 4. normalize and return
#' @rdname comat
#' @export
get_comatrix.igraph <- function(x,
                                values = seq(length(x)), # default nodes
                                discrete = equal_discrete(2),   # currently a function from factory - need to export factory
                                normalize = normalize_glcm,  # currently a function
                                nlevels = length(x), # default all
                                   ...) { #   doesn't do anything
  
  # Generate adjacency list from igraph object
  adj_list = igraph::as_adj_list(x)
  
  # discretize node values (`values`) to `nlevels` levels, equal to integers 1:nlevels
  Value <- discretize(values, nlevels)
  
  nlevels <- Value %>% unique %>% length
  
  # Initialize comat
  comat <- matrix(-1, nrow = nlevels, ncol = nlevels)
  
  # count co-occurrences
  for (i in 1:length(adj_list)) {
    for (j in 1:length(adj_list[[i]])) {
      ii = Value[i]
      jj = Value[adj_list[[i]][j]]
      comat[ii, jj] <- comat[ii,jj] + 1
    }
  }
  comat = comat+1
  
  # normalize
  comat <- normalize(comat)
  
  # return co-occurrence matrix
  return(comat)
}

# for fitness landscapes
# 1. discretize FL to `nlevels` levels
# 2. initialize co-occurrence matrix with all -1 (zero if actually processed)
# 3. count co-occurrences
# 4. normalize and return
#' @rdname comat
#' @export
get_comatrix.FitLandDF <- function(x,
                                   discrete = equal_discrete(2),             # currently a function from factory - need to export factory
                                   neighbor = manhattan(1),                  # currently a function from factory - need to export factory
                                   normalize = normalize_glcm,               # currently a function
                                   ...) { #   doesn't do anything
  # discretize FL (`x`) to `nlevels` levels, equal to integers 1:nlevels
  x$Value <- discrete(x$Value)

  # initialize co-occurrence matrix with all -1 (zero if actually processed)
  nlevels <- x$Value %>% unique %>% length
  comat <- matrix(-1, nrow = nlevels, ncol = nlevels)

  # count co-occurrences
  for (i in seq_len(nrow(comat))) {
    for (j in seq_len(ncol(comat))) {
      comat[i, j] <- count_element_occur(x, i, j, neighbor = neighbor)
    }
  }

  # make symmetric
  comat <- comat + t(comat)

  # normalize
  comat <- normalize(comat)

  # return co-occurrence matrix
  return(comat)
}
