#' Resampling schemes
#' 
#' Bootstrap a vector according to multiple resampling schemes: resampling, block resampling, Wild bootstrap. 
#' @param x A numeric vector
#' @param boot.scheme The type of resampling scheme used, see Details
#' @param seed the seed used, see \code{\link{set.seed}}
#' @param block.size for the \code{resample_block} scheme, the size of the blocks.
#'  
#' @details This function offers various bootstrap/resampling schemes:
#' \describe{
#'   \item{resample}{Resampling with replacement}
#'   \item{resample_block}{Resampling contiguous observations (blocks) with replacement. Use argument \code{block.size}}
#'   \item{wild1}{Wild bootstrap: do not resample, but add a N(0,1) distribution to each value}
#'   \item{wild12}{Wild bootstrap: same, but add instead -1 or 1. }
#' }
#' @export

resample_vec <- function(x, boot.scheme=c("resample","resample_block", "wild1", "wild2", "check"), seed = NULL, 
                   block.size=2){
  
  boot.scheme <- match.arg(boot.scheme)
  
  ## everything is matrix here
  X <-  as.matrix(x)
  t <-  nrow(X)

  if(!is.null(seed)) set.seed(seed) 
  X_boot <- switch(boot.scheme, 
                  "resample" =  X[sample(seq_len(t), replace=TRUE),,drop= FALSE], 
                  "resample_block" = X[sample.block(seq_len(t), block.size = block.size), , drop= FALSE], 
                  "wild1" = X+rnorm(t), 
                  "wild2" = X+sample(c(-1,1), size = t, replace=TRUE),
                  "check" = X)
  if(!is.matrix(x)) X_boot <- X_boot[,1]
  X_boot
}

if(FALSE) {
  head(resample_vec(x=lynx, seed = 123))
  head(resample_vec(cbind(lynx, lynx), seed = 123))
  
  head(resample_vec(cbind(lynx, lynx), seed = 123, boot.scheme = "resample_block"))
}

### for backwards compatibility
resample_vec

#### Small function to sample in block: sample.block 
sample.block.alter <- function(x, size=length(x), block.size=2){
  n <- length(x)
  n.blocks <- ceiling((n/block.size))
  
  block <- rep(1:n.blocks, each=block.size)[1:n]
  
  # sample the blocks:
  sam <- sample(1:n.blocks, size=n.blocks, replace=TRUE, prob = table(block))
  
  get_pos <- function(x) which(block %in% x)
  get_pos(sam[4])
  pos <- unlist(lapply(sam, get_pos))
  length(pos)
  
  ## match first position in each block
  m <- match(sam, block)
  m2 <- rep(m, each=block.size)[1:n]   ## repeat "first
  m3 <- m2 +rep(0:(block.size-1), n)[1:n] # match now next
  
  if(any(m3>n)) m3[m3>n] <- m3[m3>n]-(max(m3)-n)
  x[m3]
}

sample.block <- function(x, size=length(x), block.size=2){
  n <- length(x)
  n.blocks <- ceiling((n/block.size))
  
  block <- rep(1:n.blocks, each=block.size)[1:n]
  
  # sample the blocks:
  sam <- sample(1:n.blocks, size=n.blocks, replace=TRUE, prob = table(block))
  
  ## match first position in each block
  m <- match(sam, block)
  m2 <- rep(m, each=block.size)[1:n]   ## repeat "first
  m3 <- m2 +rep(0:(block.size-1), n)[1:n] # match now next
  
  if(any(m3>n)) m3[m3>n] <- m3[m3>n]-(max(m3)-n)
  x[m3]
}


if(FALSE) {
  M <-  matrix(1:20, ncol = 2)
  sample.block(x=M)
  sample.block(x=1:13, block.size = 3)
}
