
###--------------------------------------------------------------###
###        Create an artificial data from the TDBT model         ###
###--------------------------------------------------------------###

## INPUT:
# num.freq: The total number of simulated matches between entities.
# K:        The dimensionality of the "worth" vector f.
# w:        A K×1 vector representing the relative weight of each dimension in f.
# F:        A K×N matrix where each column f_i represents the worth vector of player i.

## OUTPUT:
# An N×N matrix where the (i, j) entry indicates that player i defeats player j.

Create.Artificial.Data <- function(num.freq = 1000, w0, F0) {
  N <- ncol(F0) # number of entities
  K <- nrow(F0) # dimensionality of each entity
  scores <- as.numeric(exp( w.true %*% F.true ))

  # 結果を格納する N×N 行列を初期化（対角は NA）
  result <- matrix(NA_integer_, nrow = N, ncol = N)
  rownames(result) <- colnames(result) <- 
    if (!is.null(colnames(F0))) colnames(F0) else paste0("Entity", 1:N)
  
  # Simulate for each pair (i, j) 
  set.seed(73)
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      p_ij <- scores[i] / (scores[i] + scores[j])
      win.freq <- rbinom(n = 1, size = num.freq, prob = p_ij)
      result[i, j] <- win.freq
      result[j, i] <- num.freq - win.freq
    }
  }
  return(result)
}
