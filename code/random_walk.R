random_walk_ranking <- function(AdjMatrix, seeds, r = 0.3, tol = 1e-6, max_iter = 1000) {
  n <- nrow(AdjMatrix)
  P <- AdjMatrix / rowSums(AdjMatrix)
  p0 <- rep(0, n)
  seed_index <- which(rownames(AdjMatrix) %in% seeds)
  p0[seed_index] <- 1 / length(seed_index)
  p <- p0
  for (i in 1:max_iter) {
    p_new <- (1 - r) * (t(P) %*% p) + r * p0
    if (sum(abs(p_new - p)) < tol) {
      break
    }
    p <- p_new
  }
  result <- data.frame(Gene = rownames(AdjMatrix),
                       Score = as.numeric(p))
  result <- result[order(-result$Score), ]
  return(result)
}


