DE_Score <- function(normal_sample,cancer_sample,value){

  gene_exp1 <- matrix(data=0,nrow=nrow(cancer_sample),ncol=2)
  gene_exp1[,1] <- rownames(cancer_sample)
  for(i in 1:nrow(cancer_sample))
  {
    m1 <- mean(as.numeric(normal_sample[i,]))
    m2 <- mean(as.numeric(cancer_sample[i,]))
    var1 <- sd(as.numeric(normal_sample[i,]))
    var2 <- sd(as.numeric(cancer_sample[i,]))
    var_sum <- var1^2+var2^2
    sum1 <- (m1-m2)^2/(4*var_sum)
    sum2 <-  0.5*log(var_sum/(2*var1*var2))
    gene_exp1[i,2] <- sum1+sum2
  }
  dele_index <- which(gene_exp1[,2]=='Inf')
  gene_exp1[dele_index,2] <- 0
  del <- as.numeric(gene_exp1[,2])
  nlist <- which(del>value)
  gene_exp1 <- gene_exp1[nlist,]
  gene_exp <- gene_exp1[order(gene_exp1[,2],decreasing = T),]
  return(gene_exp)
}
