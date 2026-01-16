library(readr)
library(igraph)
library(irlba)

####Step1.Processing expression data####
load( ".../data/BRCA/brca_gene_exp.RData")
R_mut <- brca_gene_exp 
mrna_exp <- R_mut 
gene_exp1 <- matrix(data=0,nrow=nrow(mrna_exp),ncol=6)
gene_exp1[,1] <- mrna_exp[,1]
colnames(gene_exp1) <- c("gene","V1","V2","V3","V4","V5")
rownames(mrna_exp) <- mrna_exp[,1]
mrna_exp <- mrna_exp[,-1]

##Expression matrix filtering
means<-apply(mrna_exp,1,mean) 
var<-apply(mrna_exp,1,var) 
sd<-apply(mrna_exp,1,sd) 
FF<-var/(1+var) 
ET<-means+2.5*sd*FF 
gene_exp1[,2] <- means
gene_exp1[,3] <- var
gene_exp1[,4] <- sd
gene_exp1[,5] <- FF
gene_exp1[,6] <- ET
write.table(gene_exp1,'.../code/gene_exp1_brca.txt',quote=FALSE, row.names=FALSE)
write.csv(mrna_exp,'.../code/mrna_exp_brca.csv')
write.table(means,'.../code/means_brca.txt',quote=F,row.names=F,col.names=F)

##Import the results of remrna_exp.m:remrna_exp_brca.txt
#For convenience in running, we have already placed this file in the folder
remrna_exp <- read.table('.../code/remrna_exp_brca.txt',sep=',')
rowna <- rownames(mrna_exp)
colna <- colnames(mrna_exp)
rownames(remrna_exp)<-rowna
colnames(remrna_exp)<-colna
mrna_col_name <- colnames(remrna_exp)
normal_exp <- grep('.11A|11B',mrna_col_name)
nor <- length(normal_exp)
cer <- nor + 1
total <- length(mrna_col_name)
brca_normal_s <- remrna_exp[, 1:nor]
brca_cancer_s <- remrna_exp[, cer:total]

##Calculate differential expression score
source(".../code/DE_Score.R")
load(".../data/BRCA/brca_cancer_si.rda")
load(".../data/BRCA/brca_normal_si.rda")
DE_05 <- DE_Score(brca_normal_s,brca_cancer_s,0.5)
mirna_value <- DE_Score(brca_normal_si,brca_cancer_si,0.5)

####Step2.Multilayer network construction####
load(".../data/BRCA/brca_mut.RData")
source(".../code/construct_layer.R")
load(".../data/network.rda")
load(".../data/protein_info.rda")
construct_layer(brca_mut,network,protein_info,DE_05) #Export 6 CSV files

####Step3.Calculate control ability score and regulatory potential score####
##Calculate control ability score
load(".../data/mir_gene_network.rda")
colnames(mir_gene_network) <- mir_gene_network[1,] 
mir_gene_network <- mir_gene_network[-1,] 
data <- read_csv(".../data/layer_2_brca.csv")
two_gene_list <- data.frame(
  Gene = data[,1],
  Score1 = 0,
  Score2 = 0,
  Score3 = 0,
  Score4 = 0,
  Score5 = 0,
  Score6 = 0,
  Score7 = 0,
  Score8 = 0,
  Score9 = 0,
  Score10 = 0,
  Score11 = 0,
  Score12 = 0,
  Score13 = 0,
  Score14 = 0,
  Score15 = 0,
  Score16 = 0,
  stringsAsFactors = FALSE
)
for (i in 1:nrow(mirna_value)) {
  g_list <- which(mir_gene_network[,2]==mirna_value[i,1])
  if(length(g_list)==0)next
  gn_list <- mir_gene_network[g_list,1]
  for (j in 1:length(gn_list)) {
    r_index <- which(two_gene_list[,1]==gn_list[j])
    if(length(r_index)==0)next
    value1 <- as.numeric(two_gene_list[r_index,4])
    two_gene_list[r_index,4] <- value1 + as.numeric(mirna_value[i,2])
  }
}
colnames(two_gene_list)[4] <- "control ability score"

##Calculate regulatory potential score
C1 <- read_csv(".../data/C1_brca.csv")
colnames(C1)[1:2] <- c("Gene", "Score")
colnames(two_gene_list)[1:2] <- c("Gene", "C1")
matched_index <- match(two_gene_list$Gene, C1$Gene)
two_gene_list$C1 <- ifelse(
  is.na(matched_index),
  0,
  C1$Score[matched_index]
)
colnames(two_gene_list)[2] <- "S_Dir"

C2 <- read_csv(".../data/C2_brca.csv")
colnames(C2)[1:2] <- c("Gene", "Score")
matched_index <- match(two_gene_list$Gene, C2$Gene)
two_gene_list$Score2 <- ifelse(
  is.na(matched_index),
  0,
  C2$Score[matched_index]
)
colnames(two_gene_list)[3] <- "S_Ind"

two_gene_list[, 14] <- pmax(two_gene_list[, 2], two_gene_list[, 3])
colnames(two_gene_list)[14] <- "regulatory potential score"

####Step4.Calculate multi-network score####
STRING <- read_csv(".../data/STRING_adj_brca.csv", show_col_types = FALSE)
colnames(STRING)[1:2] <- c("GeneA", "GeneB")
genes <- unique(c(STRING$GeneA, STRING$GeneB))
N <- length(genes)
gene_index <- setNames(1:N, genes)
PPI <- matrix(0, nrow = N, ncol = N)
rownames(PPI) <- genes
colnames(PPI) <- genes
row_ids <- gene_index[STRING$GeneA]
col_ids <- gene_index[STRING$GeneB]
PPI[cbind(row_ids, col_ids)] <- 1
PPI[cbind(col_ids, row_ids)] <- 1

##Calculate TNI
TNI<-matrix(0,nrow(PPI),ncol(PPI))
for(i in 1:7686) 
{
  for (j in 1:7686){ 
    number<-0 
    if(PPI[i,j]==1)
    {
      for (k in 1:7686) { 
        if(k!=j && k!=i && PPI[i,k]==1 && PPI[j,k]==1) 
        {
          number<-number+1
        }
      }
      if(length(which(PPI[i,]==1))>1 && length(which(PPI[j,]==1))>1)
      {
        ECC[i,j]<-number/min(length(which(PPI[i,]==1)),length(which(PPI[j,]==1)))
      }
    }
  }
}
SUMTNI<-matrix(0,nrow(PPI),1)
SUMTNI<-rowSums(TNI) 
SUMTNI_mat <- cbind(Gene = rownames(PPI), Score = as.character(SUMTNI))
colnames(two_gene_list)[1] <- "Gene"
colnames(SUMTNI_mat) <- c("Gene", "Score")
matched_index <- match(two_gene_list$Gene, SUMTNI_mat[, "Gene"])
two_gene_list[,5] <- ifelse(
  is.na(matched_index),
  0,
  as.numeric(SUMTNI_mat[matched_index, "Score"])
)
colnames(two_gene_list)[5] <- "TNI"
maxTNI<- max(two_gene_list[,5])
for(i in 1:nrow(two_gene_list))
{
  two_gene_list[i,6] <- two_gene_list[i,5]/maxTNI    
}
colnames(two_gene_list)[6] <- "NL_TNI"

##Calculate MIE
Mutation<-t(brca_mut)
SM_mut<-Mutation
SM_score=matrix(data=0,nrow=nrow(SM_mut),ncol=2)
SM_score[,1]=rownames(SM_mut)
for(i in 1:nrow(SM_mut))
{
  SM_score[i,2] <- sum(as.numeric(SM_mut[i,]))/ncol(SM_mut)  
}
vertex_W=matrix(0.1,nrow=nrow(Mutation),ncol = ncol(Mutation))
vertexes=rownames(Mutation)
colnames(vertex_W)=colnames(Mutation)
rownames(vertex_W)=rownames(Mutation) 

for(i in 1:ncol(Mutation)){ 
  v_i=intersect(vertexes[which(Mutation[,i]==1)],two_gene_list[,1]) 
  if(length(v_i)!=0){
    for (j in 1:length(v_i)){
      V=v_i[j]
      A=which(SM_score[,1]==V)
      Fdegree=-(as.numeric(SM_score[A,2]))*log2(as.numeric(SM_score[A,2]))
      vertex_W[V,i]=Fdegree
    }
  }
}
vertex_W[is.na(vertex_W)]<-0 
hypergraph=matrix(0,nrow =nrow(vertex_W),ncol=2 )
hypergraph[,1]=rownames(vertex_W)
TT=rowSums(vertex_W)
hypergraph[,2]=TT
two_gene_list[,7] <- 0
match_idx <- match(two_gene_list[,1], hypergraph[,1])
two_gene_list[!is.na(match_idx), 7] <- hypergraph[match_idx[!is.na(match_idx)], 2]
colnames(two_gene_list)[7] <- "MIE"

##Calculate MTCS
two_gene_list[,8] <- as.numeric(as.character(two_gene_list[,6])) *
  as.numeric(as.character(two_gene_list[,7]))
colnames(two_gene_list)[8] <- "MTCS"

##RWR
source(".../code/random_walk.R")
seed <- two_gene_list[, c(1, 8)]
seed <- data.frame(seed)
seed <- seed[order(-seed[,2]), ]

##Select top n genes
top_n <- ceiling(0.05 * nrow(seed))
seed <- seed[1:top_n, 1]

#PPI1：STRING
STRING_result <- random_walk_ranking(PPI, seed)
match_idx <- match(two_gene_list[,1], STRING_result[,1])
two_gene_list[!is.na(match_idx), 9] <- STRING_result[match_idx[!is.na(match_idx)], 2]
colnames(two_gene_list)[9] <- "STRING"

#PPI2:HINT
load(".../data/HINT.RData")
gene_name <- matrix(data=0,nrow=nrow(mul_edge_list),ncol=2)
for (i in 1:nrow(mul_edge_list)) {
  index=mul_edge_list[i,1]
  a=mul_gene_list[index,2]
  gene_name[i,1]=a
  }
for (i in 1:nrow(mul_edge_list)) {
  index2=mul_edge_list[i,2]
  a=mul_gene_list[index2,2]
  gene_name[i,2]=a
 }
genes <- intersect(unique(as.vector(gene_name)), two_gene_list[,1])
g <- graph_from_edgelist(as.matrix(gene_name), directed=FALSE)
g <- induced_subgraph(g, vids = genes)
HINT <- as.matrix(as_adj(g, sparse=FALSE))
HINT <- as.matrix(as_adjacency_matrix(g, sparse=FALSE))  
nonzero_idx <- which(rowSums(HINT) > 0 | colSums(HINT) > 0)
HINT <- HINT[nonzero_idx, nonzero_idx]
HINT_result <- random_walk_ranking(HINT, seed)
match_idx <- match(two_gene_list[,1], HINT_result[,1])
two_gene_list[!is.na(match_idx), 10] <- HINT_result[match_idx[!is.na(match_idx)], 2]
colnames(two_gene_list)[10] <- "HINT"

#PPI3:CPDB
load(".../data/HINT.RData")
CPDB<-PPI
common_genes <- intersect(rownames(CPDB), two_gene_list[,1])
CPDB <- CPDB[common_genes, common_genes]
nonzero_idx <- which(rowSums(CPDB) > 0 | colSums(CPDB) > 0)
CPDB <- CPDB[nonzero_idx, nonzero_idx]
CPDB_result <- random_walk_ranking(CPDB, seeds = seed)
match_idx <- match(two_gene_list[,1], CPDB_result[,1])
two_gene_list[!is.na(match_idx), 11] <- CPDB_result[match_idx[!is.na(match_idx)], 2]
colnames(two_gene_list)[11] <- "CPDB"

#PPI4:MULTINET  
load(".../data/MULTINET.RData")
MULTINET<-PPI
common_genes <- intersect(rownames(MULTINET), two_gene_list[,1])
MULTINET <- MULTINET[common_genes, common_genes]
nonzero_idx <- which(rowSums(MULTINET) > 0 | colSums(MULTINET) > 0)
MULTINET <- MULTINET[nonzero_idx, nonzero_idx]
MULTINET_result <- random_walk_ranking(MULTINET, seeds = seed)
match_idx <- match(two_gene_list[,1], MULTINET_result[,1])
two_gene_list[!is.na(match_idx), 12] <- MULTINET_result[match_idx[!is.na(match_idx)], 2]
colnames(two_gene_list)[12] <- "MULTINET"

##Calculate multi-network diffusion score
score_matrix <- two_gene_list[, c( 9, 10, 11, 12)]
rownames(score_matrix) <- two_gene_list[, 1]  
mean_score <- rowMeans(score_matrix)
min_score <- apply(score_matrix, 1, min)
max_score <- apply(score_matrix, 1, max)
max_score[max_score == 0] <- 1e-6  
#Ci
consistency_score <- min_score / max_score
fused_score <- mean_score * consistency_score
two_gene_list[, 13] <- fused_score
colnames(two_gene_list)[13] <- "multi-network diffusion score"



####Step5.Calculate the final score using SVD####
two_gene_list[,15] <- 0
two_gene_list[,16] <- 0
#Z-score 
zscore <- function(x) {
  if (sd(x) == 0) return(rep(0, length(x)))
  (x - mean(x)) / sd(x)
}
score1_norm <- zscore(two_gene_list[, 4])  
score2_norm <- zscore(two_gene_list[, 13])  
score3_norm <- zscore(two_gene_list[, 14])  
#SVD
S <- cbind(score1_norm, score2_norm, score3_norm)
res <- irlba(S, nv = 1)             
U <- res$u * res$d[1]               
TriForceScore <- abs(U)           
two_gene_list[, 15] <- TriForceScore
colnames(two_gene_list)[15] <- " TriForceScore"
write.csv(two_gene_list, file = ".../code/final_score.csv", row.names = FALSE)


