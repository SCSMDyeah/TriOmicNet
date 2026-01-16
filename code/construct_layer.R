construct_layer <- function(G_mut,network,protein_info,DE_Score){
  DEGs <- DE_Score
  edges <- network
  gene_point <- protein_info
  gene_first <- colnames(G_mut) 
  gene_third <- setdiff(DEGs[,1],gene_first)
  gene_second <- setdiff(DEGs[,1],gene_third)
  gene_point1 <- gene_point[match(gene_third,gene_point$preferred_name),1:2]
  dele_row <- which(is.na(gene_point1[,2])==TRUE)
  gene_third <- setdiff(gene_third,gene_third[dele_row])
  gene_point1 <- gene_point[match(gene_third,gene_point$preferred_name),1:2]
  n1 <- 0
  setdiff_gene <- matrix(0,1,ncol = 10000)
  for (i in 1:length(gene_third)) {
    protein_name <-  gene_point1[i,1]
    protein_list <- which(edges[,2]==protein_name)
    if(length(protein_list)==0){
      n1 <- n1 + 1
      setdiff_gene[1,n1] <- gene_third[i]
    }
  }
  setdiff_gene <- setdiff(setdiff_gene[1,],"0")
  gene_third <- setdiff(gene_third,setdiff_gene)
  gene_point1 <- gene_point[match(gene_third,gene_point$preferred_name),1:2]
  for (i in 1:length(gene_third)) {
    protein_name <-  gene_point1[i,1]
    protein_list <- which(edges[,2]==protein_name)
    edges[protein_list,2] <- gene_third[i]
  }
 
  edges1 <- edges[order(edges[,2],decreasing = T),]
  aa <- which(edges1[,2]>"A") 
  edges_third <- edges1[aa,]
  layer34_edges <- edges_third

  for (i in 1:length(gene_second)) {
    gene_name <- gene_second[i]
    p_row <- which(gene_point[,2]==gene_name)
    p_name <- gene_point[p_row,1]
    row_list <- which(layer34_edges[,1]==p_name)
    if(length(row_list)==0)next
    layer34_edges[row_list,1] <- gene_name
  }

  aa <- which(layer34_edges[,1]>"A")
  layer34_edges <- layer34_edges[aa,]

  gene_second <- unique(layer34_edges[,1])
  gene_score_second <- matrix(0,nrow = length(gene_second),ncol = 2)
  gene_score_second[,1] <- gene_second
  for (i in 1:nrow(gene_score_second)) {
    gene1 <- gene_score_second[i,1]
    row_edge_list <- which(layer34_edges[,1]==gene1)
    down_gene_list <- layer34_edges[row_edge_list,2]
    exp_score <- 0
    for (j in 1:length(down_gene_list)) {
      down_exp_index <- which(DEGs[,1]==down_gene_list[j])
      exp_score <- exp_score + as.numeric(DEGs[down_exp_index,2])
    }
    
    gene_score_second[i,2] <- exp_score
  }

  gene_first <- setdiff(gene_first,DEGs[,1])
  layer24_edges <- edges_third
  for (i in 1:length(gene_first)) {
    gene_name <- gene_first[i]
    p_row <- which(gene_point[,2]==gene_name)
    p_name <- gene_point[p_row,1]
    row_list <- which(layer24_edges[,1]==p_name)
    if(length(row_list)==0)next
    layer24_edges[row_list,1] <- gene_name
  }
  aa <- which(layer24_edges[,1]>"A")
  layer24_edges <- layer24_edges[aa,]

  gene_first1 <- unique(layer24_edges[,1])
  gene_score3 <- matrix(0,nrow = length(gene_first1),ncol = 2)
  gene_score3[,1] <- gene_first1
  for (i in 1:nrow(gene_score3)) {
    gene1 <- gene_score3[i,1]
    row_edge_list <- which(layer24_edges[,1]==gene1)
    down_gene_list <- layer24_edges[row_edge_list,2]
    exp_score <- 0
    for (j in 1:length(down_gene_list)) {
      down_exp_index <- which(DEGs[,1]==down_gene_list[j])
      exp_score <- exp_score + as.numeric(DEGs[down_exp_index,2])
    }
    gene_score3[i,2] <- exp_score
  }
  write.csv(gene_score3,".../data/C1_brca.csv",row.names = F)

  edges <- network
  gene_point1 <- gene_point[match(gene_second,gene_point$preferred_name),1:2]
  for (i in 1:length(gene_second)) {
    protein_name <-  gene_point1[i,1]
    protein_list <- which(edges[,2]==protein_name)
    edges[protein_list,2] <- gene_second[i]
  }
  aa <- which(edges[,2]>"A")
  edges_second <- edges[aa,]
  for (i in 1:length(gene_first)) {
    gene_name <- gene_first[i]
    p_row <- which(gene_point[,2]==gene_name)
    p_name <- gene_point[p_row,1]
    row_list <- which(edges_second[,1]==p_name)
    if(length(row_list)==0)next
    edges_second[row_list,1] <- gene_name
  }
  aa <- which(edges_second[,1]>"A")
  layer23_edges <- edges_second[aa,]

  gene_firs2 <- unique(layer23_edges[,1])
  gene_score2 <- matrix(0,nrow = length(gene_firs2),ncol = 2)
  gene_score2[,1] <- gene_firs2
  for (i in 1:nrow(gene_score2)) {
    gene1 <- gene_score2[i,1]
    row_edge_list <- which(layer23_edges[,1]==gene1)
    down_gene_list <- layer23_edges[row_edge_list,2]
    exp_score <- 0
    for (j in 1:length(down_gene_list)) {
      down_exp_index <- which(gene_score_second[,1]==down_gene_list[j])
      exp_score <- exp_score + as.numeric(gene_score_second[down_exp_index,2])
    }
    gene_score2[i,2] <- exp_score/length(down_gene_list)
  }
  write.csv(gene_score2,".../data/C2_brca.csv",row.names = F)

  net_gene <- union(gene_score2[,1],gene_score3[,1])
  write.csv(net_gene,".../data/layer_2_brca.csv",row.names = F)
  edges <- network
  gene_point1 <- gene_point[match(net_gene,gene_point$preferred_name),1:2]
  for (i in 1:length(net_gene)) {
    protein_name <-  gene_point1[i,1]
    protein_list <- which(edges[,2]==protein_name)
    edges[protein_list,2] <- net_gene[i]
  }
  aa <- which(edges[,2]>"A")
  net_edges <- edges[aa,]
  setdiff_gene <- matrix(0,1,ncol = 10000)
  k <- 0
  for (i in 1:length(net_gene)) {
    protein_name <-  gene_point1[i,1]
    protein_list <- which(net_edges[,1]==protein_name)
    if(length(protein_list)==0){
      k = k + 1
      setdiff_gene[1,k] <- i
    }
  }
  setdiff_gene <- setdiff(setdiff_gene[1,],0)
  net_gene <- setdiff(net_gene,net_gene[setdiff_gene])
  write.csv(net_gene,".../data/layer2_gene_list_brca.csv",row.names = F)

  gene_point1 <- gene_point[match(net_gene,gene_point$preferred_name),1:2]
  edges <- network
  for (i in 1:length(net_gene)) {
    protein_name <-  gene_point1[i,1]
    protein_list <- which(edges[,2]==protein_name)
    edges[protein_list,2] <- i
  }
  aa <- which(edges[,2]>=9607)
  bb <- which(edges[,2]<=9606)
  net_edges <- rbind(edges[aa,],edges[bb,])
  for (i in 1:length(net_gene)) {
    protein_name <-  gene_point1[i,1]
    protein_list <- which(net_edges[,1]==protein_name)
    net_edges[protein_list,1] <- i
  }
  aa <- which(net_edges[,1]>=9607)
  bb <- which(net_edges[,1]<=9606)
  net_edges <- rbind(net_edges[aa,],net_edges[bb,])
  write.csv(net_edges,".../data/STRING_PPI_brca.csv",row.names = F)

  net_edges <- read.csv(".../data/STRING_PPI_brca.csv", header = TRUE)
  net_gene <- read.csv(".../data/layer2_gene_list_brca.csv", header = TRUE)

  if (is.data.frame(net_gene)) {
    net_gene <- net_gene[,1]
  }

  id_to_gene <- data.frame(id = 1:length(net_gene), gene = net_gene, stringsAsFactors = FALSE)

  net_edges_named <- net_edges
  net_edges_named[,1] <- id_to_gene$gene[match(net_edges[,1], id_to_gene$id)]
  net_edges_named[,2] <- id_to_gene$gene[match(net_edges[,2], id_to_gene$id)]

  write.csv(net_edges_named, ".../data/STRING_adj_brca.csv", row.names = FALSE)
}