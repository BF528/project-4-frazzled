#install.packages("enrichR")
library(enrichR)
library(dplyr)
library(ggplot2)

dbs <- c("GO_Molecular_Function_2015", "GO_Cellular_Component_2015", "GO_Biological_Process_2015")
#sample_data = read.csv('/projectnb/bf528/project_4_scrnaseq/GSM2230760_marker_genes.csv')
real_data = read.csv('/projectnb/bf528/users/frazzled/project_4/analyst/Marker_Genes')
real_data <- real_data[real_data$avg_log2FC > 0.5,]
by_cluster<-split(real_data, real_data$cluster)

output_list = list()
i = 0
#enrich
for (table in by_cluster) {
  enriched <- enrichr(table$gene, databases = dbs)
  enriched[[1]]['Source'] = 'Molecular Function'
  enriched[[2]]['Source'] = 'Cellular Component'
  enriched[[3]]['Source'] = 'Biological Process'
  all_data = data.frame(rbind(enriched[[1]], enriched[[2]], enriched[[3]]))
  output_list[[toString(i)]] <-all_data %>% slice_min(Adjusted.P.value, n = 10)
  i = i+1
}

#
i = 1
for (table in output_list) {
  print(i)
  df = output_list[[i]]
  df$num_genes = lengths(strsplit(df$Genes, ';'))
  print(df$num_genes)
  print(ggplot(df, aes(y=Term, x=num_genes, fill = Adjusted.P.value)) + 
    geom_bar(stat = "identity") + scale_color_gradient(low="blue", high="red" ) 
    + ggtitle(toString(names(by_cluster[i]))))

  i = i +1
}
#filtering by p value
real_data <- real_data[real_data$p_val_adj < 0.05,]
by_cluster<-split(real_data, real_data$cluster)

output_list = list()
i = 0
#enrich
for (table in by_cluster) {
  enriched <- enrichr(table$gene, databases = dbs)
  enriched[[1]]['Source'] = 'Molecular Function'
  enriched[[2]]['Source'] = 'Cellular Component'
  enriched[[3]]['Source'] = 'Biological Process'
  all_data = data.frame(rbind(enriched[[1]], enriched[[2]], enriched[[3]]))
  output_list[[toString(i)]] <-all_data %>% slice_min(Adjusted.P.value, n = 30)
  i = i+1
}

#
i = 1
for (table in output_list) {
  print(i)
  df = output_list[[i]]
  df$num_genes = lengths(strsplit(df$Genes, ';'))
  print(df$num_genes)
  print(ggplot(df, aes(y=Term, x=num_genes, fill = Adjusted.P.value)) + 
          geom_bar(stat = "identity") + scale_color_gradient(low="blue", high="red" ) 
        + ggtitle(toString(names(by_cluster[i]))))
  
  i = i +1
}


