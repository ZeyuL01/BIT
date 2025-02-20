library(patchwork)

work_dir<-"/Users/zeyulu/Desktop/Project/BIT/revision_data/comparison/BIT/"
work_files_BIT<-list.files(work_dir,pattern="*_rank_table.csv")

work_files_BIT<-work_files_BIT[c(1,8,9,2,3,4,5,6,7)]
work_dir
work_files_BIT
#PART B
cell_names<-sapply(strsplit(work_files_BIT,"_rank",fixed=TRUE),function(x){return(x[[1]])})
cell_names<-tools::toTitleCase(cell_names)
cell_names<-c("B cell","Monocyte","NK cell","CD4+ T cell","CD8+ T cell","Dendritic cell","gdT cell","HSPC","MAIT cell")
colors<-c("#9B3A4D","#FC8002","#394A92","#70A0AC","#D2352C","#8E549E","#BAAFD1","#497EB2","#ADDB88")
plot_list<-list()
for(i in 1:9){
  table<-read.csv(paste0(work_dir,work_files_BIT[i]),row.names=1)
  data<-data.frame(TR=table$TR[1:10],
                   Score=table$BIT_score[1:10],
                   Lower=table$BIT_score_lower[1:10],
                   Upper=table$BIT_score_upper[1:10],
                   group=cell_names[i])
  data$TR<-factor(data$TR,levels=rev(data$TR))
  plot_list[[i]]<-ggplot(data, aes(x = TR, y = Score)) +
    geom_bar(stat = "identity", fill = colors[i], color = "black") +
    geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.4) +
    coord_flip() +
    theme_bw() + scale_y_continuous(limits = c(0, max(data$Upper)*1.2),breaks = scales::breaks_extended(n = 4),expand = c(0,0)) + facet_grid(.~group)+
    labs(
      title = "",
      x = "Transcriptional Regulator",
      y = "BIT Score"
    )+theme(title=element_text(size=9),plot.margin=unit(c(0.25,0.25,0.25,0.25),"cm"),axis.text.y = element_text(size=10,color="black"),
            axis.text.x = element_text(size=10,color="black"),
            axis.title.x = element_text(size=14,color="black"),
            axis.title.y = element_text(size=14,color="black"),
            strip.background = element_rect(fill="#DBD1B6"),
            strip.text = element_text(size=12, colour="black",margin=ggplot2::margin(1,1,1,1,"mm")))
}

p_combined<-plot_list[[1]]+plot_list[[2]]+ plot_list[[3]] +
  plot_list[[4]]+plot_list[[5]]+plot_list[[6]] +
  plot_list[[7]]+plot_list[[8]]+plot_list[[9]]+ plot_layout(ncol=3,guides="collect",axis_titles = "collect")
print(p_combined)



####################
library(clusterProfiler)

split_go_term <- function(term) {
  words <- strsplit(term, " ")[[1]]
  n_words <- length(words)

  # Find split points for 2-3 lines
  if(n_words <= 3) {
    split_at <- ceiling(n_words/2)
  } else {
    split_at <- c(ceiling(n_words/3), ceiling(n_words*2/3))
  }

  # Split words into groups
  groups <- split(words, cut(seq_along(words), breaks = c(0, split_at, n_words)))

  # Combine with newlines
  paste(sapply(groups, paste, collapse = " "), collapse = "\n")
}

integer_breaks <- function(x, n = 4) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(numeric())
  rng <- range(x)
  breaks <- floor(rng[1]):ceiling(rng[2])
  breaks <- unique(round(breaks))
  if (length(breaks) > n) {
    breaks <- breaks[seq(1, length(breaks), length.out = n)]
  }
  breaks
}

work_files_BIT
plot_list<-list()

for(i in 1:9){
  table<-read.csv(paste0(work_dir,work_files_BIT[i]),row.names=1)
  data<-data.frame(TR=table$TR[1:20],
                   group=cell_names[i])
  BIT_gene_ids<-bitr(data$TR,fromType = "SYMBOL",toType = "ENTREZID",OrgDb="org.Hs.eg.db")
  BIT_Results=enrichGO(gene          = BIT_gene_ids$ENTREZID[1:20],
                       OrgDb = "org.Hs.eg.db",
                       ont = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.01,
                       qvalueCutoff  = 0.05)
  GO_BIT_table<-head(BIT_Results,5)
  GO_PLOT_Table_BIT<-data.frame(GO=GO_BIT_table$Description,
                                GeneRatio=Trans_to_double(GO_BIT_table),
                                Pvalue=GO_BIT_table$pvalue,
                                Count=GO_BIT_table$Count)
  GO_PLOT_Table_BIT$GO <- sapply(GO_PLOT_Table_BIT$GO, split_go_term)
  plot_list[[i]] <- ggplot(GO_PLOT_Table_BIT, aes(x = GeneRatio, y = reorder(GO, -Pvalue), size = Count, color = Pvalue)) +
    geom_point() +
    # Fix Pvalue color legend (4 breaks, no overlap)
    scale_color_gradient(
      low = "red",
      high = "blue",
      limits = c(min(GO_BIT_table$pvalue), max(GO_BIT_table$pvalue)),
      breaks = scales::breaks_pretty(n = 4), # 4 breaks # 2 decimal places
      guide = guide_colorbar(
        order = 1,
        title.position = "top",
        barheight = unit(1.6, "cm"),
        ticks = FALSE
      )
    ) +
    # Fix Count size legend (integer breaks)
    scale_size_continuous(
      breaks = integer_breaks(GO_PLOT_Table_BIT$Count, n = 4), # 4 integer breaks
      range = c(2, 5), # Adjust point sizes
      guide = guide_legend(
        order = 2,
        title.position = "top",
        override.aes = list(color = "black")
      )
    ) +
    theme_bw() +
    labs(y = "GO", x = "Gene Ratio", color = "P-Value", size = "Count") +
    ggtitle(cell_names[i]) +
    theme(
      text = element_text(size = 12),
      legend.position = "right",
      legend.box = "vertical",
      legend.spacing.y = unit(0.1, "cm"),
      legend.margin = margin(0, 0, 0, 0),
      axis.text.y = element_text(color = "black")
    ) +
    xlim(c(0.0, 0.6))
}
p_combined<-plot_list[[1]]+plot_list[[2]]+ plot_list[[3]] +
plot_list[[4]]+plot_list[[5]]+plot_list[[6]] +
  plot_list[[7]]+plot_list[[8]]+plot_list[[9]]+ plot_layout(ncol=3,axis_titles = "collect")

print(p_combined)

GO_PLOT_Table_BIT

work_files_BIT

##############BcellComp
work_dir<-"/Users/zeyulu/Desktop/Project/BIT/revision_data/comparison/"
work_dir_HPA<-"/Users/zeyulu/Desktop/Project/BIT/revision_data/comparison/HPA/"

HPA_TFs<-c("B_cell_TF.tsv","dendritic_TF.tsv","monocytes_TF.tsv","NK_cells_TF.tsv","T_cells_TF.tsv")
Cell_Type_names<-c("B cell","Dendritic cell","Monocyte","NK cell","T cell")
HPA_TFs_df_list<-list()
for(i in 1:5){
  test_df<-read.table(paste0(work_dir_HPA,HPA_TFs[i]),sep="\t",header=TRUE)
  HPA_TFs_df_list[[Cell_Type_names[i]]]<-test_df$Gene
}

BIT_result<-read.csv(paste0(work_dir,"BIT.csv"),row.names=1)
ArchR_result<-read.csv(paste0(work_dir,"ArchR.csv"),row.names=1)
for(i in 1:9){
  ArchR_result[,i]<-sapply(strsplit(ArchR_result[,i],"_",fixed=TRUE),function(x){return(x[[1]])})
}
scbasset_result<-read.csv(paste0(work_dir,"scbasset.csv"))
scbasset_result<-scbasset_result[,c(seq(1,18,2))]
SnapATAC2_result<-read.csv(paste0(work_dir,"SnapATAC2.csv"),row.names=1)

Cell_Types<-c("B cell","CD4+ cell","CD8+ cell","Dendritic cell","gdT cell","HSPC","MAIT cell","Monocyte","NK cell")

colnames(BIT_result)<-Cell_Types
colnames(ArchR_result)<-Cell_Types
colnames(scbasset_result)<-Cell_Types
colnames(SnapATAC2_result)<-Cell_Types

Top20_list<-list("BIT"=BIT_result$`B cell`[1:20],
"ArchR"=ArchR_result$`B cell`[1:20],
"scBasset"=scbasset_result$`B cell`[1:20],
"SnapATAC2"=SnapATAC2_result$`B cell`[1:20])

Top20_list

for(i in 1:4){
  BIT_gene_ids<-bitr(Top20_list[[i]],fromType = "SYMBOL",toType = "ENTREZID",OrgDb="org.Hs.eg.db")
  BIT_Results=enrichGO(gene          = BIT_gene_ids$ENTREZID[1:20],
                       OrgDb = "org.Hs.eg.db",
                       ont = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.01,
                       qvalueCutoff  = 0.05)
  GO_BIT_table<-head(BIT_Results,5)
  GO_PLOT_Table_BIT<-data.frame(GO=GO_BIT_table$Description,
                                GeneRatio=Trans_to_double(GO_BIT_table),
                                Pvalue=GO_BIT_table$pvalue,
                                Count=GO_BIT_table$Count)
  plot_list[[i]] <- ggplot(GO_PLOT_Table_BIT, aes(x = GeneRatio, y = reorder(GO, -Pvalue), size = Count, color = Pvalue)) +
    geom_point() +
    # Fix Pvalue color legend (4 breaks, no overlap)
    scale_color_gradient(
      low = "red",
      high = "blue",
      limits = c(min(GO_BIT_table$pvalue), max(GO_BIT_table$pvalue)),
      breaks = scales::breaks_pretty(n = 4), # 4 breaks # 2 decimal places
      guide = guide_colorbar(
        order = 1,
        title.position = "top",
        barheight = unit(1.6, "cm"),
        ticks = FALSE
      )
    ) +
    # Fix Count size legend (integer breaks)
    scale_size_continuous(
      breaks = integer_breaks(GO_PLOT_Table_BIT$Count, n = 4), # 4 integer breaks
      range = c(2, 5), # Adjust point sizes
      guide = guide_legend(
        order = 2,
        title.position = "top",
        override.aes = list(color = "black")
      )
    ) +
    theme_bw() +
    labs(y = "GO", x = "Gene Ratio", color = "P-Value", size = "Count") +
    theme(
      text = element_text(size = 12),
      legend.position = "right",
      legend.box = "vertical",
      legend.spacing.y = unit(0.1, "cm"),
      legend.margin = margin(0, 0, 0, 0),
      axis.text.y = element_text(color = "black")
    ) +
    xlim(c(0.0, 0.6))
}
p_combined<-plot_list[[1]]+plot_list[[2]]+ plot_list[[3]] +
  plot_list[[4]]+ plot_layout(ncol=2,axis_titles = "collect")

print(p_combined)



BIT_gene_ids<-bitr(Top20_list[[1]],fromType = "SYMBOL",toType = "ENTREZID",OrgDb="org.Hs.eg.db")
BIT_Results=enrichGO(gene          = BIT_gene_ids$ENTREZID[1:20],
                     OrgDb = "org.Hs.eg.db",
                     ont = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05)
GO_BIT_table<-head(BIT_Results,5)
GO_PLOT_Table_BIT<-data.frame(GO=GO_BIT_table$Description,
                              GeneRatio=Trans_to_double(GO_BIT_table),
                              Pvalue=GO_BIT_table$pvalue,
                              Count=GO_BIT_table$Count)
GO_PLOT_Table_BIT$GO[5]<-"adaptive immune response based on\nsomatic recombination of immune receptors\nbuilt from immunoglobulin superfamily domains"
plot_list[[1]] <- ggplot(GO_PLOT_Table_BIT, aes(x = GeneRatio, y = reorder(GO, -Pvalue), size = Count, color = Pvalue)) +
  geom_point() +
  # Fix Pvalue color legend (4 breaks, no overlap)
  scale_color_gradient(
    low = "red",
    high = "blue",
    limits = c(min(GO_BIT_table$pvalue), max(GO_BIT_table$pvalue)),
    breaks = scales::breaks_pretty(n = 4), # 4 breaks # 2 decimal places
    guide = guide_colorbar(
      order = 1,
      title.position = "top",
      barheight = unit(1.6, "cm"),
      ticks = FALSE
    )
  ) +
  # Fix Count size legend (integer breaks)
  scale_size_continuous(
    breaks = integer_breaks(GO_PLOT_Table_BIT$Count, n = 4), # 4 integer breaks
    range = c(2, 5), # Adjust point sizes
    guide = guide_legend(
      order = 2,
      title.position = "top",
      override.aes = list(color = "black")
    )
  ) +
  theme_bw() +
  labs(y = "GO", x = "Gene Ratio", color = "P-Value", size = "Count") +
  theme(
    text = element_text(size = 12),
    legend.position = "right",
    legend.box = "vertical",
    legend.spacing.y = unit(0.1, "cm"),
    legend.margin = margin(0, 0, 0, 0),
    axis.text.y = element_text(size = 11,color = "black"),
    axis.text.x = element_text(size = 11,color = "black")
  ) +
  xlim(c(0.0, 0.6))


plot_list[[4]]
