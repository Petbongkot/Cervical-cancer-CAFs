library(limma)
library(tidyverse)
library(hgu133plus2.db)
library(oligo)
library(pd.hg.u133.plus.2)
library(ggplot2)
library(ggrepel)
library(AnnotationDbi)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(dplyr)
library(cowplot)

#Normalize and Calculate t-test
filename_use <- list.files()
Rawdata <- read.celfiles(filenames=filename_use)
Rawdata<- rma(Rawdata, normalize = F, background = F)
Expression_probe <- exprs(Rawdata)
Table <- data.frame(row.names = filename_use,group = c(rep('normal',10),rep("cancer",21)))

Table <-model.matrix(~0+group,data= Table)
Contrast.mat <- makeContrasts(contrasts = 'groupcancer-groupnormal',levels = colnames(Table))
Fit <- lmFit(Expression_probe,Table)
Fit <- contrasts.fit(Fit,Contrast.mat)
Fit <- eBayes(Fit)

Result_all <-topTable(Fit, n=Inf)
Result_all<-as.data.frame(Result_all)
Result_all <- Result_all %>% mutate(diffexpressed = ifelse(adj.P.Val < 0.05,
                                                           ifelse(logFC < -1.5, 'Down-regulated',
                                                                  ifelse(logFC > 1.5, 'Up-regulated', 'Not significant')),
                                                           'Not significant'))
probe_ids <- rownames(Result_all)
gene_symbols <- mapIds(hgu133plus2.db, keys = probe_ids, keytype = "PROBEID" ,column = "SYMBOL")
gene_symbols <- data.frame(gene_symbols)
Result_all <- Result_all %>% mutate(gene_name = gene_symbols$gene_symbols,.before=logFC)
Result_all <- as.data.frame(Result_all)
Result_all_NA <- Result_all %>% filter_all(any_vars(is.na(.)))
Result_all <-Result_all %>% filter(gene_name != 'NA')
Result_all <- Result_all %>% arrange(desc(logFC))
Dupgene<-Result_all %>% filter(duplicated(gene_name)) %>% as.data.frame()
Result_all <- Result_all %>% distinct(gene_name, .keep_all = T)


#Summary
Result_all %>% mutate(diffexpressed = diffexpressed %>% as.character())
upregulated <- Result_all[which(Result_all$diffexpressed == 'Up-regulated'),]
downregulated <- Result_all[which(Result_all$diffexpressed == 'Down-regulated'),]
upregulated
downregulated

topDEG <- c(head( Result_all[order(dplyr::desc(Result_all$logFC)),], n = 30)$gene_name,
            tail( Result_all[order(dplyr::desc(Result_all$logFC)),], n = 30)$gene_name)
Result_all$topDEG <- ifelse(Result_all$gene_name %in% topDEG,Result_all$gene_name, NA)

volcano_plot <- Result_all %>%
  ggplot(aes(x = logFC, y = -log10(adj.P.Val), col= diffexpressed,label= topDEG)) +
  geom_vline(xintercept = c(-1.5, 1.5), col = "black", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "black", linetype = 'dashed') +
  geom_point(shape =19, size =1,) +
  scale_color_manual(values = c("blue", "grey", "red"),
                     labels = c("Down-regulated", "Not significant", "Up-regulated"))+
  coord_cartesian(ylim = c(0, 15), xlim = c(-5, 5)) +
  labs(x = expression("logFC"), y = expression("-log"[10]*"adj.P.Val")) +
  scale_x_continuous(breaks = seq(-10,10,2))+
  geom_text_repel(size = 2,max.overlaps = 15)+
  theme_minimal()+
  theme(legend.position = "inside",legend.position.inside = c(0.8, 0.9), plot.title = element_text(hjust = 0.5)) 

volcano_plot <- ggdraw() +
  draw_plot(volcano_plot, x = 0, y = 0, width = 1, height = 1) +
  draw_label(
    "B",
    x = 0.01, y = 0.99,
    hjust = 0, vjust = 1,
    fontface = "bold",
    size = 18
  )

volcano_plot
