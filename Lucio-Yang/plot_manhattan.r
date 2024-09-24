rm(list=ls())
args <- commandArgs(T)
library(ggplot2)
library(ggforce)
library(tidyverse)
library(ggrepel)
library(gridExtra)
assoc_file <- args[1]
output_file <- paste(assoc_file,'.png',sep = '')
data <- read.table(assoc_file, header = T, sep = '\t')
plot_gwas <- function(data,output){
  chr_len <- data %>% group_by(chr) %>% summarise(chr_len=max(ps))
  chr_pos <- chr_len  %>% mutate(total = cumsum(as.numeric(chr_len)) - chr_len) %>%select(-chr_len)
  Snp_pos <- chr_pos %>%left_join(data, ., by="chr") %>%arrange(chr, ps) %>%mutate( BPcum = ps + total)
  X_axis <-  Snp_pos %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  data <- Snp_pos %>%
    #mutate( is_highlight=ifelse(SNP %in% snpsOfInterest, "yes", "no")) %>%
    mutate(is_annotate=ifelse(-log10(p_wald)>7, "yes", "no"))

  p <- ggplot(data, aes(x=BPcum, y=-log10(p_wald))) +
  geom_point( aes(color=as.factor(chr)), alpha=0.8, size=2.5) +
  #scale_color_manual(values = rep(c("blue4", "orange"), 14 )) +
  #scale_color_manual(values = rep(c("#C4493A", "#90B07A"), 21 )) +
  scale_color_manual(values = c("#ff0000","#fc4444","#fc6404","#fcd444","#8cc43c",
                                "#029658","#1abc9c","#5bc0de","#6454ac","#fc8c84",
                                "#ff0000","#fc6404","#8cc43c","#1abc9c"))+
  scale_x_continuous( label =X_axis$chr, breaks= X_axis$center )+
  scale_y_continuous(expand = c(0, 0) ) +
  coord_cartesian(ylim = c(0,-log10(min(data$p_wald))+5),expand = FALSE)+
  # add threshold line
  geom_hline(yintercept = c(-log10(1/nrow(Snp_pos)), -log10(0.05/nrow(Snp_pos))), color = c('green', 'red'), size = 1.2, linetype = c("twodash", "twodash"))+ 
  # spot if highlight
  #geom_point(data=subset(data, is_highlight=="yes"), color="orange", size=2) +
  # higlight label
  #geom_label_repel( data=subset(data, is_annotate=="yes"), aes(label=SNP), size=2) +
  # facet_zoom(x = BPcum >= 8828782866 & BPcum <=8840380615)+
  #theme_bw() +
  xlab('Chromosome')+ ylab('-log10(P-value)')+labs(title = '')+
  theme(panel.background = element_blank(),
        legend.position="none",
        axis.line = element_line(colour = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 16,color = 'black'),
        axis.text.x = element_text(size = 16,color = 'black'),
        axis.title.y.left = element_text(size = 18),
        axis.title.x.bottom = element_text(size = 18),
        title = element_text(size=18))
  ggsave(dpi = 300,width = 14,height = 8,filename = output,plot = p)
}
plot_gwas(data,output_file)
