library(cancereffectsizeR)
library(ggplot2)
library(stringr)
library(scales)

scientific <- function(x){
  ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", label_scientific()(x)))))
}

#Load in your epistasis data and change the object name
load("PRAD_epistasis_final.RData")
PRAD_epistasis <- epistasis_results(PRAD_epistasis)
PRAD_epistasis <- PRAD_epistasis$gene_epistasis_1

###SPOP###

#Select your gene pairs of interest
SPOP_list <- c(15, 31, 46, 60, 73, 85, 96, 106, 115, 123, 130, 136, 141, 145, 148, 151, 152)

epistatic_change_SPOP <- c()

#Decoupling the gene pairs
for(x in SPOP_list){
  gene1_after_gene2 <- unlist(c(as.character("gene1_after_gene2"), as.numeric(PRAD_epistasis[x,5] - PRAD_epistasis[x,3])))
  gene2_after_gene1 <- unlist(c(as.character("gene2_after_gene1"), as.numeric(PRAD_epistasis[x,6] - PRAD_epistasis[x,4])))
  gene1_after_gene2[1] <- str_replace(gene1_after_gene2[1], "gene1", as.character(PRAD_epistasis[x,1]))
  gene1_after_gene2[1] <- str_replace(gene1_after_gene2[1], "gene2", as.character(PRAD_epistasis[x,2]))
  gene2_after_gene1[1] <- str_replace(gene2_after_gene1[1], "gene1", as.character(PRAD_epistasis[x,1]))
  gene2_after_gene1[1] <- str_replace(gene2_after_gene1[1], "gene2", as.character(PRAD_epistasis[x,2]))
  epistatic_change_SPOP <- rbind(epistatic_change_SPOP, gene1_after_gene2, gene2_after_gene1)
}

epistatic_change_SPOP <- data.frame(gene = epistatic_change_SPOP[,1], change = as.numeric(epistatic_change_SPOP[,2]))

#Separating "Before" and "After"
epistatic_change_SPOP_before <- epistatic_change_SPOP[grep("SPOP_", epistatic_change_SPOP[,1]),]
epistatic_change_SPOP_before$time <- rep("Before", 17)
epistatic_change_SPOP_before <- epistatic_change_SPOP_before[order(-epistatic_change_SPOP_before$change),]
epistatic_change_SPOP_after <- epistatic_change_SPOP[grep("_SPOP", epistatic_change_SPOP[,1]),]
epistatic_change_SPOP_after$time <- rep("After", 17)
epistatic_change_SPOP_after <- epistatic_change_SPOP_after[order(-epistatic_change_SPOP_after$change),]

#Need to have extra underscore to have unique names
epistatic_change_SPOP_before[,1] <- sub("SPOP_after_", "", epistatic_change_SPOP_before[,1])
epistatic_change_SPOP_after[,1] <- sub("after_SPOP", "", epistatic_change_SPOP_after[,1])

#Blank spot
blank <- data.frame(gene = "BLANK", change = -5000, time = "Before")

epistatic_change_SPOP <- rbind(epistatic_change_SPOP_before, blank, epistatic_change_SPOP_after)

gene_labels <- c(epistatic_change_SPOP_after$gene, "", epistatic_change_SPOP_before$gene)
gene_labels <- sub("_", "", gene_labels)

waterfall_SPOP <- ggplot(epistatic_change_SPOP, aes(x= reorder(gene, change), y=change, fill=time)) +
  geom_bar(stat = "identity", position = "dodge") + theme_classic() +
  scale_fill_manual(values = c("#F8766D","#00BFC4"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "bottom", legend.title = element_blank()) +
  ggtitle("SPOP gene pairs") + xlab("Gene") + ylab ("Epistatic change in selection")+
  scale_x_discrete(labels = gene_labels) + scale_fill_discrete(breaks = c("Before", "After"))+
  scale_y_continuous(labels = scientific, breaks = c(-1e4, 0, 1e4), limits = c(-2e4, 2e4))

waterfall_SPOP

###TP53###

#Select your gene pairs of interest
TP53_list <- c(33, 87, 98, 125, 132, 138, 152, 153)

epistatic_change_TP53 <- c()

#Decoupling the gene pairs
for(x in TP53_list){
  gene1_after_gene2 <- unlist(c(as.character("gene1_after_gene2"), as.numeric(PRAD_epistasis[x,5] - PRAD_epistasis[x,3])))
  gene2_after_gene1 <- unlist(c(as.character("gene2_after_gene1"), as.numeric(PRAD_epistasis[x,6] - PRAD_epistasis[x,4])))
  gene1_after_gene2[1] <- str_replace(gene1_after_gene2[1], "gene1", as.character(PRAD_epistasis[x,1]))
  gene1_after_gene2[1] <- str_replace(gene1_after_gene2[1], "gene2", as.character(PRAD_epistasis[x,2]))
  gene2_after_gene1[1] <- str_replace(gene2_after_gene1[1], "gene1", as.character(PRAD_epistasis[x,1]))
  gene2_after_gene1[1] <- str_replace(gene2_after_gene1[1], "gene2", as.character(PRAD_epistasis[x,2]))
  epistatic_change_TP53 <- rbind(epistatic_change_TP53, gene1_after_gene2, gene2_after_gene1)
}

epistatic_change_TP53 <- data.frame(gene = epistatic_change_TP53[,1], change = as.numeric(epistatic_change_TP53[,2]))

#Separating "Before" and "After"
epistatic_change_TP53_before <- epistatic_change_TP53[grep("TP53_", epistatic_change_TP53[,1]),]
epistatic_change_TP53_before$time <- rep("Before", 8)
epistatic_change_TP53_before <- epistatic_change_TP53_before[order(-epistatic_change_TP53_before$change),]
epistatic_change_TP53_after <- epistatic_change_TP53[grep("_TP53", epistatic_change_TP53[,1]),]
epistatic_change_TP53_after$time <- rep("After", 8)
epistatic_change_TP53_after <- epistatic_change_TP53_after[order(-epistatic_change_TP53_after$change),]

#Need to have extra underscore to have unique names
epistatic_change_TP53_before[,1] <- sub("TP53_after_", "", epistatic_change_TP53_before[,1])
epistatic_change_TP53_after[,1] <- sub("after_TP53", "", epistatic_change_TP53_after[,1])

#Blank spot
blank <- data.frame(gene = "BLANK", change = -5000, time = "Before")

epistatic_change_TP53 <- rbind(epistatic_change_TP53_before, epistatic_change_TP53_after)

gene_labels <- c("SPOP", "CTNNB1", "KMT2D", "AR", "FOXA1", "MUC16", "KMT2D", "SYNE1", "KMT2C",
                 "FOXA1", "SPOP", "CTNNB1", "SYNE1", "AR", "KMT2C", "MUC16")

waterfall_TP53 <- ggplot(epistatic_change_TP53, aes(x= reorder(gene, change), y=change, fill=time)) +
  geom_bar(stat = "identity", position = "dodge") + theme_classic() +
  scale_fill_manual(values = c("#F8766D","#00BFC4"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "bottom", legend.title = element_blank()) +
  ggtitle("TP53 gene pairs") + xlab("Gene") + ylab ("Epistatic change in selection")+
  scale_x_discrete(labels = gene_labels) + scale_fill_discrete(breaks = c("Before", "After"))+
  scale_y_continuous(labels = scientific, breaks = c(-1e4, 0, 1e4), limits = c(-2e4, 2e4))

waterfall_TP53
