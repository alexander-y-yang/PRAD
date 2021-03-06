library(cancereffectsizeR)
library(scales)
library(stringr)
library(dplyr)
library(ggplot2)
library(cowplot)

scientific <- function(x){
  ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", label_scientific()(x)))))
}


#Load in your file
threestage_final <- load_cesa("threestage_final.rds")
threestage_results <- snv_results(threestage_final)
threestage_results <- threestage_results$selection.1

threestage_results$variant_name <- str_replace(threestage_results$variant_name, "_", " ")

aac <- threestage_results$variant_type == "aac"
threestage_results <- threestage_results[aac,]
threestage_results <- threestage_results[,c(1,3:5,7:12,19,37,39,41)]

#Separating the data in to early/late/met

threestage_results_early <- threestage_results[,c(1:2,5:6,11,12)]
all <-threestage_results_early$maf_freq_in_Early >=1
threestage_results_early <- threestage_results_early[all,]
threestage_results_early <- threestage_results_early[order(-si_1),]
colnames(threestage_results_early)[2] <- "si"
threestage_results_early$progression <- rep("Early", length(threestage_results_early$variant_name))

threestage_results_late <- threestage_results[,c(1,3,7:8,11,13)]
all <-threestage_results_late$maf_freq_in_Late >=1
threestage_results_late <- threestage_results_late[all,]
threestage_results_late <- threestage_results_late[order(-si_2),]
colnames(threestage_results_late)[2] <- "si"
threestage_results_late$progression <- rep("Late", length(threestage_results_late$variant_name))

threestage_results_met <- threestage_results[,c(1,4,9:10,11,14)]
all <-threestage_results_met$maf_freq_in_Metastasis >=1
threestage_results_met <- threestage_results_met[all,]
threestage_results_met <- threestage_results_met[order(-si_3),]
colnames(threestage_results_met)[2] <- "si"
threestage_results_met$progression <- rep("Metastasis", length(threestage_results_met$variant_name))

#Extracting recurrent variants

recurrent <- threestage_results_early$maf_freq_in_Early > 1
threestage_results_early_recur <- threestage_results_early[recurrent,]

recurrent <- threestage_results_late$maf_freq_in_Late > 1
threestage_results_late_recur <- threestage_results_late[recurrent,]

recurrent <- threestage_results_met$maf_freq_in_Metastasis > 1
threestage_results_met_recur <- threestage_results_met[recurrent,]


##########################################
#Summary function

summary_gene <- function(data) {
  data_clean <- data %>% 
    arrange(desc(si)) %>%
    filter(si > 1)
  
  # Summarise information of gene with multiple variants
  info1 <- data_clean %>% group_by(gene) %>%
    summarise(cum_si = sum(si), 
              mean_si = mean(si),
              median_si = median(si),
              sd = sd(si),
              max_si = max(si),
              n_variant = n_distinct(variant_name)) %>%
    filter(n_variant > 1)
  
  top_variant <- data_clean %>%
    group_by(gene) %>% filter(row_number() == 1)
  
  merge_info <- merge(info1, top_variant[, -3], by.x = "gene") %>%
    arrange(desc(cum_si), desc(n_variant))
  return(merge_info)
}

#Used to find genes that have at least one recurrent variants
summary_gene_recur <- function(data) {
  data_clean <- data %>% 
    arrange(desc(si)) %>%
    filter(si > 1)
  
  # Summarise information of gene with multiple variants
  info1 <- data_clean %>% group_by(gene) %>%
    summarise(cum_si = sum(si), # change sum to mean and sd
              mean_si = mean(si),
              median_si = median(si),
              sd = sd(si),
              max_si = max(si),
              n_variant = n_distinct(variant_name)) %>%
    filter(n_variant > 0)
  
  top_variant <- data_clean %>%
    group_by(gene) %>% filter(row_number() == 1)
  
  merge_info <- merge(info1, top_variant[, -3], by.x = "gene") %>%
    arrange(desc(cum_si), desc(n_variant))
  return(merge_info)
}

###############################################################################################
#Gene level SI
################################################################################################

early_data <- data.frame(variant_name = threestage_results_early$variant_name,
                         gene = threestage_results_early$gene,
                         si = threestage_results_early$si)
late_data <- data.frame(variant_name = threestage_results_late$variant_name,
                        gene = threestage_results_late$gene,
                        si = threestage_results_late$si)
met_data <- data.frame(variant_name = threestage_results_met$variant_name,
                       gene = threestage_results_met$gene,
                       si = threestage_results_met$si)

early_info <- summary_gene(early_data)
late_info <- summary_gene(late_data)
met_info <- summary_gene(met_data)

early_data_recur <- data.frame(variant_name = threestage_results_early_recur$variant_name,
                               gene = threestage_results_early_recur$gene,
                               si = threestage_results_early_recur$si)
late_data_recur <- data.frame(variant_name = threestage_results_late_recur$variant_name,
                              gene = threestage_results_late_recur$gene,
                              si = threestage_results_late_recur$si)
met_data_recur <- data.frame(variant_name = threestage_results_met_recur$variant_name,
                             gene = threestage_results_met_recur$gene,
                             si = threestage_results_met_recur$si)

early_info_recur <- summary_gene_recur(early_data_recur)
late_info_recur <- summary_gene_recur(late_data_recur)
met_info_recur <- summary_gene_recur(met_data_recur)

#Filtering out all genes that have NO recurrent variants, 
#aka filtering in genes that have at least ONE recurrent variant

early_info <- early_info[which(early_info$gene %in% early_info_recur$gene),]
late_info <- late_info[which(late_info$gene %in% late_info_recur$gene),]
met_info <- met_info[which(met_info$gene %in% met_info_recur$gene),]


fill_na <- function(x, fill = 0) {
  x = ifelse(is.na(x), fill, x)
  return(x)
}


prim_info <- merge(early_info, late_info, by = "gene", all = T, 
                   suffixes = c(".e", ".l")) %>%
  mutate_at(c("cum_si.e", "cum_si.l", 
              "mean_si.e", "mean_si.l", 
              "sd.e", "sd.l",
              "n_variant.e", "n_variant.l"), fill_na) %>%
  mutate(n_variant_prim = n_variant.e + n_variant.l, 
         mean_si_prim = (cum_si.e + cum_si.l) / n_variant_prim) %>%
  arrange(desc(n_variant_prim))

colnames(met_info) <- paste(colnames(met_info), ".m", sep = "")
colnames(met_info)[1] <- "gene"

stage_merge <- merge(prim_info, met_info, by = "gene", all = T) %>%
  mutate_at(c("cum_si.e", "cum_si.l", "cum_si.m", 
              "mean_si.e", "mean_si.l", "mean_si.m", 
              "sd.e", "sd.l", "sd.m",
              "n_variant.e", "n_variant.l", "n_variant.m"), fill_na) %>%
  mutate(n_variant_total = n_variant.e + n_variant.l + n_variant.m, 
         mean_si_total = (cum_si.e + cum_si.l + cum_si.m) / n_variant_total) %>%
  arrange(desc(n_variant_total))


########################################################################################
# Early

stage_merge_early_ordered <- stage_merge[order(-stage_merge$mean_si.e),]
selected_early_genes <- stage_merge_early_ordered$gene[1:10]

#select all variants within gene list
early_list.e <- threestage_results_early %>%
  filter(gene %in% selected_early_genes) %>%
  select(1,2,7,6,5)
early_list.l <- threestage_results_late %>%
  filter(gene %in% selected_early_genes) %>%
  select(1,2,7,6,5)
early_list.m <- threestage_results_met %>%
  filter(gene %in% selected_early_genes) %>%
  select(1,2,7,6,5)

colnames(early_list.e)[4] <- "maf_freq"
colnames(early_list.l)[4] <- "maf_freq"
colnames(early_list.m)[4] <- "maf_freq"

early_list <- rbind(early_list.e, early_list.l)
early_list <- rbind(early_list, early_list.m)

#set order of genes for plot
early_list$gene <- early_list$gene %>%
  factor(levels = selected_early_genes)

early_jitter <- ggplot(early_list, aes(x=gene, y=si, color=progression)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1))+ 
  xlab("Gene or protein") + ylab("Scaled selection coefficient") + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(color = "black"),
        panel.border = element_blank(),
        legend.position = c(0.95, 0.85))+
  scale_color_discrete(name = "Stage", labels = c("Lower-risk", "Higher-risk", "Metastasis")) +
  geom_vline(xintercept=seq(1.5, length(unique(early_list$gene))-0.5, 1), 
             lwd=.5, colour="lightgrey") +
  scale_y_continuous(labels=scientific, limits = c(-0.8e4, 1e5))

early_jitter

########################################################################################
# Late

stage_merge_late_ordered <- stage_merge[order(-stage_merge$mean_si.l),]
selected_late_genes <- stage_merge_late_ordered$gene[1:10]
  
#select all variants within gene list
late_list.e <- threestage_results_early %>%
  filter(gene %in% selected_late_genes) %>%
  select(1,2,7,6,5)
late_list.l <- threestage_results_late %>%
  filter(gene %in% selected_late_genes) %>%
  select(1,2,7,6,5)
late_list.m <- threestage_results_met %>%
  filter(gene %in% selected_late_genes) %>%
  select(1,2,7,6,5)

colnames(late_list.e)[4] <- "maf_freq"
colnames(late_list.l)[4] <- "maf_freq"
colnames(late_list.m)[4] <- "maf_freq"

late_list <- rbind(late_list.e, late_list.l)
late_list <- rbind(late_list, late_list.m)

#set order of genes for plot
late_list$gene <- late_list$gene %>%
  factor(levels = selected_late_genes)

#Dummy points

dummy_OR4N4.e <- list("OR4N4 Variant.e", as.double(0.001), "Early", "1", "OR4N4")
dummy_CHRNA6.e <- list("CHRNA6 Variant.e", as.double(0.001), "Early", "1", "CHRNA6")

late_list <- late_list %>%
  rbind(dummy_OR4N4.e) %>%
  rbind(dummy_CHRNA6.e) 

library(ggplot2)

late_jitter<- ggplot(late_list, aes(x=gene, y=si, color=progression)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1))+ 
  xlab("Gene or protein") + ylab("Scaled selection coefficient") + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(color = "black"),
        panel.border = element_blank(),
        legend.position = c(0.95, 0.85))+
  scale_color_discrete(name = "Stage", labels = c("Lower-risk", "Higher-risk", "Metastasis")) +
  geom_vline(xintercept=seq(1.5, length(unique(late_list$gene))-0.5, 1), 
             lwd=.5, colour="lightgrey") +
  scale_y_continuous(labels=scientific, limits = c(-.24e4, 3e4), breaks = c(0, 1e4, 1.5e4, 2e4, 3e4))

late_jitter

########################################################################################
# Metastasis

stage_merge_met_ordered <- stage_merge[order(-stage_merge$mean_si.m),]
selected_met_genes <- stage_merge_met_ordered$gene[1:10]

#select all variants within gene list
met_list.e <- threestage_results_early %>%
  filter(gene %in% selected_met_genes) %>%
  select(1,2,7,6,5)
met_list.l <- threestage_results_late %>%
  filter(gene %in% selected_met_genes) %>%
  select(1,2,7,6,5)
met_list.m <- threestage_results_met %>%
  filter(gene %in% selected_met_genes) %>%
  select(1,2,7,6,5)

colnames(met_list.e)[4] <- "maf_freq"
colnames(met_list.l)[4] <- "maf_freq"
colnames(met_list.m)[4] <- "maf_freq"

met_list <- rbind(met_list.e, met_list.l)
met_list <- rbind(met_list, met_list.m)

#set order of genes for plot
met_list$gene <- met_list$gene %>%
  factor(levels = selected_met_genes)

#Dummy points

dummy_ORC3.e <- list("ORC3 Variant.e", as.double(0.001), "Early", "1", "ORC3")
dummy_ORC3.l <- list("ORC3 Variant.l", as.double(0.001), "Late", "1", "ORC3")
dummy_ZNF780B.l <- list("ZNF780B Variant.l", as.double(0.001), "Late", "1", "ZNF780B")
dummy_DIMT1.e <- list("DIMT1 Variant.e", as.double(0.001), "Early", "1", "DIMT1")
dummy_DIMT1.l <- list("DIMT1 Variant.l", as.double(0.001), "Late", "1", "DIMT1")
dummy_KRTAP13_3.l <- list("KRTAP13-3 Variant.l", as.double(0.001), "Late", "1", "KRTAP13-3")
dummy_ZNF714.e <- list("ZNF714 Variant.e", as.double(0.001), "Early", "1", "ZNF714")
dummy_GRB7.e <- list("GRB7 Variant.e", as.double(0.001), "Early", "1", "GRB7")
dummy_APCS.e <- list("APCS Variant.e", as.double(0.001), "Early", "1", "APCS")

met_list <- met_list %>%
  rbind(dummy_ORC3.e) %>%
  rbind(dummy_ORC3.l) %>%
  rbind(dummy_ZNF780B.l)%>%
  rbind(dummy_DIMT1.e)%>%
  rbind(dummy_DIMT1.l)%>%
  rbind(dummy_KRTAP13_3.l)%>%
  rbind(dummy_ZNF714.e)%>%
  rbind(dummy_GRB7.e)%>%
  rbind(dummy_APCS.e)

library(ggplot2)

met_jitter <- ggplot(met_list, aes(x=gene, y=si, color=progression)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1))+
  xlab("Gene or protein") + ylab("Scaled selection coefficient") + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(color = "black"),
        panel.border = element_blank(),
        legend.position = c(0.95, 0.85))+
  scale_color_discrete(name = "Stage", labels = c("Lower-risk", "Higher-risk", "Metastasis")) +
  geom_vline(xintercept=seq(1.5, length(unique(met_list$gene))-0.5, 1), 
             lwd=.5, colour="lightgrey") +
  scale_y_continuous(labels=scientific, limits = c(-0.9e4, 1.051e5))

met_jitter
##############################################################################
#Wilcoxon test

early.e <- list()
early.l <- list()
early.m <- list()

for (x in unique(early_list$gene)) {
  values <- early_list %>%
    filter(gene == x) %>%
    filter(progression == "Early") %>%
    pull(si)
  early.e <- c(early.e, list(values))
}
names(early.e) <- unique(early_list$gene)

for (x in unique(early_list$gene)) {
  values <- early_list %>%
    filter(gene == x) %>%
    filter(progression == "Late") %>%
    pull(si)
  early.l <- c(early.l, list(values))
}
names(early.l) <- unique(early_list$gene)

for (x in unique(early_list$gene)) {
  values <- early_list %>%
    filter(gene == x) %>%
    filter(progression == "Metastasis") %>%
    pull(si)
  early.m <- c(early.m, list(values))
}
names(early.m) <- unique(early_list$gene)

late.e <- list()
late.l <- list()
late.m <- list()

for (x in unique(late_list$gene)) {
  values <- late_list %>%
    filter(gene == x) %>%
    filter(progression == "Early") %>%
    pull(si)
  late.e <- c(late.e, list(values))
}
names(late.e) <- unique(late_list$gene)

for (x in unique(late_list$gene)) {
  values <- late_list %>%
    filter(gene == x) %>%
    filter(progression == "Late") %>%
    pull(si)
  late.l <- c(late.l, list(values))
}
names(late.l) <- unique(late_list$gene)

for (x in unique(late_list$gene)) {
  values <- late_list %>%
    filter(gene == x) %>%
    filter(progression == "Metastasis") %>%
    pull(si)
  late.m <- c(late.m, list(values))
}
names(late.m) <- unique(late_list$gene)

met.e <- list()
met.l <- list()
met.m <- list()

for (x in unique(met_list$gene)) {
  values <- met_list %>%
    filter(gene == x) %>%
    filter(progression == "Early") %>%
    pull(si)
  met.e <- c(met.e, list(values))
}
names(met.e) <- unique(met_list$gene)

for (x in unique(met_list$gene)) {
  values <- met_list %>%
    filter(gene == x) %>%
    filter(progression == "Late") %>%
    pull(si)
  met.l <- c(met.l, list(values))
}
names(met.l) <- unique(met_list$gene)

for (x in unique(met_list$gene)) {
  values <- met_list %>%
    filter(gene == x) %>%
    filter(progression == "Metastasis") %>%
    pull(si)
  met.m <- c(met.m, list(values))
}
names(met.m) <- unique(met_list$gene)


#Tests for Lower-risk
wilcox_early.e_l <- c()
for (x in 1:10){
  wilcox <- wilcox.test(early.e[[x]], early.l[[x]])
  wilcox_early.e_l <- c(wilcox_early.e_l, wilcox$p.value)
}

wilcox_early.e_m <- c()
for (x in 1:10){
  wilcox <- wilcox.test(early.e[[x]], early.m[[x]])
  wilcox_early.e_m <- c(wilcox_early.e_m, wilcox$p.value)
}

wilcox_early.l_m <- c()
for (x in 1:10){
  wilcox <- wilcox.test(early.l[[x]], early.m[[x]])
  wilcox_early.l_m <- c(wilcox_early.l_m, wilcox$p.value)
}

early.wilcox <- data.frame(early_late = wilcox_early.e_l,
                           early_met = wilcox_early.e_m,
                           late_met = wilcox_early.l_m)
row.names(early.wilcox) <- unique(early_list$gene)

#Tests for Higher-risk
wilcox_late.e_l <- c()
for (x in 1:10){
  wilcox <- wilcox.test(late.e[[x]], late.l[[x]])
  wilcox_late.e_l <- c(wilcox_late.e_l, wilcox$p.value)
}

wilcox_late.e_m <- c()
for (x in 1:10){
  wilcox <- wilcox.test(late.e[[x]], late.m[[x]])
  wilcox_late.e_m <- c(wilcox_late.e_m, wilcox$p.value)
}

wilcox_late.l_m <- c()
for (x in 1:10){
  wilcox <- wilcox.test(late.l[[x]], late.m[[x]])
  wilcox_late.l_m <- c(wilcox_late.l_m, wilcox$p.value)
}

late.wilcox <- data.frame(early_late = wilcox_late.e_l,
                           early_met = wilcox_late.e_m,
                           late_met = wilcox_late.l_m)
row.names(late.wilcox) <- unique(late_list$gene)

#Tests for Metastasis
wilcox_met.e_l <- c()
for (x in 1:10){
  wilcox <- wilcox.test(met.e[[x]], met.l[[x]])
  wilcox_met.e_l <- c(wilcox_met.e_l, wilcox$p.value)
}

wilcox_met.e_m <- c()
for (x in 1:10){
  wilcox <- wilcox.test(met.e[[x]], met.m[[x]])
  wilcox_met.e_m <- c(wilcox_met.e_m, wilcox$p.value)
}

wilcox_met.l_m <- c()
for (x in 1:10){
  wilcox <- wilcox.test(met.l[[x]], met.m[[x]])
  wilcox_met.l_m <- c(wilcox_met.l_m, wilcox$p.value)
}

met.wilcox <- data.frame(early_late = wilcox_met.e_l,
                          early_met = wilcox_met.e_m,
                          late_met = wilcox_met.l_m)
row.names(met.wilcox) <- unique(met_list$gene)

# write.table(early.wilcox, file = "wilcox_early.txt", sep = "\t",
#             row.names = TRUE, col.names = TRUE, quote = FALSE)
# write.table(late.wilcox, file = "wilcox_late.txt", sep = "\t",
#             row.names = TRUE, col.names = TRUE, quote = FALSE)
# write.table(met.wilcox, file = "wilcox_met.txt", sep = "\t",
#             row.names = TRUE, col.names = TRUE, quote = FALSE)

###########################################################################
#Extracting data for Bayesian approach
#Extracting SI + 95% CI for each variant in each of the featured top 10 genes

#Early
early_list.e <- threestage_results_early %>%
  filter(gene %in% selected_early_genes) %>%
  select(1,2,3,4,7,6,5)
early_list.l <- threestage_results_late %>%
  filter(gene %in% selected_early_genes) %>%
  select(1,2,3,4,7,6,5)
early_list.m <- threestage_results_met %>%
  filter(gene %in% selected_early_genes) %>%
  select(1,2,3,4,7,6,5)

colnames(early_list.e)[3] <- "ci_low_95"
colnames(early_list.l)[3] <- "ci_low_95"
colnames(early_list.m)[3] <- "ci_low_95"

colnames(early_list.e)[4] <- "ci_high_95"
colnames(early_list.l)[4] <- "ci_high_95"
colnames(early_list.m)[4] <- "ci_high_95"

colnames(early_list.e)[6] <- "maf_freq"
colnames(early_list.l)[6] <- "maf_freq"
colnames(early_list.m)[6] <- "maf_freq"

early_list <- rbind(early_list.e, early_list.l)
early_list <- rbind(early_list, early_list.m)

#Late
late_list.e <- threestage_results_early %>%
  filter(gene %in% selected_late_genes) %>%
  select(1,2,3,4,7,6,5)
late_list.l <- threestage_results_late %>%
  filter(gene %in% selected_late_genes) %>%
  select(1,2,3,4,7,6,5)
late_list.m <- threestage_results_met %>%
  filter(gene %in% selected_late_genes) %>%
  select(1,2,3,4,7,6,5)

colnames(late_list.e)[3] <- "ci_low_95"
colnames(late_list.l)[3] <- "ci_low_95"
colnames(late_list.m)[3] <- "ci_low_95"

colnames(late_list.e)[4] <- "ci_high_95"
colnames(late_list.l)[4] <- "ci_high_95"
colnames(late_list.m)[4] <- "ci_high_95"

colnames(late_list.e)[6] <- "maf_freq"
colnames(late_list.l)[6] <- "maf_freq"
colnames(late_list.m)[6] <- "maf_freq"

late_list <- rbind(late_list.e, late_list.l)
late_list <- rbind(late_list, late_list.m)

#Metastasis
met_list.e <- threestage_results_early %>%
  filter(gene %in% selected_met_genes) %>%
  select(1,2,3,4,7,6,5)
met_list.l <- threestage_results_late %>%
  filter(gene %in% selected_met_genes) %>%
  select(1,2,3,4,7,6,5)
met_list.m <- threestage_results_met %>%
  filter(gene %in% selected_met_genes) %>%
  select(1,2,3,4,7,6,5)

colnames(met_list.e)[3] <- "ci_low_95"
colnames(met_list.l)[3] <- "ci_low_95"
colnames(met_list.m)[3] <- "ci_low_95"

colnames(met_list.e)[4] <- "ci_high_95"
colnames(met_list.l)[4] <- "ci_high_95"
colnames(met_list.m)[4] <- "ci_high_95"

colnames(met_list.e)[6] <- "maf_freq"
colnames(met_list.l)[6] <- "maf_freq"
colnames(met_list.m)[6] <- "maf_freq"

met_list <- rbind(met_list.e, met_list.l)
met_list <- rbind(met_list, met_list.m)

############################################################################

write.table(early_list, file = "lower-risk_gene_variants.txt", sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(late_list, file = "higher-risk_gene_variants.txt", sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(met_list, file = "metastatic_gene_variants.txt", sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)

##############################################################################
# Variants

early_top10 <- threestage_results_early_recur[1:10,]
late_top10 <- threestage_results_late_recur[1:10,]
met_top10 <- threestage_results_met_recur[1:10,]

unique_jitter_early_all <- ggplot(data = early_top10, aes(x = 1, y = si, color=progression)) +
  geom_jitter(position=position_jitter(0.0, seed = 5))+
  theme(axis.ticks.x = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(color = "black"),
        legend.position="none")+
  xlab("Top 10 variants") + ylab("") +
  scale_color_manual(values = "#F8766D") + 
  scale_x_continuous(limits = c(0.975, 1.25)) +
  scale_y_continuous(labels=scientific, limits = c(-0.8e4, 1e5))

unique_jitter_late_all<- ggplot(data = late_top10, aes(x = 1, y = si, color=progression)) +
  geom_jitter(position=position_jitter(0.0, seed = 5))+
  theme(axis.ticks.x = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(color = "black"),
        legend.position="none")+
  xlab("Top 10 variants") + ylab("") +
  scale_color_manual(values = "#00BA38") + 
  scale_x_continuous(limits = c(0.975, 1.25)) +
  scale_y_continuous(labels=scientific, limits = c(-.24e4, 3e4), breaks = c(0, 1e4, 1.5e4, 2e4, 3.0e4))

unique_jitter_met_all<- ggplot(data = met_top10, aes(x = 1, y = si, color=progression)) +
  geom_jitter(position=position_jitter(0.0, seed = 5))+
  theme(axis.ticks.x = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(color = "black"),
        legend.position="none")+
  xlab("Top 10 variants") + ylab("") +
  scale_color_manual(values = "#619CFF") + 
  scale_x_continuous(limits = c(0.975, 1.25)) +
  scale_y_continuous(labels=scientific, limits = c(-0.9e4, 1.051e5))

#########################################
early_title <- ggdraw() + 
  draw_label(
    "Lower-risk",
    fontface = 'bold',
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 49)
  )

late_title <- ggdraw() + 
  draw_label(
    "Higher-risk",
    fontface = 'bold',
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 49)
  )

met_title <- ggdraw() + 
  draw_label(
    "Metastases",
    fontface = 'bold',
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 49)
  )


#######################################
#Combining all the graphs + titles

early_combined <- plot_grid(early_jitter, unique_jitter_early_all, 
                            align = "h", axis = "t", nrow = 1, ncol = 2, scale = 1, rel_widths = c(5,1),
                            labels = c("A", "B"), label_size = 10)
early_combined_title <- plot_grid(early_title, early_combined, 
                                  align = "h", axis = "t", nrow = 2, ncol = 1, scale = 1, rel_heights = c(0.1,1))
late_combined <- plot_grid(late_jitter, unique_jitter_late_all, 
                           align = "h", axis = "t", nrow = 1, ncol = 2, scale = 1, rel_widths = c(5,1),
                           labels = c("C", "D"), label_size = 10)
late_combined_title <- plot_grid(late_title, late_combined, 
                                 align = "h", axis = "t", nrow = 2, ncol = 1, scale = 1, rel_heights = c(0.1,1))
met_combined <- plot_grid(met_jitter, unique_jitter_met_all, 
                          align = "h", axis = "t", nrow = 1, ncol = 2, scale = 1, rel_widths = c(5,1),
                          labels = c("E", "F"), label_size = 10)
met_combined_title <- plot_grid(met_title, met_combined, 
                                align = "h", axis = "t", nrow = 2, ncol = 1, scale = 1, rel_heights = c(0.1,1))
bar_box_combined <- plot_grid(early_combined_title, late_combined_title, met_combined_title,
                              align = "h", axis = "t", nrow = 3, ncol = 1, scale = 1)

bar_box_combined
ggsave("PRAD_figures/bar_jitter_test.png", width = 12.5, height = 12.5)

