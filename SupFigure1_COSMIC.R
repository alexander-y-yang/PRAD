library(cancereffectsizeR
library(ggplot2)
library(cowplot)

twostage_final <- load_cesa("twostage_final.rds")

PRAD_analysis <- twostage_final
PRAD_results <- snv_results(PRAD_analysis)
PRAD_results <- PRAD_results$selection.1

#seperating into prim and met
tumor_names <- unique(PRAD_analysis@maf$Unique_Patient_Identifier)

tumor_prim <- c()
tumor_met <- c()

samples <- PRAD_analysis@samples

for (x in tumor_names){
  if (samples$group[which(samples$Unique_Patient_Identifier == x)] == "Primary"){
    tumor_prim <- c(tumor_prim, x)
  }
  else if (samples$group[which(samples$Unique_Patient_Identifier == x)] == "Metastasis"){
    tumor_met <- c(tumor_met, x)
  }
}

signature_table <- PRAD_analysis@trinucleotide_mutation_weights$signature_weight_table

tumor_prim_signature <- data.frame(matrix(ncol = 72, nrow = 0))
colnames(tumor_prim_signature) <- colnames(signature_table)[5:76]

tumor_met_signature <- data.frame(matrix(ncol = 72, nrow = 0))
colnames(tumor_met_signature) <- colnames(signature_table)[5:76]

#Separating primary and tumors and extracting their cosmic signatures
for(x in tumor_prim){
  tumor_prim_signature <- rbind(tumor_prim_signature, signature_table[which(signature_table$Unique_Patient_Identifier == x),5:76])
}


for(x in tumor_met){
  tumor_met_signature <- rbind(tumor_met_signature, signature_table[which(signature_table$Unique_Patient_Identifier == x),5:76])
}

#Remove columns that do not have any weight
#Remove SBS10a, SBS10b, SBS25, SBS27, SBS31, SBS32, SBS35, SBS43, SBS45-60
remove_col <- which(colnames(tumor_prim_signature) == "SBS4"|
                      colnames(tumor_prim_signature) == "SBS7a"|
                      colnames(tumor_prim_signature) == "SBS7b"|
                      colnames(tumor_prim_signature) == "SBS7c"|
                      colnames(tumor_prim_signature) == "SBS7d"|
                      colnames(tumor_prim_signature) == "SBS9"|
                      colnames(tumor_prim_signature) == "SBS10a"|
                      colnames(tumor_prim_signature) == "SBS10b"|
                      colnames(tumor_prim_signature) == "SBS11"|
                      colnames(tumor_prim_signature) == "SBS14"|
                      colnames(tumor_prim_signature) == "SBS15"|
                      colnames(tumor_prim_signature) == "SBS16"|
                      colnames(tumor_prim_signature) == "SBS17a"|
                      colnames(tumor_prim_signature) == "SBS17b"|
                      colnames(tumor_prim_signature) == "SBS19"|
                      colnames(tumor_prim_signature) == "SBS20"|
                      colnames(tumor_prim_signature) == "SBS21"|
                      colnames(tumor_prim_signature) == "SBS22"|
                      colnames(tumor_prim_signature) == "SBS23"|
                      colnames(tumor_prim_signature) == "SBS24"|
                      colnames(tumor_prim_signature) == "SBS25"|
                      colnames(tumor_prim_signature) == "SBS26"|
                      colnames(tumor_prim_signature) == "SBS27"|
                      colnames(tumor_prim_signature) == "SBS28"|
                      colnames(tumor_prim_signature) == "SBS29"|
                      colnames(tumor_prim_signature) == "SBS30"|
                      colnames(tumor_prim_signature) == "SBS31"|
                      colnames(tumor_prim_signature) == "SBS32"|
                      colnames(tumor_prim_signature) == "SBS34"|
                      colnames(tumor_prim_signature) == "SBS35"|
                      colnames(tumor_prim_signature) == "SBS36"|
                      colnames(tumor_prim_signature) == "SBS38"|
                      colnames(tumor_prim_signature) == "SBS42"|
                      colnames(tumor_prim_signature) == "SBS43"|
                      colnames(tumor_prim_signature) == "SBS44"|
                      colnames(tumor_prim_signature) == "SBS45"|
                      colnames(tumor_prim_signature) == "SBS46"|
                      colnames(tumor_prim_signature) == "SBS47"|
                      colnames(tumor_prim_signature) == "SBS48"|
                      colnames(tumor_prim_signature) == "SBS49"|
                      colnames(tumor_prim_signature) == "SBS50"|
                      colnames(tumor_prim_signature) == "SBS51"|
                      colnames(tumor_prim_signature) == "SBS52"|
                      colnames(tumor_prim_signature) == "SBS53"|
                      colnames(tumor_prim_signature) == "SBS54"|
                      colnames(tumor_prim_signature) == "SBS55"|
                      colnames(tumor_prim_signature) == "SBS56"|
                      colnames(tumor_prim_signature) == "SBS57"|
                      colnames(tumor_prim_signature) == "SBS58"|
                      colnames(tumor_prim_signature) == "SBS59"|
                      colnames(tumor_prim_signature) == "SBS60"|
                      colnames(tumor_prim_signature) == "SBS84"|
                      colnames(tumor_prim_signature) == "SBS85"|
                      colnames(tumor_prim_signature) == "SBS88")


tumor_prim_signature <- subset(tumor_prim_signature,
                               select = -remove_col)

tumor_met_signature <- subset(tumor_met_signature,
                              select = -remove_col)

library(stringr)
colnames(tumor_prim_signature) <- str_replace(colnames(tumor_prim_signature), "SBS", "")
colnames(tumor_met_signature) <- str_replace(colnames(tumor_met_signature), "SBS", "")

library(reshape2)
#Need to convert data frame from wide to long format
tumor_prim_signature <- as.data.frame(t(as.matrix(tumor_prim_signature)))
tumor_prim_signature$group <- row.names(tumor_prim_signature)
tumor_prim_signature.m <- melt(tumor_prim_signature, id.vars = "group")

tumor_met_signature <- as.data.frame(t(as.matrix(tumor_met_signature)))
tumor_met_signature$group <- row.names(tumor_met_signature)
tumor_met_signature.m <- melt(tumor_met_signature, id.vars = "group")

#Ordering the signatures for ggplot2
tumor_prim_signature.m$group <- factor(tumor_prim_signature.m$group, levels=unique(as.character(tumor_prim_signature.m$group)))
tumor_met_signature.m$group <- factor(tumor_met_signature.m$group, levels=unique(as.character(tumor_met_signature.m$group)))

save(tumor_prim_signature.m, file="PRAD_tumor_prim_signature.RData")
save(tumor_met_signature.m, file="PRAD_tumor_met_signature.RData")

tumor_prim_signature_plot <- ggplot(data=tumor_prim_signature.m, aes(x=group, y=value, fill=group)) +
  #  geom_dotplot(binaxis="y", binwidth=.001, stackdir="center") +
  geom_violin(scale="width") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text (hjust = 1, angle = 45)) +
  theme(panel.background = element_blank()) +
  theme(panel.grid.major.y = element_line(color="lightgrey"), panel.grid.minor.y = element_line(color="lightgrey")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  xlab("COSMIC Signature") +
  ylab("Signature Weight") + ylim(0, 0.7) +
  ggtitle("Primary")

tumor_met_signature_plot <- ggplot(data=tumor_met_signature.m, aes(x=group, y=value, fill=group)) +
  geom_violin(scale="width") +
  theme(plot.title = element_text(hjust = 0), axis.text.x = element_text (hjust = 0.5)) +
  theme(panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  xlab("COSMIC Signature") +
  ylab("Signature Weight") + ylim(0, 1) +
  ggtitle("Metastasis")

tumor_overlay_signature_plot <- ggplot() +
  geom_violin(data=tumor_met_signature.m, aes(x=group, y=value, fill=group), scale="width", alpha=0.3)+
  geom_violin(data=tumor_prim_signature.m, aes(x=group, y=value, fill=group), scale="width", width = 0.4) +
  theme(plot.title = element_text(hjust = 0), axis.text.x = element_text (hjust = 0.5)) +
  theme(panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  xlab("COSMIC signature") +
  ylab("Signature weight") + ylim(0, 1) +
  ggtitle("Primary vs. metastasis")


tumor_met_signature_plot
ggsave("PRAD_figures/PRAD_cosmic_signatures_met.png", width = 12, height = 8)

tumor_overlay_signature_plot
ggsave("PRAD_figures/PRAD_cosmic_signatures_overlay.png", width = 12, height = 5)
#####################################################################################################################
threestage_final <- load_cesa("threestage_final.rds")
threestage_results <- snv_results(threestage_final)
threestage_results <- threestage_results$selection.1

#seperating into prim and met
tumor_names_threestage <- unique(threestage_final@maf$Unique_Patient_Identifier)

tumor_early <- c()
tumor_late <- c()

samples_threestage <- threestage_final@samples
signature_table_threestage <- threestage_final@trinucleotide_mutation_weights$signature_weight_table

for (x in tumor_names_threestage){
  if (samples_threestage$group[which(samples_threestage$Unique_Patient_Identifier == x)] == "Early"){
    tumor_early <- c(tumor_early, x)
  }
  else if (samples_threestage$group[which(samples_threestage$Unique_Patient_Identifier == x)] == "Late"){
    tumor_late <- c(tumor_late, x)
  }
}

tumor_early_signature <- data.frame(matrix(ncol = 72, nrow = 0))
colnames(tumor_early_signature) <- colnames(signature_table_threestage)[5:76]

tumor_late_signature <- data.frame(matrix(ncol = 72, nrow = 0))
colnames(tumor_late_signature) <- colnames(signature_table_threestage)[5:76]

#Separating primary and tumors and extracting their cosmic signatures
for(x in tumor_early){
  tumor_early_signature <- rbind(tumor_early_signature, signature_table_threestage[which(signature_table_threestage$Unique_Patient_Identifier == x),5:76])
}


for(x in tumor_late){
  tumor_late_signature <- rbind(tumor_late_signature, signature_table_threestage[which(signature_table_threestage$Unique_Patient_Identifier == x),5:76])
}

#Remove columns that do not have any weight
#Remove SBS10a, SBS10b, SBS25, SBS27, SBS31, SBS32, SBS35, SBS43, SBS45-60
remove_col <- which(colnames(tumor_early_signature) == "SBS4"|
                      colnames(tumor_early_signature) == "SBS7a"|
                      colnames(tumor_early_signature) == "SBS7b"|
                      colnames(tumor_early_signature) == "SBS7c"|
                      colnames(tumor_early_signature) == "SBS7d"|
                      colnames(tumor_early_signature) == "SBS9"|
                      colnames(tumor_early_signature) == "SBS10a"|
                      colnames(tumor_early_signature) == "SBS10b"|
                      colnames(tumor_early_signature) == "SBS11"|
                      colnames(tumor_early_signature) == "SBS14"|
                      colnames(tumor_early_signature) == "SBS15"|
                      colnames(tumor_early_signature) == "SBS16"|
                      colnames(tumor_early_signature) == "SBS17a"|
                      colnames(tumor_early_signature) == "SBS17b"|
                      colnames(tumor_early_signature) == "SBS19"|
                      colnames(tumor_early_signature) == "SBS20"|
                      colnames(tumor_early_signature) == "SBS21"|
                      colnames(tumor_early_signature) == "SBS22"|
                      colnames(tumor_early_signature) == "SBS23"|
                      colnames(tumor_early_signature) == "SBS24"|
                      colnames(tumor_early_signature) == "SBS25"|
                      colnames(tumor_early_signature) == "SBS26"|
                      colnames(tumor_early_signature) == "SBS27"|
                      colnames(tumor_early_signature) == "SBS28"|
                      colnames(tumor_early_signature) == "SBS29"|
                      colnames(tumor_early_signature) == "SBS30"|
                      colnames(tumor_early_signature) == "SBS31"|
                      colnames(tumor_early_signature) == "SBS32"|
                      colnames(tumor_early_signature) == "SBS34"|
                      colnames(tumor_early_signature) == "SBS35"|
                      colnames(tumor_early_signature) == "SBS36"|
                      colnames(tumor_early_signature) == "SBS38"|
                      colnames(tumor_early_signature) == "SBS42"|
                      colnames(tumor_early_signature) == "SBS43"|
                      colnames(tumor_early_signature) == "SBS44"|
                      colnames(tumor_early_signature) == "SBS45"|
                      colnames(tumor_early_signature) == "SBS46"|
                      colnames(tumor_early_signature) == "SBS47"|
                      colnames(tumor_early_signature) == "SBS48"|
                      colnames(tumor_early_signature) == "SBS49"|
                      colnames(tumor_early_signature) == "SBS50"|
                      colnames(tumor_early_signature) == "SBS51"|
                      colnames(tumor_early_signature) == "SBS52"|
                      colnames(tumor_early_signature) == "SBS53"|
                      colnames(tumor_early_signature) == "SBS54"|
                      colnames(tumor_early_signature) == "SBS55"|
                      colnames(tumor_early_signature) == "SBS56"|
                      colnames(tumor_early_signature) == "SBS57"|
                      colnames(tumor_early_signature) == "SBS58"|
                      colnames(tumor_early_signature) == "SBS59"|
                      colnames(tumor_early_signature) == "SBS60"|
                      colnames(tumor_early_signature) == "SBS84"|
                      colnames(tumor_early_signature) == "SBS85"|
                      colnames(tumor_early_signature) == "SBS88")


tumor_early_signature <- subset(tumor_early_signature,
                                select = -remove_col)

tumor_late_signature <- subset(tumor_late_signature,
                               select = -remove_col)

library(stringr)
colnames(tumor_early_signature) <- str_replace(colnames(tumor_early_signature), "SBS", "")
colnames(tumor_late_signature) <- str_replace(colnames(tumor_late_signature), "SBS", "")

library(reshape2)
#Need to convert data frame from wide to long format
tumor_early_signature <- as.data.frame(t(as.matrix(tumor_early_signature)))
tumor_early_signature$group <- row.names(tumor_early_signature)
tumor_early_signature.m <- melt(tumor_early_signature, id.vars = "group")

tumor_late_signature <- as.data.frame(t(as.matrix(tumor_late_signature)))
tumor_late_signature$group <- row.names(tumor_late_signature)
tumor_late_signature.m <- melt(tumor_late_signature, id.vars = "group")

#Ordering the signatures for ggplot2
tumor_early_signature.m$group <- factor(tumor_early_signature.m$group, levels=unique(as.character(tumor_early_signature.m$group)))
tumor_late_signature.m$group <- factor(tumor_late_signature.m$group, levels=unique(as.character(tumor_late_signature.m$group)))

save(tumor_early_signature.m, file="PRAD_tumor_early_signature.RData")
save(tumor_late_signature.m, file="PRAD_tumor_late_signature.RData")

tumor_early_signature_plot <- ggplot(data=tumor_early_signature.m, aes(x=group, y=value, fill=group)) +
  geom_violin(scale="width") +
  theme(plot.title = element_text(hjust = 0), axis.text.x = element_text (hjust = 0.5)) +
  theme(panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  xlab("COSMIC Signature") +
  ylab("Signature Weight") + ylim(0, 1) +
  ggtitle("Lower-risk")

tumor_late_signature_plot <- ggplot(data=tumor_late_signature.m, aes(x=group, y=value, fill=group)) +
  geom_violin(scale="width") +
  theme(plot.title = element_text(hjust = 0), axis.text.x = element_text (hjust = 0.5)) +
  theme(panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  xlab("COSMIC Signature") +
  ylab("Signature Weight") + ylim(0, 1) +
  ggtitle("Higher-risk")

#combining all three plots

combined_PRAD_signature <- plot_grid(tumor_early_signature_plot, tumor_late_signature_plot, tumor_met_signature_plot,
                                     ncol=1, nrow=3, labels = c("A", "B", "C"), label_size = 12)
combined_PRAD_signature
ggsave("PRAD_figures/combined_PRAD_signature.png", width = 6, height = 9)
