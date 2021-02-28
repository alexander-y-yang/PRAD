library(cancereffectsizeR)
library(ggplot2)
library(cowplot)

twostage_final <- load_cesa("twostage_final.rds")

PRAD_analysis <- twostage_final
PRAD_results <- snv_results(twostage_final)
PRAD_results <- PRAD_results$selection.1

#Using Vincent's code from the HNSCC paper as a framework (not very efficient, but it works)
load("trinuc_data_HNSC_HPVpos.RData")

#Extracting trinuc data

PRAD_trinuc <- PRAD_analysis@trinucleotide_mutation_weights$trinuc_proportion_matrix

#Seperating trinuc data into prim and met
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

tumor_trinuc_prim <- data.frame(matrix(ncol = 96, nrow = 0))
tumor_trinuc_met <- data.frame(matrix(ncol = 96, nrow = 0))
colnames(tumor_trinuc_met) <- colnames(PRAD_trinuc)

for (x in tumor_prim){
  tumor_trinuc_prim <- rbind(tumor_trinuc_prim, PRAD_trinuc[x,])
}
colnames(tumor_trinuc_prim) <- colnames(PRAD_trinuc)

for(x in tumor_met){
  tumor_trinuc_met <- rbind(tumor_trinuc_met, PRAD_trinuc[x,])
}
colnames(tumor_trinuc_met) <- colnames(PRAD_trinuc)

#Averaging trinuc data
tumor_trinuc_prim_avg <- apply(tumor_trinuc_prim, 2, mean)
tumor_trinuc_met_avg <- apply(tumor_trinuc_met, 2, mean)

tumor_trinuc_prim_avg_ordered <- c()
tumor_trinuc_met_avg_ordered <- c()

#Reordering tumor_trinuc_prim_avg
#N_A
for(x in 1:96){
  if (grepl("A", substr(colnames(tumor_trinuc_prim)[x], 1, 1)) & grepl("A", substr(colnames(tumor_trinuc_prim)[x], 7, 7)) == TRUE){
    tumor_trinuc_prim_avg_ordered <- c(tumor_trinuc_prim_avg_ordered, tumor_trinuc_prim_avg[x])
  }
}
for(x in 1:96){
  if (grepl("C", substr(colnames(tumor_trinuc_prim)[x], 1, 1)) & grepl("A", substr(colnames(tumor_trinuc_prim)[x], 7, 7)) == TRUE){
    tumor_trinuc_prim_avg_ordered <- c(tumor_trinuc_prim_avg_ordered, tumor_trinuc_prim_avg[x])
  }
}
for(x in 1:96){
  if (grepl("G", substr(colnames(tumor_trinuc_prim)[x], 1, 1)) & grepl("A", substr(colnames(tumor_trinuc_prim)[x], 7, 7)) == TRUE){
    tumor_trinuc_prim_avg_ordered <- c(tumor_trinuc_prim_avg_ordered, tumor_trinuc_prim_avg[x])
  }
}
for(x in 1:96){
  if (grepl("T", substr(colnames(tumor_trinuc_prim)[x], 1, 1)) & grepl("A", substr(colnames(tumor_trinuc_prim)[x], 7, 7)) == TRUE){
    tumor_trinuc_prim_avg_ordered <- c(tumor_trinuc_prim_avg_ordered, tumor_trinuc_prim_avg[x])
  }
}
#N_C
for(x in 1:96){
  if (grepl("A", substr(colnames(tumor_trinuc_prim)[x], 1, 1)) & grepl("C", substr(colnames(tumor_trinuc_prim)[x], 7, 7)) == TRUE){
    tumor_trinuc_prim_avg_ordered <- c(tumor_trinuc_prim_avg_ordered, tumor_trinuc_prim_avg[x])
  }
}
for(x in 1:96){
  if (grepl("C", substr(colnames(tumor_trinuc_prim)[x], 1, 1)) & grepl("C", substr(colnames(tumor_trinuc_prim)[x], 7, 7)) == TRUE){
    tumor_trinuc_prim_avg_ordered <- c(tumor_trinuc_prim_avg_ordered, tumor_trinuc_prim_avg[x])
  }
}
for(x in 1:96){
  if (grepl("G", substr(colnames(tumor_trinuc_prim)[x], 1, 1)) & grepl("C", substr(colnames(tumor_trinuc_prim)[x], 7, 7)) == TRUE){
    tumor_trinuc_prim_avg_ordered <- c(tumor_trinuc_prim_avg_ordered, tumor_trinuc_prim_avg[x])
  }
}
for(x in 1:96){
  if (grepl("T", substr(colnames(tumor_trinuc_prim)[x], 1, 1)) & grepl("C", substr(colnames(tumor_trinuc_prim)[x], 7, 7)) == TRUE){
    tumor_trinuc_prim_avg_ordered <- c(tumor_trinuc_prim_avg_ordered, tumor_trinuc_prim_avg[x])
  }
}
#N_G
for(x in 1:96){
  if (grepl("A", substr(colnames(tumor_trinuc_prim)[x], 1, 1)) & grepl("G", substr(colnames(tumor_trinuc_prim)[x], 7, 7)) == TRUE){
    tumor_trinuc_prim_avg_ordered <- c(tumor_trinuc_prim_avg_ordered, tumor_trinuc_prim_avg[x])
  }
}
for(x in 1:96){
  if (grepl("C", substr(colnames(tumor_trinuc_prim)[x], 1, 1)) & grepl("G", substr(colnames(tumor_trinuc_prim)[x], 7, 7)) == TRUE){
    tumor_trinuc_prim_avg_ordered <- c(tumor_trinuc_prim_avg_ordered, tumor_trinuc_prim_avg[x])
  }
}
for(x in 1:96){
  if (grepl("G", substr(colnames(tumor_trinuc_prim)[x], 1, 1)) & grepl("G", substr(colnames(tumor_trinuc_prim)[x], 7, 7)) == TRUE){
    tumor_trinuc_prim_avg_ordered <- c(tumor_trinuc_prim_avg_ordered, tumor_trinuc_prim_avg[x])
  }
}
for(x in 1:96){
  if (grepl("T", substr(colnames(tumor_trinuc_prim)[x], 1, 1)) & grepl("G", substr(colnames(tumor_trinuc_prim)[x], 7, 7)) == TRUE){
    tumor_trinuc_prim_avg_ordered <- c(tumor_trinuc_prim_avg_ordered, tumor_trinuc_prim_avg[x])
  }
}
#N_T
for(x in 1:96){
  if (grepl("A", substr(colnames(tumor_trinuc_prim)[x], 1, 1)) & grepl("T", substr(colnames(tumor_trinuc_prim)[x], 7, 7)) == TRUE){
    tumor_trinuc_prim_avg_ordered <- c(tumor_trinuc_prim_avg_ordered, tumor_trinuc_prim_avg[x])
  }
}
for(x in 1:96){
  if (grepl("C", substr(colnames(tumor_trinuc_prim)[x], 1, 1)) & grepl("T", substr(colnames(tumor_trinuc_prim)[x], 7, 7)) == TRUE){
    tumor_trinuc_prim_avg_ordered <- c(tumor_trinuc_prim_avg_ordered, tumor_trinuc_prim_avg[x])
  }
}
for(x in 1:96){
  if (grepl("G", substr(colnames(tumor_trinuc_prim)[x], 1, 1)) & grepl("T", substr(colnames(tumor_trinuc_prim)[x], 7, 7)) == TRUE){
    tumor_trinuc_prim_avg_ordered <- c(tumor_trinuc_prim_avg_ordered, tumor_trinuc_prim_avg[x])
  }
}
for(x in 1:96){
  if (grepl("T", substr(colnames(tumor_trinuc_prim)[x], 1, 1)) & grepl("T", substr(colnames(tumor_trinuc_prim)[x], 7, 7)) == TRUE){
    tumor_trinuc_prim_avg_ordered <- c(tumor_trinuc_prim_avg_ordered, tumor_trinuc_prim_avg[x])
    
  }
}

#Reordering tumor_trinuc_met_avg
#N_A
for(x in 1:96){
  if (grepl("A", substr(colnames(tumor_trinuc_met)[x], 1, 1)) & grepl("A", substr(colnames(tumor_trinuc_met)[x], 7, 7)) == TRUE){
    tumor_trinuc_met_avg_ordered <- c(tumor_trinuc_met_avg_ordered, tumor_trinuc_met_avg[x])
  }
}
for(x in 1:96){
  if (grepl("C", substr(colnames(tumor_trinuc_met)[x], 1, 1)) & grepl("A", substr(colnames(tumor_trinuc_met)[x], 7, 7)) == TRUE){
    tumor_trinuc_met_avg_ordered <- c(tumor_trinuc_met_avg_ordered, tumor_trinuc_met_avg[x])
  }
}
for(x in 1:96){
  if (grepl("G", substr(colnames(tumor_trinuc_met)[x], 1, 1)) & grepl("A", substr(colnames(tumor_trinuc_met)[x], 7, 7)) == TRUE){
    tumor_trinuc_met_avg_ordered <- c(tumor_trinuc_met_avg_ordered, tumor_trinuc_met_avg[x])
  }
}
for(x in 1:96){
  if (grepl("T", substr(colnames(tumor_trinuc_met)[x], 1, 1)) & grepl("A", substr(colnames(tumor_trinuc_met)[x], 7, 7)) == TRUE){
    tumor_trinuc_met_avg_ordered <- c(tumor_trinuc_met_avg_ordered, tumor_trinuc_met_avg[x])
  }
}
#N_C
for(x in 1:96){
  if (grepl("A", substr(colnames(tumor_trinuc_met)[x], 1, 1)) & grepl("C", substr(colnames(tumor_trinuc_met)[x], 7, 7)) == TRUE){
    tumor_trinuc_met_avg_ordered <- c(tumor_trinuc_met_avg_ordered, tumor_trinuc_met_avg[x])
  }
}
for(x in 1:96){
  if (grepl("C", substr(colnames(tumor_trinuc_met)[x], 1, 1)) & grepl("C", substr(colnames(tumor_trinuc_met)[x], 7, 7)) == TRUE){
    tumor_trinuc_met_avg_ordered <- c(tumor_trinuc_met_avg_ordered, tumor_trinuc_met_avg[x])
  }
}
for(x in 1:96){
  if (grepl("G", substr(colnames(tumor_trinuc_met)[x], 1, 1)) & grepl("C", substr(colnames(tumor_trinuc_met)[x], 7, 7)) == TRUE){
    tumor_trinuc_met_avg_ordered <- c(tumor_trinuc_met_avg_ordered, tumor_trinuc_met_avg[x])
  }
}
for(x in 1:96){
  if (grepl("T", substr(colnames(tumor_trinuc_met)[x], 1, 1)) & grepl("C", substr(colnames(tumor_trinuc_met)[x], 7, 7)) == TRUE){
    tumor_trinuc_met_avg_ordered <- c(tumor_trinuc_met_avg_ordered, tumor_trinuc_met_avg[x])
  }
}
#N_G
for(x in 1:96){
  if (grepl("A", substr(colnames(tumor_trinuc_met)[x], 1, 1)) & grepl("G", substr(colnames(tumor_trinuc_met)[x], 7, 7)) == TRUE){
    tumor_trinuc_met_avg_ordered <- c(tumor_trinuc_met_avg_ordered, tumor_trinuc_met_avg[x])
  }
}
for(x in 1:96){
  if (grepl("C", substr(colnames(tumor_trinuc_met)[x], 1, 1)) & grepl("G", substr(colnames(tumor_trinuc_met)[x], 7, 7)) == TRUE){
    tumor_trinuc_met_avg_ordered <- c(tumor_trinuc_met_avg_ordered, tumor_trinuc_met_avg[x])
  }
}
for(x in 1:96){
  if (grepl("G", substr(colnames(tumor_trinuc_met)[x], 1, 1)) & grepl("G", substr(colnames(tumor_trinuc_met)[x], 7, 7)) == TRUE){
    tumor_trinuc_met_avg_ordered <- c(tumor_trinuc_met_avg_ordered, tumor_trinuc_met_avg[x])
  }
}
for(x in 1:96){
  if (grepl("T", substr(colnames(tumor_trinuc_met)[x], 1, 1)) & grepl("G", substr(colnames(tumor_trinuc_met)[x], 7, 7)) == TRUE){
    tumor_trinuc_met_avg_ordered <- c(tumor_trinuc_met_avg_ordered, tumor_trinuc_met_avg[x])
  }
}
#N_T
for(x in 1:96){
  if (grepl("A", substr(colnames(tumor_trinuc_met)[x], 1, 1)) & grepl("T", substr(colnames(tumor_trinuc_met)[x], 7, 7)) == TRUE){
    tumor_trinuc_met_avg_ordered <- c(tumor_trinuc_met_avg_ordered, tumor_trinuc_met_avg[x])
  }
}
for(x in 1:96){
  if (grepl("C", substr(colnames(tumor_trinuc_met)[x], 1, 1)) & grepl("T", substr(colnames(tumor_trinuc_met)[x], 7, 7)) == TRUE){
    tumor_trinuc_met_avg_ordered <- c(tumor_trinuc_met_avg_ordered, tumor_trinuc_met_avg[x])
  }
}
for(x in 1:96){
  if (grepl("G", substr(colnames(tumor_trinuc_met)[x], 1, 1)) & grepl("T", substr(colnames(tumor_trinuc_met)[x], 7, 7)) == TRUE){
    tumor_trinuc_met_avg_ordered <- c(tumor_trinuc_met_avg_ordered, tumor_trinuc_met_avg[x])
  }
}
for(x in 1:96){
  if (grepl("T", substr(colnames(tumor_trinuc_met)[x], 1, 1)) & grepl("T", substr(colnames(tumor_trinuc_met)[x], 7, 7)) == TRUE){
    tumor_trinuc_met_avg_ordered <- c(tumor_trinuc_met_avg_ordered, tumor_trinuc_met_avg[x])
    
  }
}

PRAD_trinuc_heatmap_data <- data.frame(mutation=trinuc.mutation_data$section_labels,
                                       Upstream=trinuc.mutation_data$Upstream,
                                       Downstream=trinuc.mutation_data$Downstream,
                                       mutated_from=trinuc.mutation_data$mutated_from,
                                       mutated_to=trinuc.mutation_data$mutated_to,
                                       trinuc_context=rep(NA, 96),
                                       trinuc_prim=tumor_trinuc_prim_avg_ordered,
                                       trinuc_met=tumor_trinuc_met_avg_ordered)

library(stringr)
PRAD_trinuc_heatmap_data$mutation <- str_replace(PRAD_trinuc_heatmap_data$mutation, "%->%", "\u2192")

for (x in 1:96){
  PRAD_trinuc_heatmap_data$trinuc_context[x] <- paste(PRAD_trinuc_heatmap_data$Upstream[x], 
                                                      PRAD_trinuc_heatmap_data$mutated_from[x],
                                                      PRAD_trinuc_heatmap_data$Downstream[x],
                                                      sep = "")  
}

levels(PRAD_trinuc_heatmap_data$mutation) <- c("C\u2192A", "C\u2192G", "C\u2192T", "T\u2192A", "T\u2192C", "T\u2192G")

# save(PRAD_trinuc_heatmap_data, file="PRAD_trinuc_heatmap_data.RData")
# 
# load("PRAD_trinuc_heatmap_data.RData")

PRAD_trinuc_bar_met <- ggplot(data=PRAD_trinuc_heatmap_data, aes(x = trinuc_context, y= trinuc_met*100, fill = mutation)) +
  geom_bar(stat = "identity") + theme_bw()+
  facet_wrap(.~mutation, nrow = 1, scales = "free_x")+
  theme(axis.line = element_line(color = 'white'), axis.ticks.y = element_blank(), panel.border = element_blank(),
        axis.text.y = element_text(), axis.text.x = element_text(angle=90, vjust = 0.5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank())+ 
  theme(legend.position="none")+
  xlab("") + ylab("Mutation type probability")+ ggtitle("Metastasis")+
  scale_y_continuous(labels = function(x) paste0(x, "%"))

PRAD_trinuc_bar_met
ggsave("PRAD_figures/PRAD_trinuc_bar_met.png", width=12, height=3)


########################################################################################################
#Early vs Late

threestage_final <- load_cesa("threestage_final.rds")
threestage_results <- snv_results(threestage_final)
threestage_results <- threestage_results$selection.1

#Extracting trinuc data

threestage_trinuc <- threestage_final@trinucleotide_mutation_weights$trinuc_proportion_matrix

#Seperating trinuc data into prim and met
tumor_names_threestage <- unique(threestage_final@maf$Unique_Patient_Identifier)

tumor_early <- c()
tumor_late <- c()

samples_threestage <- threestage_final@samples


for (x in tumor_names_threestage){
  if (samples_threestage$group[which(samples_threestage$Unique_Patient_Identifier == x)] == "Early"){
    tumor_early <- c(tumor_early, x)
  }
  else if (samples_threestage$group[which(samples_threestage$Unique_Patient_Identifier == x)] == "Late"){
    tumor_late <- c(tumor_late, x)
  }
}

tumor_trinuc_early <- data.frame(matrix(ncol = 96, nrow = 0))
tumor_trinuc_late <- data.frame(matrix(ncol = 96, nrow = 0))

for (x in tumor_early){
  tumor_trinuc_early <- rbind(tumor_trinuc_early, threestage_trinuc[x,])
}
colnames(tumor_trinuc_early) <- colnames(threestage_trinuc)

for(x in tumor_late){
  tumor_trinuc_late <- rbind(tumor_trinuc_late, threestage_trinuc[x,])
}
colnames(tumor_trinuc_late) <- colnames(threestage_trinuc)

#Averaging trinuc data
tumor_trinuc_early_avg <- apply(tumor_trinuc_early, 2, mean)
tumor_trinuc_late_avg <- apply(tumor_trinuc_late, 2, mean)

tumor_trinuc_early_avg_ordered <- c()
tumor_trinuc_late_avg_ordered <- c()


#Reordering tumor_trinuc_early_avg
#N_A
for(x in 1:96){
  if (grepl("A", substr(colnames(tumor_trinuc_early)[x], 1, 1)) & grepl("A", substr(colnames(tumor_trinuc_early)[x], 7, 7)) == TRUE){
    tumor_trinuc_early_avg_ordered <- c(tumor_trinuc_early_avg_ordered, tumor_trinuc_early_avg[x])
  }
}
for(x in 1:96){
  if (grepl("C", substr(colnames(tumor_trinuc_early)[x], 1, 1)) & grepl("A", substr(colnames(tumor_trinuc_early)[x], 7, 7)) == TRUE){
    tumor_trinuc_early_avg_ordered <- c(tumor_trinuc_early_avg_ordered, tumor_trinuc_early_avg[x])
  }
}
for(x in 1:96){
  if (grepl("G", substr(colnames(tumor_trinuc_early)[x], 1, 1)) & grepl("A", substr(colnames(tumor_trinuc_early)[x], 7, 7)) == TRUE){
    tumor_trinuc_early_avg_ordered <- c(tumor_trinuc_early_avg_ordered, tumor_trinuc_early_avg[x])
  }
}
for(x in 1:96){
  if (grepl("T", substr(colnames(tumor_trinuc_early)[x], 1, 1)) & grepl("A", substr(colnames(tumor_trinuc_early)[x], 7, 7)) == TRUE){
    tumor_trinuc_early_avg_ordered <- c(tumor_trinuc_early_avg_ordered, tumor_trinuc_early_avg[x])
  }
}
#N_C
for(x in 1:96){
  if (grepl("A", substr(colnames(tumor_trinuc_early)[x], 1, 1)) & grepl("C", substr(colnames(tumor_trinuc_early)[x], 7, 7)) == TRUE){
    tumor_trinuc_early_avg_ordered <- c(tumor_trinuc_early_avg_ordered, tumor_trinuc_early_avg[x])
  }
}
for(x in 1:96){
  if (grepl("C", substr(colnames(tumor_trinuc_early)[x], 1, 1)) & grepl("C", substr(colnames(tumor_trinuc_early)[x], 7, 7)) == TRUE){
    tumor_trinuc_early_avg_ordered <- c(tumor_trinuc_early_avg_ordered, tumor_trinuc_early_avg[x])
  }
}
for(x in 1:96){
  if (grepl("G", substr(colnames(tumor_trinuc_early)[x], 1, 1)) & grepl("C", substr(colnames(tumor_trinuc_early)[x], 7, 7)) == TRUE){
    tumor_trinuc_early_avg_ordered <- c(tumor_trinuc_early_avg_ordered, tumor_trinuc_early_avg[x])
  }
}
for(x in 1:96){
  if (grepl("T", substr(colnames(tumor_trinuc_early)[x], 1, 1)) & grepl("C", substr(colnames(tumor_trinuc_early)[x], 7, 7)) == TRUE){
    tumor_trinuc_early_avg_ordered <- c(tumor_trinuc_early_avg_ordered, tumor_trinuc_early_avg[x])
  }
}
#N_G
for(x in 1:96){
  if (grepl("A", substr(colnames(tumor_trinuc_early)[x], 1, 1)) & grepl("G", substr(colnames(tumor_trinuc_early)[x], 7, 7)) == TRUE){
    tumor_trinuc_early_avg_ordered <- c(tumor_trinuc_early_avg_ordered, tumor_trinuc_early_avg[x])
  }
}
for(x in 1:96){
  if (grepl("C", substr(colnames(tumor_trinuc_early)[x], 1, 1)) & grepl("G", substr(colnames(tumor_trinuc_early)[x], 7, 7)) == TRUE){
    tumor_trinuc_early_avg_ordered <- c(tumor_trinuc_early_avg_ordered, tumor_trinuc_early_avg[x])
  }
}
for(x in 1:96){
  if (grepl("G", substr(colnames(tumor_trinuc_early)[x], 1, 1)) & grepl("G", substr(colnames(tumor_trinuc_early)[x], 7, 7)) == TRUE){
    tumor_trinuc_early_avg_ordered <- c(tumor_trinuc_early_avg_ordered, tumor_trinuc_early_avg[x])
  }
}
for(x in 1:96){
  if (grepl("T", substr(colnames(tumor_trinuc_early)[x], 1, 1)) & grepl("G", substr(colnames(tumor_trinuc_early)[x], 7, 7)) == TRUE){
    tumor_trinuc_early_avg_ordered <- c(tumor_trinuc_early_avg_ordered, tumor_trinuc_early_avg[x])
  }
}
#N_T
for(x in 1:96){
  if (grepl("A", substr(colnames(tumor_trinuc_early)[x], 1, 1)) & grepl("T", substr(colnames(tumor_trinuc_early)[x], 7, 7)) == TRUE){
    tumor_trinuc_early_avg_ordered <- c(tumor_trinuc_early_avg_ordered, tumor_trinuc_early_avg[x])
  }
}
for(x in 1:96){
  if (grepl("C", substr(colnames(tumor_trinuc_early)[x], 1, 1)) & grepl("T", substr(colnames(tumor_trinuc_early)[x], 7, 7)) == TRUE){
    tumor_trinuc_early_avg_ordered <- c(tumor_trinuc_early_avg_ordered, tumor_trinuc_early_avg[x])
  }
}
for(x in 1:96){
  if (grepl("G", substr(colnames(tumor_trinuc_early)[x], 1, 1)) & grepl("T", substr(colnames(tumor_trinuc_early)[x], 7, 7)) == TRUE){
    tumor_trinuc_early_avg_ordered <- c(tumor_trinuc_early_avg_ordered, tumor_trinuc_early_avg[x])
  }
}
for(x in 1:96){
  if (grepl("T", substr(colnames(tumor_trinuc_early)[x], 1, 1)) & grepl("T", substr(colnames(tumor_trinuc_early)[x], 7, 7)) == TRUE){
    tumor_trinuc_early_avg_ordered <- c(tumor_trinuc_early_avg_ordered, tumor_trinuc_early_avg[x])
    
  }
}


#Reordering tumor_trinuc_late_avg
#N_A
for(x in 1:96){
  if (grepl("A", substr(colnames(tumor_trinuc_late)[x], 1, 1)) & grepl("A", substr(colnames(tumor_trinuc_late)[x], 7, 7)) == TRUE){
    tumor_trinuc_late_avg_ordered <- c(tumor_trinuc_late_avg_ordered, tumor_trinuc_late_avg[x])
  }
}
for(x in 1:96){
  if (grepl("C", substr(colnames(tumor_trinuc_late)[x], 1, 1)) & grepl("A", substr(colnames(tumor_trinuc_late)[x], 7, 7)) == TRUE){
    tumor_trinuc_late_avg_ordered <- c(tumor_trinuc_late_avg_ordered, tumor_trinuc_late_avg[x])
  }
}
for(x in 1:96){
  if (grepl("G", substr(colnames(tumor_trinuc_late)[x], 1, 1)) & grepl("A", substr(colnames(tumor_trinuc_late)[x], 7, 7)) == TRUE){
    tumor_trinuc_late_avg_ordered <- c(tumor_trinuc_late_avg_ordered, tumor_trinuc_late_avg[x])
  }
}
for(x in 1:96){
  if (grepl("T", substr(colnames(tumor_trinuc_late)[x], 1, 1)) & grepl("A", substr(colnames(tumor_trinuc_late)[x], 7, 7)) == TRUE){
    tumor_trinuc_late_avg_ordered <- c(tumor_trinuc_late_avg_ordered, tumor_trinuc_late_avg[x])
  }
}
#N_C
for(x in 1:96){
  if (grepl("A", substr(colnames(tumor_trinuc_late)[x], 1, 1)) & grepl("C", substr(colnames(tumor_trinuc_late)[x], 7, 7)) == TRUE){
    tumor_trinuc_late_avg_ordered <- c(tumor_trinuc_late_avg_ordered, tumor_trinuc_late_avg[x])
  }
}
for(x in 1:96){
  if (grepl("C", substr(colnames(tumor_trinuc_late)[x], 1, 1)) & grepl("C", substr(colnames(tumor_trinuc_late)[x], 7, 7)) == TRUE){
    tumor_trinuc_late_avg_ordered <- c(tumor_trinuc_late_avg_ordered, tumor_trinuc_late_avg[x])
  }
}
for(x in 1:96){
  if (grepl("G", substr(colnames(tumor_trinuc_late)[x], 1, 1)) & grepl("C", substr(colnames(tumor_trinuc_late)[x], 7, 7)) == TRUE){
    tumor_trinuc_late_avg_ordered <- c(tumor_trinuc_late_avg_ordered, tumor_trinuc_late_avg[x])
  }
}
for(x in 1:96){
  if (grepl("T", substr(colnames(tumor_trinuc_late)[x], 1, 1)) & grepl("C", substr(colnames(tumor_trinuc_late)[x], 7, 7)) == TRUE){
    tumor_trinuc_late_avg_ordered <- c(tumor_trinuc_late_avg_ordered, tumor_trinuc_late_avg[x])
  }
}
#N_G
for(x in 1:96){
  if (grepl("A", substr(colnames(tumor_trinuc_late)[x], 1, 1)) & grepl("G", substr(colnames(tumor_trinuc_late)[x], 7, 7)) == TRUE){
    tumor_trinuc_late_avg_ordered <- c(tumor_trinuc_late_avg_ordered, tumor_trinuc_late_avg[x])
  }
}
for(x in 1:96){
  if (grepl("C", substr(colnames(tumor_trinuc_late)[x], 1, 1)) & grepl("G", substr(colnames(tumor_trinuc_late)[x], 7, 7)) == TRUE){
    tumor_trinuc_late_avg_ordered <- c(tumor_trinuc_late_avg_ordered, tumor_trinuc_late_avg[x])
  }
}
for(x in 1:96){
  if (grepl("G", substr(colnames(tumor_trinuc_late)[x], 1, 1)) & grepl("G", substr(colnames(tumor_trinuc_late)[x], 7, 7)) == TRUE){
    tumor_trinuc_late_avg_ordered <- c(tumor_trinuc_late_avg_ordered, tumor_trinuc_late_avg[x])
  }
}
for(x in 1:96){
  if (grepl("T", substr(colnames(tumor_trinuc_late)[x], 1, 1)) & grepl("G", substr(colnames(tumor_trinuc_late)[x], 7, 7)) == TRUE){
    tumor_trinuc_late_avg_ordered <- c(tumor_trinuc_late_avg_ordered, tumor_trinuc_late_avg[x])
  }
}
#N_T
for(x in 1:96){
  if (grepl("A", substr(colnames(tumor_trinuc_late)[x], 1, 1)) & grepl("T", substr(colnames(tumor_trinuc_late)[x], 7, 7)) == TRUE){
    tumor_trinuc_late_avg_ordered <- c(tumor_trinuc_late_avg_ordered, tumor_trinuc_late_avg[x])
  }
}
for(x in 1:96){
  if (grepl("C", substr(colnames(tumor_trinuc_late)[x], 1, 1)) & grepl("T", substr(colnames(tumor_trinuc_late)[x], 7, 7)) == TRUE){
    tumor_trinuc_late_avg_ordered <- c(tumor_trinuc_late_avg_ordered, tumor_trinuc_late_avg[x])
  }
}
for(x in 1:96){
  if (grepl("G", substr(colnames(tumor_trinuc_late)[x], 1, 1)) & grepl("T", substr(colnames(tumor_trinuc_late)[x], 7, 7)) == TRUE){
    tumor_trinuc_late_avg_ordered <- c(tumor_trinuc_late_avg_ordered, tumor_trinuc_late_avg[x])
  }
}
for(x in 1:96){
  if (grepl("T", substr(colnames(tumor_trinuc_late)[x], 1, 1)) & grepl("T", substr(colnames(tumor_trinuc_late)[x], 7, 7)) == TRUE){
    tumor_trinuc_late_avg_ordered <- c(tumor_trinuc_late_avg_ordered, tumor_trinuc_late_avg[x])
    
  }
}

threestage_trinuc_heatmap_data <- data.frame(mutation=trinuc.mutation_data$mutation,
                                             Upstream=trinuc.mutation_data$Upstream,
                                             Downstream=trinuc.mutation_data$Downstream,
                                             mutated_from=trinuc.mutation_data$mutated_from,
                                             mutated_to=trinuc.mutation_data$mutated_to,
                                             trinuc_context=rep(NA, 96),
                                             trinuc_early=tumor_trinuc_early_avg_ordered,
                                             trinuc_late=tumor_trinuc_late_avg_ordered)


library(stringr)
threestage_trinuc_heatmap_data$mutation <- str_replace(threestage_trinuc_heatmap_data$mutation, "to", "\u2192")

for (x in 1:96){
  threestage_trinuc_heatmap_data$trinuc_context[x] <- paste(threestage_trinuc_heatmap_data$Upstream[x], 
                                                            threestage_trinuc_heatmap_data$mutated_from[x],
                                                            threestage_trinuc_heatmap_data$Downstream[x],
                                                            sep = "")  
}

levels(threestage_trinuc_heatmap_data$mutation) <- c("C\u2192A", "C\u2192G", "C\u2192T", "T\u2192A", "T\u2192C", "T\u2192G")

# save(threestage_trinuc_heatmap_data, file="threestage_trinuc_heatmap_data.RData")
# 
# load("threestage_trinuc_heatmap_data.RData")


PRAD_trinuc_bar_early <- ggplot(data=threestage_trinuc_heatmap_data, aes(x = trinuc_context, y= trinuc_early*100, fill = mutation)) +
  geom_bar(stat = "identity") + theme_bw()+
  facet_wrap(.~mutation, nrow = 1, scales = "free_x")+
  theme(axis.line = element_line(color = 'white'), axis.ticks.y = element_blank(), panel.border = element_blank(),
        axis.text.y = element_text(), axis.text.x = element_text(angle=90, vjust = 0.5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank())+ 
  theme(legend.position="none")+
  xlab("") + ylab("Mutation type probability")+ ggtitle("Lower-risk")+
  scale_y_continuous(labels = function(x) paste0(x, "%"))

PRAD_trinuc_bar_early

PRAD_trinuc_bar_late <- ggplot(data=threestage_trinuc_heatmap_data, aes(x = trinuc_context, y= trinuc_late*100, fill = mutation)) +
  geom_bar(stat = "identity") + theme_bw()+
  facet_wrap(.~mutation, nrow = 1, scales = "free_x")+
  theme(axis.line = element_line(color = 'white'), axis.ticks.y = element_blank(), panel.border = element_blank(),
        axis.text.y = element_text(), axis.text.x = element_text(angle=90, vjust = 0.5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank())+ 
  theme(legend.position="none")+
  xlab("") + ylab("Mutation type probability")+ ggtitle("Higher-risk")+
  scale_y_continuous(labels = function(x) paste0(x, "%"))

PRAD_trinuc_bar_late

#Combined graphs

combined_trinuc_bar <-plot_grid(PRAD_trinuc_bar_early, PRAD_trinuc_bar_late, PRAD_trinuc_bar_met,
                                labels = c("A", "B", "C"), label_size = 12,
                                align="h", axis="t", nrow=3, ncol=1, rel_heights = c(1,1,1))

combined_trinuc_bar

ggsave("PRAD_figures/combined_trinuc_bar.png", width=12, height=8)
