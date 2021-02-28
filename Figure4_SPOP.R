library(cancereffectsizeR)
library(scales)
library(stringr)
library(ggplot2)

scientific <- function(x){
  ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", label_scientific()(x)))))
}

common.text.size <- 4
###############

twostage_final <- load_cesa("twostage_final.rds")

PRAD_analysis <- twostage_final
PRAD_results <- snv_results(PRAD_analysis)
PRAD_results <- PRAD_results$selection.1

gene_name <- word(PRAD_results$variant_name, sep = "_")
PRAD_results$gene <- gene_name

PRAD_results$variant_name <- str_replace(PRAD_results$variant_name, "SPOP_", "")

aac <- PRAD_results$variant_type == "aac"
PRAD_results <- PRAD_results[aac,]

recurrent <- PRAD_results$maf_freq_in_Primary > 1
PRAD_results_recurrent <- PRAD_results[recurrent,]

all <- PRAD_results$maf_freq_in_Primary > 0
PRAD_results <- PRAD_results[all,]

PRAD_results <- PRAD_results[order(-si_1),]
PRAD_results_recurrent <- PRAD_results_recurrent[order(-si_1),]

rm(twostage_final)

#Selection Intesities for SPOP in Prim

spop_true <- PRAD_results$gene == "SPOP"
PRAD_results <- PRAD_results[spop_true,]

spop_true <- PRAD_results_recurrent$gene == "SPOP"
PRAD_results_recurrent <- PRAD_results_recurrent[spop_true,]

#########################################################################

bargraph_SPOP_SI <- ggplot(data=PRAD_results_recurrent, aes(x=reorder(variant_name, -si_1), y=si_1, fill=si_1))+
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=ci_low_95_si_1, ymax=ci_high_95_si_1), width=0.333) +
  theme(axis.text.x = element_text (hjust = 1, angle = 45)) +
  ylab("Scaled selection coefficient") +
  scale_fill_gradient(low="gold", high="red2") +
  theme(legend.position = "none")+
  theme(panel.background = element_blank(), axis.title.x = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  geom_text(aes(label=maf_freq_in_Primary, y=-5000), size = common.text.size) +
  scale_y_continuous(labels=scientific, breaks = c(0, 5e4, 1e5, 1.5e5, 2e5))

bargraph_SPOP_SI
ggsave("PRAD_figures/SPOP/SPOP_recurrent_SI.png", width=9, height=6)