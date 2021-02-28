library(cancereffectsizeR)
library(scales)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(RColorBrewer)

scientific <- function(x){
  ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", label_scientific()(x)))))
}

# twostage_final <- load_cesa("twostage_final.rds")
# 
# mut_rate <- data.frame(gene=twostage_final@mutrates$gene,
#                                  prim_rate=twostage_final@mutrates$rate_grp_1,
#                                  met_rate=twostage_final@mutrates$rate_grp_2)
# 
# save(mut_rate, file="PRAD_mutrate.RData")

load("PRAD_mutrate.RData")

# threestage_final <- load_cesa("threestage_final.rds")
# 
# mut_rate_earlylate <- data.frame(gene=threestage_final@mutrates$gene,
#                        early_rate=threestage_final@mutrates$rate_grp_1,
#                        late_rate=threestage_final@mutrates$rate_grp_2)
# 
# save(mut_rate_earlylate, file="PRAD_mutrate_earlylate.RData")

load("PRAD_mutrate_earlylate.RData")
# 
# 
# PRAD_analysis_stageless <- load_cesa("stageless_final.rds")
# 
# mut_rate_stageless <- data.frame(gene=PRAD_analysis_stageless@mutrates$gene,
#                                  mutation_rate=PRAD_analysis_stageless@mutrates$rate_grp_1)
# 
# save(mut_rate_stageless, file="PRAD_mutrate_stageless.RData")

load("PRAD_mutrate_stageless.RData")


#applicable to both stageless and prim vs met
highlight <- rep(FALSE, length(mut_rate_stageless$gene))
highlight[c(16063, 6373, 1221, 12567, 17456)] <- TRUE

#####################################################################################################################

#scatter plot of gene-level mutation rates, prim vs met
selected_colors <- brewer.pal(n = length(highlight[highlight==TRUE]), name = 'Dark2')
selected_colors[2] <- "#027cd9"


mutrates_stageless_plot <- ggplot() +
  geom_jitter(data = mut_rate_stageless[!highlight,], aes(x=1, y= mutation_rate, color = mutation_rate), shape=16, position=position_jitter(0.05), size = .75, alpha=0.6) +
  scale_colour_gradient(low="#FF9A9A", high="#F08080") +
  ylab("Mutation rate") + xlab("") +
  geom_jitter(data = mut_rate_stageless[highlight,], aes(x=1, y= mutation_rate), shape=16, position=position_jitter(0.05, seed = 5), size = 3, color=selected_colors, alpha=1) +
  geom_text_repel(data = mut_rate_stageless[highlight,], aes(x=1, y= mutation_rate, label = gene), position=position_jitter(0.05, seed = 5)) +
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(color = "black"),
        legend.position="none") +
  theme(plot.margin = margin(0,0,0,45, "pt")) +
  coord_flip() + scale_y_continuous(labels=scientific) +
  expand_limits(y=c(0, 2.25e-6))

mutrates_prim_met_plot <- ggplot()+
  geom_point(data = mut_rate[!highlight,], aes(x=prim_rate, y=met_rate), size=1, color="lightcoral", alpha=0.5, shape = 16) +
  geom_point(data = mut_rate[highlight,], aes(x=prim_rate, y=met_rate), size=2.5, color=selected_colors, alpha=1) +
  geom_text_repel(data = mut_rate[highlight,], aes(x=prim_rate, y=met_rate, label=gene)) +
  theme_bw() +
  theme(panel.border = element_blank(), plot.title = element_text(hjust=0.5)) + 
  theme(panel.border = element_blank(),
        axis.line = element_line(color = 'black'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  xlab("Mutation rate in primary tumors") +
  ylab("Mutation rate in metastatic tumors") +
  geom_smooth(method="lm", color="navyblue") +
  geom_abline(slope=1, intercept=0, color = "darkred", linetype = "dashed") +
  expand_limits(x = c(0, 1.25e-6)) + 
  scale_x_continuous(labels=scientific, breaks=c(0, 5e-7, 1e-6)) + 
  scale_y_continuous(labels=scientific, breaks=c(0, 1e-6, 2e-6))

mutrates_early_late_plot <- ggplot()+
  geom_point(data = mut_rate_earlylate[!highlight,], aes(x=early_rate, y=late_rate), size=1, color="lightcoral", alpha=0.5, shape = 16) +
  geom_point(data = mut_rate_earlylate[highlight,], aes(x=early_rate, y=late_rate), size=2.5, color=selected_colors, alpha=1) +
  geom_text_repel(data = mut_rate_earlylate[highlight,], aes(x=early_rate, y=late_rate, label=gene)) +
  theme_bw() +
  theme(panel.border = element_blank(), plot.title = element_text(hjust=0.5)) + 
  theme(panel.border = element_blank(),
        axis.line = element_line(color = 'black'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  xlab("Mutation rate in lower-risk tumors") +
  ylab("Mutation rate in higher-risk tumors") +
  geom_smooth(method="lm", color="navyblue") + 
  geom_abline(slope=1, intercept=0, color = "darkred", linetype = "dashed") +
  scale_x_continuous(labels=scientific, limits=c(0,1.25e-6), breaks=c(0, 5e-7, 1e-6)) +
  scale_y_continuous(labels=scientific, breaks=c(0, 1e-6, 2e-6))

library(cowplot)


combined_gene_mutrate <-plot_grid(mutrates_stageless_plot, mutrates_early_late_plot, mutrates_prim_met_plot, 
                                  labels = c("A", "B", "C"), label_size = 12,
                                  align="h", axis="t", nrow=3, ncol=1, rel_heights = c(0.3,1,1))

combined_gene_mutrate
