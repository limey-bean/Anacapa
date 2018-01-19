library(ranacapa)
library(dplyr)
library(tibble)
library(reshape2)
library(phyloseq)
library(ggplot2)
library(vegan)
library(plotly)
library(optparse)
args <- commandArgs(trailingOnly=TRUE)


input_biom_path <- args[1]
input_meta_path <- args[2]

anacapa_biom <- read.table(input_biom_path, header = 1, sep = "\t", stringsAsFactors = F)
metadata <- read.table(input_meta_path, header = 1, sep = "\t", stringsAsFactors = F)
output_path <- args[3]

# anacapa_biom <- readRDS("~/grad/ranacapa/inst/shiny-examples/explore-anacapa-output/data/demo_biomdata.Rds")
# metadata <- readRDS("~/grad/ranacapa/inst/shiny-examples/explore-anacapa-output/data/demo_metadata.Rds")
# output_path <- "~/Desktop/test-ranacapa-graphs"


if(!dir.exists(output_path)) {
  dir.create(output_path)
}

physeq_obj <- convert_anacapa_to_phyloseq(anacapa_biom, metadata)


# Perform rarefaction and make grapsh ----------------
p <- ggrare(physeq_obj,color = colnames(metadata[1]),step = 1000, se=FALSE) + theme_bw() + theme_ranacapa()
ggsave(plot = p, filename = file.path(output_path, "rarefaction_overall.pdf"))

rarefaction_depth = 2000
rarefaction_reps  = 2

physeq_obj_rare <- custom_rarefaction(physeq_obj, sample_size = rarefaction_depth, replicates = rarefaction_reps)

p <- ggrare(physeq_obj_rare, color = colnames(metadata[1]),step = 1000, se=FALSE) + theme_bw() + theme_ranacapa()
ggsave(plot = p, filename = file.path(output_path, "rarefaction_output.pdf"))



# Alpha diversity graphs and Anova tables -----------
sink(file = file.path(output_path, paste0("AlphaDiv_AnovaSummaries.txt")), append = F)
for(current_variable in colnames(metadata)) {
  p <- plot_richness(physeq_obj_rare, x = (current_variable), measures = c("Observed", "Shannon")) +
    geom_boxplot(aes_string(fill = current_variable, alpha=0.2, show.legend = F)) + theme_bw() +
    theme_ranacapa()
  ggsave(plot = p, file = file.path(output_path, paste0("AlphaDiv_by_", current_variable,".pdf")))
  
  alpha.diversity <- estimate_richness(physeq_obj_rare, measures = c("Observed", "Shannon"))
  data <- cbind(sample_data(physeq_obj_rare), alpha.diversity)
  cat(paste("Variable:", current_variable,"\n"))
  cat("Observed:\n")
  obj <- aov(as.formula(paste("Observed", "~" , current_variable)), data)
  print(broom::tidy(obj))
  print(broom::tidy(TukeyHSD(obj)))
  cat("\nShannon Diversity:\n")
  obj <- aov(as.formula(paste("Shannon", "~" , current_variable)), data)
  print(broom::tidy(obj))
  print(broom::tidy(TukeyHSD(obj)))
  cat("\n\n----------------\n\n")
}
sink()


# Beta diversity begins ------------

# Jaccard

# Beta diversity based on Jaccard
dissimMethod = "jaccard"
d <- distance(physeq_obj_rare, method=dissimMethod)
sampledf <- data.frame(sample_data(physeq_obj_rare))
veganComm <- vegan_otu(physeq_obj_rare)
ord <- ordinate(physeq_obj_rare, method = "MDS", distance = d)

# Big loop across all variables
sink(file = file.path(output_path, paste0(dissimMethod,"_beta_statistics.txt")))
for (current_variable in colnames(metadata)){
  
  # NMDS ordination
  
  nmdsplot <- plot_ordination(physeq_obj_rare, ord, current_variable, color = current_variable) +
    theme_bw() +  # stat_ellipse(type = "t", geom = "polygon", alpha = 0.2) +
    theme_ranacapa() + ggtitle(paste(current_variable, "NMDS; dissimilarity method:",
                                     tools::toTitleCase(dissimMethod))) 
  ggsave(plot = nmdsplot, file = file.path(output_path, paste0(dissimMethod,"_nmds_by_",current_variable,".pdf")))
  
  # Network
  ig <- make_network(physeq_obj_rare,
                     distance=function(x){vegan::vegdist(x, dissimMethod)}, max.dist=0.9)
  
  igp <-  plot_network(ig, physeq_obj_rare, color=current_variable, line_weight=0.9, label=NULL) +
    theme_bw(base_size = 18)
  ggsave(plot = igp, file = file.path(output_path,paste0(dissimMethod,"network_by_",current_variable,".pdf")))
  
  # Beta diversity stats 
  cat(paste("Adonis tables for", current_variable,"\n"))
  print((adonis(as.formula(paste("d~", current_variable)), data = sampledf)$aov.tab))
  print(pairwise.adonis(veganComm,getElement(sampledf, current_variable),sim.method = dissimMethod))
  cat("\n\n")
  
  cat("Multivariate homogeneity of groups dispersions for ", current_variable,":\n")
  print(betadisper(d, getElement(sampledf, current_variable)))
  print(broom::tidy(TukeyHSD(betadisper(d, getElement(sampledf, current_variable)))))
  cat("\n\n ------------------ \n\n")
}
sink()

wcluster <- as.dendrogram(hclust(d, method = "ward.D2"))

wl <- ggdendro::ggdendrogram(wcluster, theme_dendro = FALSE, color = "red")  +
  theme_bw(base_size = 18) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(plot = wl, file = file.path(output_path, paste0(dissimMethod,"_ward_linkage.pdf")))

heat <- ggplot(data = melt(as.matrix(d)), aes(x=Var1, y=Var2, fill=value)) +
  geom_tile() + theme_bw(base_size = 18) + theme_ranacapa() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(plot = heat, file = file.path(output_path, paste0(dissimMethod,"_heatmap.pdf")))

# Bray

dissimMethod = "bray"

sink(file = file.path(output_path, paste0(dissimMethod,"_beta_statistics.txt")))
for (current_variable in colnames(metadata)){
  
  # NMDS ordination
  
  nmdsplot <- plot_ordination(physeq_obj_rare, ord, current_variable, color = current_variable) +
    theme_bw() +  # stat_ellipse(type = "t", geom = "polygon", alpha = 0.2) +
    theme_ranacapa() + ggtitle(paste(current_variable, "NMDS; dissimilarity method:",
                                     tools::toTitleCase(dissimMethod))) 
  ggsave(plot = nmdsplot, file = file.path(output_path, paste0(dissimMethod,"_nmds_by_",current_variable,".pdf")))
  
  # Network
  ig <- make_network(physeq_obj_rare,
                     distance=function(x){vegan::vegdist(x, dissimMethod)}, max.dist=0.9)
  
  igp <-  plot_network(ig, physeq_obj_rare, color=current_variable, line_weight=0.9, label=NULL) +
    theme_bw(base_size = 18)
  ggsave(plot = igp, file = file.path(output_path,paste0(dissimMethod,"_network_by_",current_variable,".pdf")))
  
  # Beta diversity stats 
  cat(paste("Adonis tables for", current_variable,"\n"))
  print((adonis(as.formula(paste("d~", current_variable)), data = sampledf)$aov.tab))
  print(pairwise.adonis(veganComm,getElement(sampledf, current_variable),sim.method = dissimMethod))
  cat("\n\n")
  
  cat("Multivariate homogeneity of groups dispersions for ", current_variable,":\n")
  print(betadisper(d, getElement(sampledf, current_variable)))
  print(broom::tidy(TukeyHSD(betadisper(d, getElement(sampledf, current_variable)))))
  cat("\n\n ------------------ \n\n")
}
sink()