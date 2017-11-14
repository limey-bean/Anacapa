## Load Necessary Files
library("qiimer", lib.loc="/Library/Frameworks/R.framework/Versions/3.3/Resources/library")
library(reshape2)
library(ggplot2)
library(ggrepel)
library(vegan)
library(dplyr)
library(wesanderson) # I just want to give a shout out to karthik for being a true legend and fellow Wes Anderson afficiando, as well as to Gaurav for sharing with me this amazing package. Truly life altering
library(cluster)
library(biomformat)
library('phyloseq')
library('gridBase')
library("optparse")


##Source https://github.com/mahendra-mariadassou/phyloseq-extended

# source("/Users/zackgold/Documents/UCLA_phd/phyloseq-extended-master/graphical_methods.R")
source("https://raw.githubusercontent.com/mahendra-mariadassou/phyloseq-extended/master/graphical_methods.R")
# source("/Users/zackgold/Documents/UCLA_phd/phyloseq-extended-master/tree_methods.R")
source("https://raw.githubusercontent.com/mahendra-mariadassou/phyloseq-extended/master/tree_methods.R")
# source("/Users/zackgold/Documents/UCLA_phd/phyloseq-extended-master/plot_merged_trees.R")
source("https://raw.githubusercontent.com/mahendra-mariadassou/phyloseq-extended/master/plot_merged_trees.R")
# source("/Users/zackgold/Documents/UCLA_phd/phyloseq-extended-master/specificity_methods.R")
source("https://raw.githubusercontent.com/mahendra-mariadassou/phyloseq-extended/master/specificity_methods.R")
# source("/Users/zackgold/Documents/UCLA_phd/phyloseq-extended-master/ternary_plot.R")
source("https://raw.githubusercontent.com/mahendra-mariadassou/phyloseq-extended/master/ternary_plot.R")
# source("/Users/zackgold/Documents/UCLA_phd/phyloseq-extended-master/richness.R")
source("https://raw.githubusercontent.com/mahendra-mariadassou/phyloseq-extended/master/richness.R")
# source("/Users/zackgold/Documents/UCLA_phd/phyloseq-extended-master/edgePCA.R")
source("https://raw.githubusercontent.com/mahendra-mariadassou/phyloseq-extended/master/edgePCA.R")
# source("/Users/zackgold/Documents/UCLA_phd/phyloseq-extended-master/copy_number_correction.R")
source("https://raw.githubusercontent.com/mahendra-mariadassou/phyloseq-extended/master/copy_number_correction.R")
# source("/Users/zackgold/Documents/UCLA_phd/phyloseq-extended-master/import_frogs.R")
source("https://raw.githubusercontent.com/mahendra-mariadassou/phyloseq-extended/master/copy_number_correction.R")
# source("/Users/zackgold/Documents/UCLA_phd/phyloseq-extended-master/prevalence.R")
source("https://raw.githubusercontent.com/mahendra-mariadassou/phyloseq-extended/master/prevalence.R")
# source("/Users/zackgold/Documents/UCLA_phd/phyloseq-extended-master/compute_niche.R")
source("https://raw.githubusercontent.com/mahendra-mariadassou/phyloseq-extended/master/compute_niche.R")
# source('/Users/zackgold/Documents/UCLA_phd/phyloseq-extended-master/pairwise_adonis.R')




#Define Arguments Needed for R script

args = commandArgs(trailingOnly=TRUE)


# Manage Null Arguments
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[3] = "out.txt"
}

#Import Biom Table
## This Biom table should have already been decontaminated and rarified through the Anacapa pipeline

#File Path to Biom Table
file_path <- args[1]
#file_path <- "/Users/zackgold/Documents/UCLA_phd/Projects/anacapa/march_seq_test/catalina_test/12s/moorea_kraken/12S_F_moorea_rarefied_6000/summarize_taxa_mor/rarefaction_6000__metadata_json_fixed.biom"

#Read Biom Table
dat <- read_biom(file_path)

#Prepare data for PhyloSeq
#Create otu table
otu <- as.matrix(biom_data(dat))

#Change OTU the column names

###Need to make changes to the R code to handle any sample list
####This problem is due to stupid naming scheme that start with numbers
##Specific to Zack's Mo'orea because of R naming convention
colnames(otu) <- c("DeepRF2.","PN1.","BRF3.","KM1_1.","LB2.","DeepRF1.","BRF1.", "PN2.","TMI2.","TMI3.", "TMI1.","PN3.","KM4_2.","M500_2.","KM1_3.","KM2_1.","MRB1.","M500_3.","LF1.","BRF2.","LB3.","KM2_2.","MidRF1.","KM4_3.","ShallowRF3.","MidRF3.","LF2.","MidRF2.","M500_1.","MRB2.","ShallowRF1.","Mor.poscntl.","MRB3.","DeepRF3.","LB1.","LF3.","KM1_2.","ShallowRF2.","KM2_3.","KM4_1.")   
#order samples alphabetically
otu<- otu[ , order(colnames(otu))]
#Phyloseq OTU type matrix
OTU = otu_table(otu, taxa_are_rows = TRUE)

#Create a taxonomy table
tax_phy <- as.matrix(observation_metadata(dat))
TAX = tax_table(tax_phy)
#Change the TAX column names
colnames(TAX) <- c("Domain","Kingdom","Phylum","Class","Order","Family","Genus","Species")

#MetaData
#Import metadata
#args[2]="/Users/zackgold/Documents/UCLA_phd/Projects/anacapa/input/moorea_1_seq_map__Project_Moorea__updated_again_3.txt"

#Convert to sample_data type file for phyloseq
mor_map <- read_qiime_mapping_file(args[2])
mor_map <- as.data.frame(mor_map)
rownames(mor_map) <- mor_map[,1]
#order samples alphabetically
mor_map <- mor_map[order(mor_map$SampleID),]
sampledata <- sample_data(mor_map)

#Build Phyloseq dataframe structure
physeq = phyloseq(OTU, TAX)
physeq1 = merge_phyloseq(physeq, sampledata)
physeq1

#Remove Columns to Analyze
heads <- colnames(sampledata)
heads <- heads[-1:-5]
heads

## Rarefaction curve, ggplot style
ggrare <- function(physeq, step = 10, label = NULL, color = NULL, plot = TRUE, parallel = FALSE, se = TRUE) {
  ## Args:
  ## - physeq: phyloseq class object, from which abundance data are extracted
  ## - step: Step size for sample size in rarefaction curves
  ## - label: Default `NULL`. Character string. The name of the variable
  ##          to map to text labels on the plot. Similar to color option
  ##          but for plotting text.
  ## - color: (Optional). Default ‘NULL’. Character string. The name of the
  ##          variable to map to colors in the plot. This can be a sample
  ##          variable (among the set returned by
  ##          ‘sample_variables(physeq)’ ) or taxonomic rank (among the set
  ##          returned by ‘rank_names(physeq)’).
  ##
  ##          Finally, The color scheme is chosen automatically by
  ##          ‘link{ggplot}’, but it can be modified afterward with an
  ##          additional layer using ‘scale_color_manual’.
  ## - color: Default `NULL`. Character string. The name of the variable
  ##          to map to text labels on the plot. Similar to color option
  ##          but for plotting text.
  ## - plot:  Logical, should the graphic be plotted.
  ## - parallel: should rarefaction be parallelized (using parallel framework)
  ## - se:    Default TRUE. Logical. Should standard errors be computed. 
  ## require vegan
  x <- as(otu_table(physeq), "matrix")
  if (taxa_are_rows(physeq)) { x <- t(x) }
  
  ## This script is adapted from vegan `rarecurve` function
  tot <- rowSums(x)
  S <- rowSums(x > 0)
  nr <- nrow(x)
  
  rarefun <- function(i) {
    cat(paste("rarefying sample", rownames(x)[i]), sep = "\n")
    n <- seq(1, tot[i], by = step)
    if (n[length(n)] != tot[i]) {
      n <- c(n, tot[i])
    }
    y <- rarefy(x[i, ,drop = FALSE], n, se = se)
    if (nrow(y) != 1) {
      rownames(y) <- c(".S", ".se")
      return(data.frame(t(y), Size = n, Sample = rownames(x)[i]))
    } else {
      return(data.frame(.S = y[1, ], Size = n, Sample = rownames(x)[i]))
    }
  }
  if (parallel) {
    out <- mclapply(seq_len(nr), rarefun, mc.preschedule = FALSE)
  } else {
    out <- lapply(seq_len(nr), rarefun)
  }
  df <- do.call(rbind, out)
  
  ## Get sample data 
  if (!is.null(sample_data(physeq, FALSE))) {
    sdf <- as(sample_data(physeq), "data.frame")
    sdf$Sample <- rownames(sdf)
    data <- merge(df, sdf, by = "Sample")
    labels <- data.frame(x = tot, y = S, Sample = rownames(x))
    labels <- merge(labels, sdf, by = "Sample")
  }
  
  ## Add, any custom-supplied plot-mapped variables
  if( length(color) > 1 ){
    data$color <- color
    names(data)[names(data)=="color"] <- deparse(substitute(color))
    color <- deparse(substitute(color))
  }
  if( length(label) > 1 ){
    labels$label <- label
    names(labels)[names(labels)=="label"] <- deparse(substitute(label))
    label <- deparse(substitute(label))
  }
  
  p <- ggplot(data = data, aes_string(x = "Size", y = ".S", group = "Sample", color = color))
  p <- p + labs(x = "Sequence Sample Size", y = "Species Richness")
  if (!is.null(label)) {
    p <- p + geom_text(data = labels, aes_string(x = "x", y = "y", label = label, color = color),
                       size = 4, hjust = 0)
  }
  p <- p + geom_line()
  if (se) { ## add standard error if available
    p <- p + geom_ribbon(aes_string(ymin = ".S - .se", ymax = ".S + .se", color = NULL, fill = color), alpha = 0.2)
  }
  if (plot) {
    plot(p)
  }
  invisible(p)
}




##Alpha Diversity
#Rarefaction Curves
for (i in 1:length(heads)){
  #subset data for non NA
  physeq1_notna <- physeq1 %>% subset_samples(!is.na(get_variable(physeq1, heads[i])))
  #Conduct rarefaction of samples
  p <- ggrare(physeq1_notna, step = 1000, color = heads[i], se=FALSE)
  
  #Plot rarefaction Curve
  head_clean <- as.name(heads[i])
  p <- p + facet_wrap(as.formula(paste("~", head_clean)), ncol = 4)
  p <- p + theme_bw() + theme(panel.grid.minor.y=element_blank(),panel.grid.minor.x=element_blank(),panel.grid.major.y=element_blank(),panel.grid.major.x=element_blank())
  
  #Save Plots
  name <- paste0(head_clean,"_rarefaction.png")
  ggsave(name, width = 20, height = 10)
}


#Alpha Diversity Boxplots
#currently only focussed on Observed and Shannon, could alter for any type of diversity analysis
for (i in 1:length(heads)){
  #subset data for non NA
  head_clean <- as.name(heads[i])
  head_clean
  physeq1_notna <- physeq1 %>% subset_samples(!is.na(get_variable(physeq1, heads[i])))
  
  #Plot Richness Box plot
  richness <- plot_richness(physeq1_notna, x=heads[i], color=heads[i], measures= c("Observed", "Shannon")) + geom_boxplot(aes(fill=paste(head_clean)),alpha=0.2) +theme_bw() +  theme(panel.grid.minor.y=element_blank(),panel.grid.minor.x=element_blank(),panel.grid.major.y=element_blank(),panel.grid.major.x=element_blank())
  richness
  
  #Save Plot
  name <- paste0(head_clean,"_richness.png")
  ggsave(name,width = 20, height = 10)
}

#Alpha Diversity Statistics
alpha.diversity <- estimate_richness(physeq1, measures = c("Observed", "Shannon"))
data <- cbind(sample_data(physeq1), alpha.diversity)

#Observed Statistics

for (i in 1:length(heads)){
  sink("Observerd_Species_Diversity_Statistics.txt", append=TRUE)
  head_clean <- as.name(heads[i])
  print(head_clean)
  physeq1.alpha.anova <- aov(as.formula(paste("Observed~", head_clean)), data)
  log_1 <- summary(physeq1.alpha.anova )
  log_2<-TukeyHSD(physeq1.alpha.anova)
  print(log_1)
  print(log_2)
  cat("=============================\n")
  sink()
}

#Shannon Statistics
for (i in 1:length(heads)){
  sink("Shannon_Diversity_Statistics.txt", append=TRUE)
  head_clean <- as.name(heads[i])
  print(head_clean)
  physeq1.alpha.anova <- aov(as.formula(paste("Shannon~", head_clean)), data)
  log_1 <- summary(physeq1.alpha.anova )
  log_2<-TukeyHSD(physeq1.alpha.anova)
  print(log_1)
  print(log_2)
  cat("=============================\n")
  sink()
}


#Beta Diversity Statistics
#Bray-Ward Heat Map
##Distance Calculation
d = distance(physeq1, method='bray')

#Plot Ward Map
jpeg("Bray-Ward_Heat_Map.jpg")
plot(hclust(d, method="ward.D"))
heatmap(as.matrix(d))
dev.off()

#Jaccard-Ward Heat Map
##Distance Calculation
d = distance(physeq1, method='jaccard')

#Plot Ward Map
jpeg("Jaccard-Ward_Heat_Map.jpg")
plot(hclust(d, method="ward.D"))
heatmap(as.matrix(d))
dev.off()

#Jaccard Ordination

for (i in 1:length(heads)){
  #Subset Data from non NA
  physeq1_notna <- physeq1 %>% subset_samples(!is.na(get_variable(physeq1, heads[i])))
  #Calculate Distance Matrix
  dist.jc <- distance(physeq1_notna, method = "jaccard")
  ord <- ordinate(physeq1_notna, method = "MDS", distance=dist.jc )
  head_clean <- as.name(heads[i])
  #Ordination
  (ordplot <- plot_ordination(physeq1_notna, ord, heads[i], color=heads[i])) + theme_bw() 
  #Plot NMDS
  nmds <- ordplot + 
    stat_ellipse(type = "t") + theme_bw() +
    ggtitle('Jaccard NMDS') + theme(plot.title = element_text(hjust = 0.5)) +
    theme(panel.grid.minor.y=element_blank(),panel.grid.minor.x=element_blank(),panel.grid.major.y=element_blank(),panel.grid.major.x=element_blank()) +
    geom_point(size = 3)
  #Save Plot
  name <- paste0(head_clean,"_jaccard_nmds.png")
  ggsave(name,width = 20, height = 10)
}

#Bray Ordination

for (i in 1:length(heads)){
  #Subset Data from non NA
  physeq1_notna <- physeq1 %>% subset_samples(!is.na(get_variable(physeq1, heads[i])))
  #Calculate Distance Matrix
  dist.jc <- distance(physeq1_notna, method = "bray")
  ord <- ordinate(physeq1_notna, method = "MDS", distance=dist.jc )
  head_clean <- as.name(heads[i])
  # Ordination
  (ordplot <- plot_ordination(physeq1_notna, ord, heads[i], color=heads[i])) + theme_bw() 
  #Plot NMDS 
  nmds <- ordplot + 
    stat_ellipse(type = "t") + theme_bw() +
    ggtitle('Bray NMDS') + theme(plot.title = element_text(hjust = 0.5)) +
    theme(panel.grid.minor.y=element_blank(),panel.grid.minor.x=element_blank(),panel.grid.major.y=element_blank(),panel.grid.major.x=element_blank()) +
    geom_point(size = 3)
  #Save Plot
  name <- paste0(head_clean,"_bray_nmds.png")
  ggsave(name,width = 20, height = 10)
}

# Labeled Ward Linkage Maps

#Color Generation
col.grand <- wes_palette(name = "GrandBudapest2", n=4, type="discrete")
col.gr <- wes_palette(name = "Zissou", n=5, type="discrete")
col2.gr <- wes_palette(name = "Darjeeling2", n=5, type="discrete")
col.gr_extra <- wes_palette("Darjeeling", n=5, type="discrete")
col.rocket <- wes_palette("BottleRocket", n=7, type="discrete")
col.rush <- wes_palette("Rushmore", n=5, type="discrete")
col.moon1 <- wes_palette("Moonrise1", n=4, type="discrete")
col.moon2 <- wes_palette("Moonrise2", n=4, type="discrete")
col.moon3 <- wes_palette("Moonrise3", n=5, type="discrete")
col.chev <- wes_palette("Chevalier", n=4, type="discrete")

palleteall <- c(col.grand,col.gr,col2.gr,col.gr_extra,col.rocket,col.rush,col.moon1,col.moon2,col.moon3,col.chev)

#Plot Ward Linkage Maps
for (i in 1:length(heads)){
  #Subset Data for non NA
  physeq1_notna <- physeq1 %>% subset_samples(!is.na(get_variable(physeq1, heads[i])))
  #Calculate Distance Matrix
  dist.jc <- distance(physeq1_notna, method = "jaccard")
  #Cluster Data
  clustering <- as.phylo(hclust(dist.jc, method = "ward.D2"))
  #Select Environmental Variables
  envtype <- get_variable(physeq1_notna,heads[i])
  head_clean <- as.name(heads[i])
  #Color
  palette <- palleteall
  #Color Tips
  tipColor <- col_factor(palette, levels = levels(envtype))(envtype)
  name <- paste0(head_clean,"_ward_linkage.jpg")
  #Plot and Save Figure
  jpeg(name)
  par(mar=c(0,0,2,0))
  plot(clustering,tip.color = tipColor, direction ="downwards",main = "Ward Linkage Map", cex=1.5)
  dev.off()
}

#Network Map
for (i in 1:length(heads)){
  #Subset Sample for non NA
  physeq1_notna <- physeq1 %>% subset_samples(!is.na(get_variable(physeq1, heads[i])))
  #Network Generation
  ig <- make_network(physeq1_notna, distance=function(x){vegan::vegdist(x, "jaccard")}, max.dist=0.9)

  head_clean <- as.name(heads[i])
  #Generate and Save Plot
  name <- paste0(head_clean,"_network_map.jpg")
  jpeg(name)
  par(mar=c(0,0,2,0))
  plot_network(ig, physeq1_notna, color=heads[i], line_weight=0.9, label=NULL)
  dev.off()
}


#Beta Diversity Statistics

# make a data frame from the sample_data

# Adonis test
# Homogeneity of dispersion test
vegan_otu <- function(physeq) {
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

#Adonis, Paired Adonis, Betadisp from Vegan
for (i in 1:length(heads)){
  #subset data for non NA
  physeq1_notna <- physeq1 %>% subset_samples(!is.na(get_variable(physeq1, heads[i])))
  #Sample Data in data frame
  sampledf <- data.frame(sample_data(physeq1_notna))
  #Make Vegan OTU Table 
  vanilla <- vegan_otu(physeq1_notna)
  #Sample Data for Vegan
  sdf <- as(sample_data(physeq1_notna), "data.frame")
  # Calculate Jaccard distance matrix
  mor_jac <- phyloseq::distance(physeq1_notna, method = "jaccard")
  #Open .txt to record data
  sink("jaccard_beta_diversity_Statistics.txt", append=TRUE)
  head_clean <- as.name(heads[i])
  print(head_clean)
  #Adonis Test
  log_1 <- adonis(as.formula(paste("mor_jac~", head_clean)), data = sampledf)
  #Pairwise Adonis
  log_2<- pairwise.adonis(vanilla,getElement(sdf, heads[i]),sim.method = 'jaccard')
  #Print
  print(log_1)
  print(log_2)
  #Betadispersion Test
  beta <- betadisper(mor_jac, getElement(sdf, heads[i]))
  log_3 <- permutest(beta)
  log_4 <- TukeyHSD(beta)
  #Print
  print(log_3)
  print(log_4)
  cat("=============================================================================================\n")
  #Close .txt
  sink()
}
