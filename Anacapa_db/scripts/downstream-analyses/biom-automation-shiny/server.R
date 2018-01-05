library(shiny)
library(ggplot2)
library(qiimer)
library(reshape2)
library(ggrepel)
library(vegan)
library(dplyr)
library(wesanderson)
library(cluster)
library(biomformat)
library(phyloseq)
library(gridBase)
library(broom)
library(gridExtra)
library(plotly)
source("https://raw.githubusercontent.com/mahendra-mariadassou/phyloseq-extended/master/graphical_methods.R")
source("https://raw.githubusercontent.com/mahendra-mariadassou/phyloseq-extended/master/tree_methods.R")
source("https://raw.githubusercontent.com/mahendra-mariadassou/phyloseq-extended/master/plot_merged_trees.R")
source("https://raw.githubusercontent.com/mahendra-mariadassou/phyloseq-extended/master/specificity_methods.R")
source("https://raw.githubusercontent.com/mahendra-mariadassou/phyloseq-extended/master/ternary_plot.R")
source("https://raw.githubusercontent.com/mahendra-mariadassou/phyloseq-extended/master/richness.R")
source("https://raw.githubusercontent.com/mahendra-mariadassou/phyloseq-extended/master/edgePCA.R")
source("https://raw.githubusercontent.com/mahendra-mariadassou/phyloseq-extended/master/copy_number_correction.R")
source("https://raw.githubusercontent.com/mahendra-mariadassou/phyloseq-extended/master/copy_number_correction.R")
source("https://raw.githubusercontent.com/mahendra-mariadassou/phyloseq-extended/master/prevalence.R")
source("https://raw.githubusercontent.com/mahendra-mariadassou/phyloseq-extended/master/compute_niche.R")

custom_rarefaction <- function(physeq, sample.size = 100000, replicates = 10, sampledata,...) {
  reps  <- replicate(replicates, rarefy_even_depth(physeq, sample.size = sample.size))
  
  dfs <- lapply(reps, function(x) as.data.frame(x@otu_table@.Data))
  
  dfs <- lapply(dfs, function(x) rownames_to_column(x, var = "taxonomy"))
  dfs <- do.call(rbind.data.frame, dfs)
  
  otu <- dfs %>% group_by(taxonomy) %>% summarize_all(funs(mean)) %>% 
    mutate_if(is.numeric, funs(round)) %>% data.frame %>%
    column_to_rownames("taxonomy") %>% as.matrix
  
  OTU <- otu_table(otu, taxa_are_rows = T)
  
  TAX <- physeq@tax_table
  physeq_to_return <- phyloseq(OTU, TAX)
  physeq_to_return <- merge_phyloseq(physeq_to_return, physeq@sam_data)
  
  return(physeq_to_return)
}

options(digits = 5)
vegan_otu <- function(physeq) {
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

server <- function(input, output)({
  
  # Set up the input biom table ---------------
  dat <- reactive({
    req(input$in_biom)
    read_biom(input$in_biom$datapath)
  })
  
  otu <- reactive({
    req(input$in_biom)
    as.matrix(biom_data(dat()))
    otu_table(as.matrix(biom_data(dat())), taxa_are_rows = TRUE)
  })
  
  OTU_p <- reactive({
    req(input$in_biom)
    otu_table(otu(), taxa_are_rows = TRUE)
  })
  
  tax_phy <- reactive({
    req(input$in_biom)
    tp <- as.matrix(colsplit(rownames(otu()),";", names = c("Domain","Kingdom","Phylum","Class","Order","Family","Genus","Species")))
    rownames(tp) <- rownames(otu())
    tp
  })  
  
  TAX_p <- reactive({
    req(input$in_biom)
    tp <- tax_table(tax_phy())
    colnames(tp) <- c("Domain","Kingdom","Phylum","Class","Order","Family","Genus","Species")
    tp
  })
  

  # BiomTable to print
  output$print_biom <- renderDataTable({
    # req(input$in_biom)
    req(input$in_metadata)
    OTU_p()
  })
  

  # Set up the biom table metadata -----------------
  sampledata <- reactive({
    req(input$in_metadata)
    md <- read_qiime_mapping_file(input$in_metadata$datapath)
    md <- as.data.frame(md)
    rownames(md) <- md[,1]
    #order samples alphabetically
    md[order(md$SampleID),]
    sample_data(md)
  })
  
  heads <- reactive({
    req(input$in_metadata)
    colnames(sampledata())#[-1:-5]
  })
  
  # UI inputs ------
  output$which_variable_r <- renderUI({
    req(input$in_metadata)
    selectInput("var", "Select the variable", choices = heads())
  })
  output$which_divtype <- renderUI({
    req(input$in_metadata)
    radioButtons("divtype", label = "Observed or Shannon diversity?", 
                 choices = c("Observed", "Shannon"))
  })
  output$which_dissim <- renderUI({
    req(input$in_metadata)
    req(input$in_biom)
    radioButtons("dissimMethod", "Which type of distance metric would you like?", 
                 choices = c("bray", "jaccard"))
  })
  output$which_taxon_level <- renderUI({
    req(input$in_metadata)
    req(input$in_biom)
    radioButtons("taxon_level", "Pick the taxonomic level for making the plot", 
                 choices = c("Kingdom", "Phylum", "Class", "Order"))
    
  })
  output$rare_depth <- renderUI({
    sliderInput("rarefaction_depth", label = "Select a depth of rarefaction", min = 10000, max = 100000, step = 1000, value = 100000)
  })
  output$rare_reps <- renderUI({
    sliderInput("rarefaction_reps", label = "Select the number of times to rarefy", min = 2, max = 20, value = 2)
  })
  
  # Make physeq object ----------
  physeq <- reactive({
    req(input$in_biom)
    req(input$in_metadata)
    ps <- phyloseq(OTU_p(), TAX_p())
    merge_phyloseq(ps, sampledata())
  })
  
  
  output$ugly_bar <- renderPlotly({
    req(input$in_biom)
    req(input$in_metadata)
    plot_bar(physeq(), fill = input$taxon_level)
    ggplotly()
  })
  # Time to make some graphs -----------
  
  data_subset_unrare <- reactive({
    req(input$in_biom)
    req(input$in_metadata)
    
    p2 <- physeq()
    sample_data(p2) <- physeq() %>% sample_data %>% subset(., !is.na(get(input$var)))
    p2
  })
  
  data_subset <- reactive({
    req(input$in_biom)
    req(input$in_metadata)
    
    data_subset_unrare <- data_subset_unrare()
    custom_rarefaction(data_subset_unrare, sample.size = input$rarefaction_depth, replicates = input$rarefaction_reps)
  })
  
  output$testdf <- renderTable({
    req(input$in_biom)
    req(input$in_metadata)
    
    sample_data(data_subset()) %>% data.frame
  })
  
  output$classtest <- renderText({class})
  # Rarefaction curve after rarefaction -----------
  output$rarefaction_ur <- renderPlot({
    req(input$in_biom)
    req(input$in_metadata)
    p <- ggrare(data_subset_unrare(), step = 1000, se=FALSE, color = input$var)  
    q <- p + # facet_wrap(as.formula(paste("~", input$var))) + 
      theme_bw() + 
      theme(panel.grid.minor.y=element_blank(),panel.grid.minor.x=element_blank(),
            panel.grid.major.y=element_blank(),panel.grid.major.x=element_blank())
    q
  })

  output$rarefaction_r <- renderPlotly({
    req(input$in_biom)
    req(input$in_metadata)
    p <- ggrare(data_subset(), step = 1000, se=FALSE, color = input$var)  
    q <- p + facet_wrap(as.formula(paste("~", input$var))) + 
      theme_bw() + 
      theme(panel.grid.minor.y=element_blank(),panel.grid.minor.x=element_blank(),
            panel.grid.major.y=element_blank(),panel.grid.major.x=element_blank())
    ggplotly(tooltip = c("Sample", input$var))
  })
  
  # Alpha diverstity boxplots ----------
  output$alpharichness <- renderPlotly({
    req(input$in_biom)
    req(input$in_metadata)
    p <- plot_richness(data_subset(), x = input$var,  measures= input$divtype)
    q <- p + geom_boxplot(aes_string(fill = input$var, alpha=0.2, show.legend = F)) + theme_bw() + 
      xlab(paste(input$divtype, "Diversity")) +
      theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank(),
            panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank())
    ggplotly(tooltip = c("x", "value"))
  
  })
  
  # Alpha diversity aov generation
  physeq.alpha.anova <- reactive({
    req(input$in_biom)
    req(input$in_metadata)
    alpha.diversity <- estimate_richness(physeq(), measures = c("Observed", "Shannon"))
    data <- cbind(sample_data(physeq()), alpha.diversity)
    
    aov(as.formula(paste(input$divtype, "~" , input$var)), data)
  })
  
  # Alpha diversity AOV print --------
  output$alphaDivAOV <- renderTable({
    req(input$in_biom)
    req(input$in_metadata)
    broom::tidy(physeq.alpha.anova())
  }, digits = 4)
  
  # Alpha Diversity tukey --------
  output$alphaDivTukey <- renderTable({
    req(input$in_biom)
    req(input$in_metadata)
    broom::tidy(TukeyHSD(physeq.alpha.anova()))
  }, digits = 4)

  # NMDS plotly-----------
  output$betanmdsplotly <- renderPlotly({
    
    req(input$in_biom)
    req(input$in_metadata)
    d <- distance(data_subset(), method=input$dissimMethod)

    ord <- ordinate(data_subset(), method = "MDS", distance = d)

    nmdsplot <- plot_ordination(data_subset(), ord, input$var, color = input$var) +
      theme_bw() +  # stat_ellipse(type = "t", geom = "polygon", alpha = 0.2) +
      ggtitle(paste(input$var, "NMDS; dissimilarity method:",
                    tools::toTitleCase(input$dissimMethod))) + theme(plot.title = element_text(hjust = 0.5)) +
      theme(panel.grid.minor.y=element_blank(),panel.grid.minor.x=element_blank(),
            panel.grid.major.y=element_blank(),panel.grid.major.x=element_blank())# + theme_bw(base_size = 19)

    ggplotly(tooltip = c(input$var, "x", "y")) %>% layout(hovermode = 'closest') 
  })
  
  
  
  # Other beta diversity plots --------
  output$dissimMap <- renderPlot({
    
    req(input$in_biom)
    req(input$in_metadata)
    d <- distance(data_subset(), method=input$dissimMethod)

    # Network map ---------
    ig <- make_network(data_subset(), 
                       distance=function(x){vegan::vegdist(x, input$dissimMethod)}, max.dist=0.9)
    
    igp <-  plot_network(ig, data_subset(), color=input$var, line_weight=0.9, label=NULL) + 
      theme_bw(base_size = 18) 
    
    # Ward linkage map --------
    wcluster <- as.dendrogram(hclust(d, method = "ward.D2"))
    envtype <- get_variable(data_subset(), input$var)
    tipColor <- col_factor(rainbow(10), levels = levels(envtype))(envtype)
    wl <- ggdendro::ggdendrogram(wcluster, theme_dendro = FALSE, color = "red")  + 
      theme_bw(base_size = 18) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    # Big old heat map -----------
    heat <- ggplot(data = melt(as.matrix(d)), aes(x=Var1, y=Var2, fill=value)) +
      geom_tile() + theme_bw(base_size = 18) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))


    # Plot them all out ------
    gridExtra::grid.arrange(igp, wl, heat, ncol = 1, heights = c(3,3,3,4))
    
    
  }, height = 2000, width = 1000 )
  
  output$adonisTable <- renderTable ({
    sampledf <- data.frame(sample_data(data_subset()))
    mor_jac <- phyloseq::distance(data_subset(), method = input$dissimMethod)
    broom::tidy(adonis(as.formula(paste("mor_jac~", input$var)), data = sampledf)$aov.tab)
    
  }, digits = 5)
  
  output$permTestTable <- renderPrint({
    sdf <- as(sample_data(data_subset()), "data.frame")
    mor_jac <- phyloseq::distance(data_subset(), method = input$dissimMethod)
    (betadisper(mor_jac, getElement(sdf, input$var)))
  })
  
  output$betaTukey <- renderTable({
    sdf <- as(sample_data(data_subset()), "data.frame")
    mor_jac <- phyloseq::distance(data_subset(), method = input$dissimMethod)
    broom::tidy(TukeyHSD(betadisper(mor_jac, getElement(sdf, input$var))))
  }, digits = 5)
  
  output$heatmap <- renderD3heatmap ({
    req(input$in_biom)
    req(input$in_metadata)
    otu_table <- data_subset()@otu_table@.Data %>% data.frame %>% rownames_to_column("full_tax") %>% colsplit("full_tax", ";")
    
  }) 
  
  
})  
  