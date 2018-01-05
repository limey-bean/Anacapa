library(plotly)
library(shiny)
shinyUI(pageWithSidebar(
  headerPanel("Exploring output from Anacapa pipeline"),
  sidebarPanel(
    
    ## conditionalPanel() functions for selected tab
    conditionalPanel(condition="input.tabselected==1",
                     fileInput("in_biom", "Select a BIOM table to upload"),
                     fileInput("in_metadata", "Select metadata")),

    conditionalPanel(condition="input.tabselected == 3 | input.tabselected == 4 | input.tabselected == 5 | input.tabselected == 6", uiOutput("which_variable_r")),
    conditionalPanel(condition="input.tabselected == 4", uiOutput("which_divtype")),
    conditionalPanel(condition="input.tabselected == 5 | input.tabselected == 6", uiOutput("which_dissim")),
    conditionalPanel(condition="input.tabselected == 7", uiOutput("which_taxon_level")),
    conditionalPanel(condition="input.tabselected == 3", uiOutput("rare_depth")),
    conditionalPanel(condition="input.tabselected == 3", uiOutput("rare_reps"))
    
  ),
  
  mainPanel(
    tabsetPanel(
      tabPanel("About", value=1, helpText("Select a biom table and a metadata file")),
      tabPanel("View OTU table", value=2, helpText("Here's the OTU table (taxon name not displayed)"), dataTableOutput("print_biom")),
      tabPanel("Rarefaction curve", value = 3, plotOutput("rarefaction_ur"),
               plotlyOutput("rarefaction_r")),
      tabPanel("Alpha Diversity exploration", value = 4, plotlyOutput("alpharichness"), 
               tableOutput("alphaDivAOV"), tableOutput("alphaDivTukey")),
      tabPanel("Beta Diversity exploration", value = 5, plotlyOutput("betanmdsplotly"), plotOutput("dissimMap")),
      tabPanel("Beta Diversity stats", value = 6, tableOutput("adonisTable"),
               verbatimTextOutput("permTestTable"),
               tableOutput("betaTukey")),
      tabPanel("ugly-bar", value = 7, plotlyOutput("ugly_bar")),
      
      
      id = "tabselected"
    )
  )
))