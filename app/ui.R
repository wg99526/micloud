library(seqinr)
library(shiny)
library(shinydashboard)
library(dashboardthemes)
library(tidyverse)
library(phyloseq)
library(plotly)
library(shinyWidgets)
library(shinyjs)
library(googleVis)
library(xtable)
library(DT)
library(htmltools)
library(biomformat)
library(phangorn)
library(bios2mds)
library(zip)

library(zCompositions) 
library(dplyr) 
library(forestplot) 
library(quantreg) 
library(fossil) 
library(picante) 
library(entropart) 

library(lme4) 
library(lmerTest) 
library(broom.mixed) 
library(gee) 
library(geepack) 

library(dirmult) 
library(robustbase) 
library(robCompositions) 
library(BiasedUrn) 
library(CompQuadForm) 
library(GUniFrac)  
library(ecodist)  
library(MiRKAT) 
library(GLMMMiRKAT)
library(devtools)
library(gridExtra) 
library(ggplot2) 
library(patchwork) 
library(ggthemes) 
library(erer) 
library(DiagrammeR) 
library(stringr) 
library(devtools)
library(betareg) 
library(reticulate) 
library(NBZIMM)
library(nlme) 
library(glmmTMB) 
library(glmm) 

{
  TITLE = p("MiCloud: Statistical Analysis of Microbiome Data", style = "font-size:18pt")
  HOME_COMMENT = p("MiCloud is a web application for statistical analysis of microbiome data. 
                   MiCloud provides step-by-step web environments for a breadth of data processing, analytic and graphical procedures. 
                   MiCloud can conduct both ecological (alpha- and beta-diversity) and taxonomical (phylum, class, order, family, genus, species) analyses for various types of host phenotypes (or disease status) and study designs with or without covariate adjustment(s). More details are as follows.", style = "font-size:13pt")
  HOME_COMMENT1 = ("Interactive procedures for various data inputs (.rdata, .rds, .biom, .txt, .csv, .tsv, .tre), quality controls (kingdom, library size, mean proportion) and data transformations (alpha- and beta-diversity, rarefaction, proportion, centered log-ratio).")
  HOME_COMMENT2 = ("Comparative/association analysis for both ecological (alpha- and beta-diversity) indices and taxonomical (phylum, class, order, family, genus, species) measurements.")
  HOME_COMMENT3 = ("Comparative/association analysis for both binary and continuous traits of host phenotypes (or medical interventions, disease status or environmental/behavioral factors).")
  HOME_COMMENT4 = ("Comparative/association analysis with or without covariate (e.g., age, gender) adjustment(s) for either cross-sectional or longitudinal/family-based microbiome studies.")
  HOME_COMMENT5 = ("Customizable/publishable data, graphs and tables.")
  INPUT_PHYLOSEQ_COMMENT1 = p("Description:", br(), br(), "This should be an '.Rdata' or '.rds' file, and the data should be in 'phyloseq' format (see ", 
                              a(tags$u("https://bioconductor.org/packages/release/bioc/html/phyloseq.html"), style = "color:red3"),
                              "). The phyloseq object should contain all the four necessary data, OTU/feature table, taxonomic table,
                              phylogenetic tree and sample data.", br(), br(), "Details:", br(), br(), 
                              "1) The OTU/feature table should contain counts, where rows are OTUs and columns are subjects 
                              (row names are OTU/feature IDs and column names are subject IDs).", br(),
                              "2) The taxonomic table should contain taxonomic names, where rows are OTUs and columns are seven taxonomic ranks 
                              (row names are OTU/feature IDs and column names are 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species' or 'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species').", br(),
                              "3) The phylogenetic tree should be a rooted tree. Otherwise, MiCloud automatically roots the tree through midpoint rooting (phangorn::midpoint). 
                              The tip labels of the phylogenetic tree are OTU/feature IDs.", br(), 
                              "4) The sample data should contain variables for the subjects about host phenotypes, medical interventions, disease status or environmental/behavioral factors, 
                              where rows are subjects and columns are variables (row names are subject IDs, and column names are variable names).", br(), br(), 
                              "* The OTUs should be matched and identical across OTU/feature table, taxonomic table and phylogenetic tree. 
                              The subjects should be matched and identical between OTU/feature table and sample data. MiCloud will analyze only the matched OTUs and subjects."
                              , style = "font-size:11pt")
  INPUT_PHYLOSEQ_COMMENT2 = p("You can download example microbiome data 'biom.Rdata' in 'phyloseq' format. The name of the 
                            phyloseq object should be 'biom'. For more details about 'phyloseq', see ", 
                              a(tags$u("https://bioconductor.org/packages/release/bioc/html/phyloseq.html"), style = "color:red3"),
                              br(), br(), "> setwd('/yourdatadirectory/')", br(), br(), "> load(file = 'biom.Rdata')", br(), br(), 
                              "> library(phyloseq)", br(), br(), " > otu.tab <- otu_table(biom)", br(), " > tax.tab <- 
                            tax_table(biom)", br(), " > tree <- phy_tree(biom)", br(), " > sam.dat <- sample_data(biom)", br(),
                              br(), "You can check if the OTUs are matched and identical across OTU/feature table, taxonomic table and
                            phylogenetic tree, and the subjects are matched and identical between OTU/feature table and sample data 
                            using following code.", br(), br(), " > identical(rownames(otu.tab), rownames(tax.tab))", br(), 
                              " > identical(rownames(otu.tab), tree$tip.label)", br(), " > identical(colnames(otu.tab), 
                            rownames(sam.dat))", style = "font-size:11pt")
  INPUT_INDIVIDUAL_DATA_COMMENT = p("Description:", br(), br(), 
                                    "1) The OTU/feature table (.txt or .csv) should contain counts, where rows are OTUs and columns are subjects (row names
                                  are OTU/feature IDs and column names are subject IDs). Alternatively, you can upload .biom file processed by QIIME", br(), 
                                    "2) The taxonomic table (.txt) should contain taxonomic names, where rows are OTUs and columns are seven taxonomic ranks (row names
                                  are OTU/feature IDs and column names are 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus',
                                  'Species' or 'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'). Alternatively, you can upload .tsv file processed by QIIME", br(), 
                                    "3) The phylogenetic tree (.tre or .nwk) should be a rooted tree. Otherwise, MiCloud
                                  automatically roots the tree through midpoint rooting (phangorn::midpoint). The tip labels of the
                                  phylogenetic tree are OTU/feature IDs.", br(), 
                                    "4) The sample data (.txt or .csv) should contain variables for the subjects about host phenotypes, medical interventions, disease status or environmental
                                  /behavioral factors, where rows are subjects and columns are variables (row names are subject IDs, and
                                  column names are variable names).", br(), br(), "* The OTUs should be matched and identical across OTU
                                  table, taxonomic table and phylogenetic tree. The subjects should be matched and identical between OTU
                                  table and sample data. MiCloud will analyze only the matched OTUs and subjects.", 
                                    style = "font-size:11pt")
}


{
  ui = dashboardPage(title="MiCloud",
                     #skin = "blue",
                     dashboardHeader(title=span(TITLE, 
                                                style = "font-size: 20px"), titleWidth = 550),
                     dashboardSidebar(
                       tags$script(JS("document.getElementsByClassName('sidebar-toggle')[0].style.visibility = 'hidden';")),
                       sidebarMenu(id = "side_menu",
                                   menuItem("Home", tabName = "home", icon = icon("home")),
                                   menuItem("Data Processing",  icon = icon("file-text-o"),
                                            menuSubItem("Data Input", tabName = "step1", icon = icon("mouse")),
                                            menuSubItem("Quality Control", tabName = "step2", icon = icon("chart-bar"))),
                                   menuItem("Ecological Analysis",  icon = icon("chart-pie"),
                                            menuSubItem("Diversity Calculation", tabName = "divCalculation", icon = icon("calculator")),
                                            menuSubItem("Alpha Diversity", tabName = "alphaDivanalysis", icon = icon("font")),
                                            menuSubItem("Beta Diversity", tabName = "betaDivanalysis", icon = icon("bold"))),
                                   menuItem("Taxonomical Analysis",  icon = icon("disease"),
                                            menuSubItem("Data Transformation", tabName = "dataTransform", icon = icon("th-large")),
                                            menuSubItem("Comparison/Association", tabName = "taxaAnalysis", icon = icon("align-left"))
                                   ))
                     ),
                     dashboardBody(
                       tags$head(tags$style(HTML(".content { padding-top: 2px;}"))),
                       tags$script(src="fileInput_text.js"),
                       useShinyjs(),
                       shinyDashboardThemes(theme = "onenote"),
                       uiOutput("themes"),
                       tabItems(
                         ##### HOME ####
                         tabItem(tabName = "home",
                                 #uiOutput("showhomepage"),
                                 div(id = "homepage", br(),
                                     strong("About", style = "font-size:17pt"), br(), HOME_COMMENT,
                                     tags$ol(
                                       tags$li(HOME_COMMENT1),
                                       tags$li(HOME_COMMENT2),
                                       tags$li(HOME_COMMENT3),
                                       tags$li(HOME_COMMENT4),
                                       tags$li(HOME_COMMENT5),
                                       style = "font-size:13pt"
                                     ),
                                     fluidRow(
                                       column(3,
                                              selectInput("selectTheme", strong("Select Theme", style = "font-size:14pt"), 
                                                          c("Choose one" = "", 
                                                            "Blue Gradient", "Flat Red", "Gray Dark", "Gray Light",
                                                            "Onenote (Default)", "Poor Mans Flatly", "Purple Gradient"
                                                          )))))
                         ),
                         ##### DATA INPUT ####
                         tabItem(tabName = "step1", br(),
                                 # div(strong("Step 1: Data Input", style = "font-size:17pt")),
                                 column(width = 6, style='padding-left:0px',
                                        box(
                                          width = NULL, status = "primary", solidHeader = TRUE,
                                          title = strong("Data Input", style = "color:black"),
                                          selectInput("inputOption", h4(strong("Data Type?")), 
                                                      c("Choose one" = "", "Phyloseq", "Individual Data"),
                                                      width = '30%'),
                                          div(id = "optionsInfo",
                                              tags$p("You can choose phyloseq or individual data.", 
                                                     style = "font-size:11pt"),
                                              style = "margin-top: -15px"),
                                          uiOutput("moreOptions")
                                        )
                                 ),
                                 column(width = 6, style='padding-left:0px',
                                        uiOutput("addDownloadinfo")
                                 )
                                 
                         ),
                         #),
                         ##### QC ####
                         tabItem(tabName = "step2", br(), 
                                 # strong("Step 2: Quality Control", style = "font-size:17pt"),
                                 sidebarLayout(
                                   position = "left",
                                   sidebarPanel(
                                     width = 3, 
                                     textInput("kingdom", h4(strong("Bacteria?")), value = "Bacteria"),
                                     p("A microbial kingdom to be analyzed. Default is 'Bacteria' for 16S data.
                                      Alternatively, you can type 'Fungi' for ITS data or any other kingdom of interest for
                                        shotgun metagenomic data.", style = "font-size:11pt"),
                                     tags$style(type = 'text/css', '#slider1 .irs-grid-text {font-size: 1px}'),
                                     tags$style(type = 'text/css', '#slider2 .irs-grid-text {font-size: 1px}'), 
                                     
                                     sliderInput("slider1", h4(strong("Library Size?")), min=0, max=10000, value = 3000, step = 1000),
                                     
                                     p("Remove subjects that have low library sizes (total read counts). Default is 3,000.", style = "font-size:11pt"),
                                     p("Library size: The total read count per subject.", style = "font-size:11pt"),
                                     sliderInput("slider2", h4(strong("Mean Proportion?")), min = 0, max = 0.1,
                                                 value = 0.002, step = 0.001,  post  = " %"),
                                     
                                     p("Remove OTUs/features that have low mean relative abundances (Unit: %). Default is 0.002%.",style = "font-size:11pt"),
                                     p("Mean proportion: The average of relative abundances (i.e., proportions) per OTU/feature.", style = "font-size:11pt"),             
                                     actionButton("run", (strong("Run!")), class = "btn-info"), 
                                     # actionButton("skip", (strong("Skip QC, Rarefy only")), class = "btn-info",
                                     #              style="background-color: #E1E0DD; border-color: #E1E0DD"), 
                                     br(), br(),
                                     uiOutput("moreControls")
                                   ),
                                   mainPanel(
                                     width = 9,
                                     fluidRow(width = 12,
                                              status = "primary", solidHeader = TRUE, 
                                              valueBoxOutput("sample_Size", width = 3),
                                              valueBoxOutput("OTUs_Size", width = 3),
                                              valueBoxOutput("phyla", width = 3),
                                              valueBoxOutput("classes", width = 3),
                                     ),
                                     fluidRow(width = 12, 
                                              status = "primary", solidHeader = TRUE,
                                              valueBoxOutput("orders", width = 3),
                                              valueBoxOutput("families", width = 3),
                                              valueBoxOutput("genera", width = 3),
                                              valueBoxOutput("species", width = 3)
                                     ),
                                     fluidRow(style = "position:relative",
                                              tabBox(width = 6,
                                                     title = strong("Library Size", style = "color:black"), 
                                                     tabPanel("Histogram",
                                                              plotlyOutput("hist"),
                                                              sliderInput("binwidth", "# of Bins:",min = 0, max = 100, value = 50, width = "100%"),
                                                              chooseSliderSkin("Round", color = "#112446")
                                                     ),
                                                     tabPanel("Box Plot",
                                                              plotlyOutput("boxplot")
                                                     )
                                              ),
                                              tabBox(width = 6,
                                                     title = strong("Mean Proportion", style = "color:black"), 
                                                     tabPanel("Histogram",
                                                              plotlyOutput("hist2"),
                                                              sliderInput("binwidth2", "# of Bins:",min = 0, max = 100, value = 50, 
                                                                          width = "100%"),
                                                              chooseSliderSkin("Round", color = "#112446")
                                                     ),
                                                     tabPanel("Box Plot",
                                                              plotlyOutput("boxplot2")))))
                                 )
                         ),
                         ##### DIVERSITY Calculation ####
                         tabItem(tabName = "divCalculation", br(),
                                 # div(strong("Diversity Calculation", style = "font-size:17pt")),
                                 column(
                                   width = 6, style='padding-left:0px',
                                   box(title = strong("Diversity Calculation", style = "color:black"), 
                                       width = NULL, status = "primary", solidHeader = TRUE,
                                       p("Calculate alpha-diversity indices: Richness (Observed), Shannon (Shannon, 1948), Simpson (Simpson, 1949), Inverse Simpson (Simpson, 1949), 
                                       Fisher (Fisher et al., 1943), Chao1 (Chao, 1984), ACE (Chao and Lee, 1992), ICE (Lee and Chao, 1994), PD (Faith, 1992)."), 
                                       p("Calculate beta-diversity indices: Jaccard dissimilarity (Jaccard, 1912), Bray-Curtis dissimilarity (Bray and Curtis, 1957), 
                                       Unweighted UniFrac distance (Lozupone and Knight, 2005), Generalized UniFrac distance (Chen et al., 2012), Weighted UniFrac distance (Lozupone et al., 2007)."), 
                                       actionButton("divCalcRun", (strong("Run!")), class = "btn-info"), 
                                   ),
                                   uiOutput("divCalcDownload")),
                                 column(
                                   width = 6, style='padding-left:0px',
                                   box(title = strong("References", style = "color:black"), 
                                       width = NULL, status = "primary", solidHeader = TRUE,
                                       p("Alpha Diversity", style = "font-size:12pt"),
                                       p("1. Chao A, Lee S. Estimating the number of classes via sample coverage. J Am Stat Assoc. 1992:87:210-217."),
                                       p("2. Chao A. Non-parametric estimation of the number of classes in a population. Scand J Stat. 1984:11:265-270."),
                                       p("3. Faith DP. Conservation evaluation and phylogenetic diversity. Biol Conserv. 1992:61:1-10."),
                                       p("4. Fisher RA, Corbet AS, Williams CB. The relation between the number of species and the number of individuals in a random sample of an animal population. J Anim Ecol. 1943:12:42-58"),
                                       p("5. Lee S, Chao A. Estimating population size via sample coverage for closed capture-recapture models. Biometrics. 1994:50:1:88-97."),      
                                       p("6. Shannon CE. A mathematical theory of communication. Bell Syst Tech J. 1948:27:379-423 & 623-656."),
                                       p("7. Simpson EH. Measurement of diversity. Nature 1949:163:688."), br(),
                                       p("Beta Diversity", style = "font-size:12pt"),
                                       p("1. Bray JR, Curtis JT. An ordination of the upland forest communities of Southern Wisconsin. Ecol Monogr. 1957;27(32549)."),
                                       p("2. Chen J, Bittinger K, Charlson ES, Hoffmann C, Lewis J, Wu GD., et al. Associating microbiome composition with environmental covariates using generalized UniFrac distances. Bioinformatics. 2012;28(16):2106-13."),
                                       p("3. Jaccard P. The distribution of the flora in the alpine zone. New Phytol. 1912;11(2):37-50."),
                                       p("4. Lozupone CA, Hamady M, Kelley ST, Knight R. Quantitative and qualitative 棺-diversity measures lead to different insights into factors that structure microbial communities. Appl Environ Microbiol. 2007;73(5):1576-85."),
                                       p("5. Lozupone CA, Knight R. UniFrac: A new phylogenetic method for comparing microbial communities. Appl Environ Microbiol. 2005;71(12):8228-35.")
                                   ))
                         ),
                         ##### ALPHA DIVERSITY ####
                         tabItem(tabName = "alphaDivanalysis", br(),
                                 # strong("Alpha Diversity", style = "font-size:17pt"),
                                 fluidRow(
                                   tabBox(width = 12,
                                          tabPanel(
                                            title ="Cross-Sectional",
                                            sidebarLayout(
                                              position = "left",
                                              sidebarPanel(width = 3,
                                                           uiOutput("primvars"),
                                                           uiOutput("prim_vars_types"),
                                                           uiOutput("covariates"), 
                                                           br(), 
                                                           uiOutput("alpha_downloadTable"),
                                                           uiOutput("alpha_references")
                                              ),
                                              mainPanel(width = 9,
                                                        fluidRow(width = 12, 
                                                                 uiOutput("alpha_display_results"))))
                                          ),
                                          tabPanel(
                                            title ="Longitudinal", 
                                            sidebarLayout(
                                              position = "left",
                                              sidebarPanel(
                                                width = 3, 
                                                uiOutput("primvars_long"),
                                                uiOutput("prim_vars_types_long"),
                                                uiOutput("covariates_long"), 
                                                br(), 
                                                uiOutput("alpha_downloadTablelong"),
                                                uiOutput("alpha_references_long")
                                              ),
                                              mainPanel(width = 9,
                                                        fluidRow(width = 12, 
                                                                 uiOutput("alpha_display_resultslong")))))))
                         ),
                         ##### BETA DIVERSITY ####
                         tabItem(tabName = "betaDivanalysis", br(),
                                 # strong("Beta Diversity", style = "font-size:17pt"),
                                 fluidRow(
                                   tabBox(width = 12,
                                          tabPanel(title ="Cross-Sectional",
                                                   sidebarLayout(
                                                     position = "left",
                                                     sidebarPanel(
                                                       width = 3,
                                                       #uiOutput("beta_expandCSoptions_croos"),
                                                       uiOutput("beta_primvar_cross"),
                                                       uiOutput("beta_prim_vars_types_cross"),
                                                       uiOutput("beta_covariates_cross"), 
                                                       br(), 
                                                       uiOutput("beta_downloadTable"),
                                                       uiOutput("beta_references")
                                                     ),
                                                     mainPanel(width = 9,
                                                               fluidRow(width = 12, 
                                                                        uiOutput("beta_display_results_cross")
                                                               )
                                                     )
                                                   )
                                          ),
                                          tabPanel(title ="Longitudinal",
                                                   sidebarLayout(
                                                     position = "left",
                                                     sidebarPanel(
                                                       width = 3,
                                                       #uiOutput("beta_expandCSoptionslong"),
                                                       uiOutput("beta_primvars_long"),
                                                       uiOutput("beta_prim_vars_types_long"),
                                                       uiOutput("beta_covariates_long"), 
                                                       br(), 
                                                       uiOutput("beta_downloadTablelong"),
                                                       uiOutput("beta_references_long")
                                                     ),
                                                     mainPanel(width = 9,
                                                               fluidRow(width = 12, 
                                                                        uiOutput("beta_display_resultslong")
                                                               )
                                                     )
                                                   )
                                          )
                                   )
                                 )
                         ),
                         ##### Data Transformation ####
                         tabItem(tabName = "dataTransform", br(),
                                 # div(strong("Data Transformation", style = "font-size:17pt")),
                                 column(
                                   width = 6, style='padding-left:0px',
                                   box(title = strong("Data Transformation", style = "color:black"), 
                                       width = NULL, status = "primary", solidHeader = TRUE,
                                       p("Transform the data into four different formats 1) count, 2) count (rarefied), 3) proportion,
                                         4) CLR (centered log ratio) (Aitchison, 1982) for each taxonomic rank (phylum, class, order, familiy, genus, species)."),
                                       actionButton("datTransRun", (strong("Run!")), class = "btn-info") 
                                   ),
                                   uiOutput("datTransDownload")),
                                 column(
                                   width = 6, style='padding-left:0px',
                                   box(title = strong("References", style = "color:black"), 
                                       width = NULL, status = "primary", solidHeader = TRUE,
                                       p("1. Aitchison J. The statistical analysis of compositional data. J R Statist Soc B. 1982;44:2:139-77")
                                   ))
                         ),
                         ##### Taxa Analysis ####
                         tabItem(tabName = "taxaAnalysis", br(),
                                 #h2(strong("Comparison / Association", style = "font-family: 'Verdana'; font-size:20pt")),
                                 fluidRow(
                                   tabBox(width = 12,
                                          tabPanel(
                                            title = "Cross-Sectional",
                                            sidebarLayout( 
                                              position = "left",
                                              sidebarPanel(width = 3,
                                                           uiOutput("primvars_taxa"),
                                                           uiOutput("morePrimvar_opt_taxa"),
                                                           uiOutput("covariates_taxa"),
                                                           br(),
                                                           uiOutput("downloadTable_taxa"),
                                                           uiOutput("taxa_references")
                                                           
                                              ),
                                              mainPanel(width = 9,
                                                        fluidRow(width = 12, 
                                                                 div(style='height:800px;overflow-y: scroll;', uiOutput("taxa_display")),
                                                                 br(),br(),
                                                                 uiOutput("taxa_display_dend")
                                                        )
                                              )
                                            )
                                          ),
                                          tabPanel(
                                            title = "Longitudinal",
                                            sidebarLayout( 
                                              position = "left",
                                              sidebarPanel(width = 3,
                                                           uiOutput("primvars_taxa.long"),
                                                           uiOutput("morePrimvar_opt_taxa.long"),
                                                           uiOutput("covariates_taxa.long"),
                                                           br(),
                                                           uiOutput("downloadTable_taxalong"),
                                                           uiOutput("taxa_references_long")
                                                           
                                              ),
                                              mainPanel(width = 9,
                                                        fluidRow(width = 12,
                                                                 div(style='height:800px;overflow-y: scroll;', uiOutput("taxa_displaylong")),
                                                                 br(),br(),
                                                                 uiOutput("taxa_displaylongd")
                                                        )))))))
                         
                       )))
}
