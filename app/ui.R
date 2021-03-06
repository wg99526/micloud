library(devtools)
library(tidyverse)
library(phangorn)
library(gridGraphics)

library(seqinr)
library(shiny)
library(shinydashboard)
library(dashboardthemes)

library(phyloseq)
library(plotly)
library(shinyWidgets)
library(shinyjs)
library(googleVis)
library(xtable)
library(DT)
library(htmltools)
library(biomformat)

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

library(gridExtra) 
library(ggplot2) 
library(patchwork) 
library(ggthemes) 
library(erer) 
library(DiagrammeR) 
library(stringr) 
library(betareg) 
library(reticulate) 
library(NBZIMM)
library(nlme) 
library(glmmTMB) 
library(glmm) 

{
  TITLE = p("MiCloud: A Unified Web Platform for Comprehensive Microbiome Data Analysis", style = "font-size:18pt")
  HOME_COMMENT = p(strong("MiCloud", style = "font-size:15pt"), "is a unified web platform for comprehensive microbiome data analysis. 
                   MiCloud provides step-by-step user-friendly web environments for a breadth of data processing, analytic and graphical procedures. 
                   MiCloud performs both ecological (alpha- and beta-diversity) and taxonomical (phylum, class, order, family, genus, 
                   species) analyses for various types of host phenotypes (or disease status) and study designs with or without covariate 
                   adjustment(s). More details are as follows.", style = "font-size:13pt")
  HOME_COMMENT1 = ("Interactive procedures for various data inputs (.rdata, .rds, .biom, .txt, .csv, .tsv, .tre), quality controls (kingdom, 
                   library size, mean proportion, taxonomic name) and data transformations (alpha- and beta-diversity, rarefaction, 
                   proportion, centered log-ratio).")
  HOME_COMMENT2 = ("Comparative analysis for both ecological (alpha- and beta-diversity) indices and microbial taxa in relative abundance 
                   (phylum, class, order, family, genus, species).")
  HOME_COMMENT3 = ("Comparative analysis for both binary and continuous traits of host phenotypes 
                   (or medical interventions, disease status or environmental/behavioral factors).")
  HOME_COMMENT4 = ("Comparative analysis with or without covariate (e.g., age, gender) adjustment(s) for either 
                   cross-sectional or longitudinal/family-based microbiome study design.")
  HOME_COMMENT5 = ("Adjustable/downloadable/publishable data, tables and graphs.")
  HOME_COMMENT6 = p("Reference: Gu, W., Moon, J., Chisina, C., Kang, B., Park, T., Koh, H. MiCloud: A unified web platform for comprehensive microbiome data analysis. (Under review)", style = "font-size:13pt")
  
  INPUT_PHYLOSEQ_COMMENT1 = p("Description:", br(), br(), "This should be an '.Rdata' or '.rds' file, and the data should be in 'phyloseq' format (see ", 
                              a(tags$u("https://bioconductor.org/packages/release/bioc/html/phyloseq.html"), style = "color:red3"),
                              "). The phyloseq object should contain all the four necessary data, feature (OTU or ASV) table, taxonomic table, 
                              metadata/sample information, and phylogenetic tree.", 
                              br(), br(), "Details:", br(), br(), 
                              "1) The feature table should contain counts, where rows are features (OTUs or ASVs) and columns are subjects 
                              (row names are feature IDs and column names are subject IDs).", br(),
                              "2) The taxonomic table should contain taxonomic names, where rows are features and columns are seven taxonomic ranks 
                              (row names are feature IDs and column names are 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species' or 
                              'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species').", br(),
                              "3) The metadata/sample information should contain variables for the subjects about host phenotypes, medical interventions, 
                              disease status or environmental/behavioral factors, where rows are subjects and columns are variables 
                              (row names are subject IDs, and column names are variable names).", br(), 
                              "4) The phylogenetic tree should be a rooted tree. Otherwise, MiCloud automatically roots the tree through midpoint rooting (phangorn::midpoint). 
                              The tip labels of the phylogenetic tree are feature IDs.", br(), br(), 
                              "* The features should be matched and identical across feature table, taxonomic table and phylogenetic tree. 
                              The subjects should be matched and identical between feature table and metadata/sample information. 
                              MiCloud will analyze only the matched features and subjects."
                              , style = "font-size:11pt")
  INPUT_PHYLOSEQ_COMMENT2 = p("You can download example microbiome data 'biom.Rdata' in 'phyloseq' format. The name of the 
                              phyloseq object should be 'biom'. For more details about 'phyloseq', see ", 
                              a(tags$u("https://bioconductor.org/packages/release/bioc/html/phyloseq.html"), style = "color:red3"), br(), br(), 
                              "> setwd('/yourdatadirectory/')", br(), br(), 
                              "> load(file = 'biom.Rdata')", br(), br(), 
                              "> library(phyloseq)", br(), br(), 
                              " > otu.tab <- otu_table(biom)", br(), 
                              " > tax.tab <- tax_table(biom)", br(), 
                              " > tree <- phy_tree(biom)", br(), 
                              " > sam.dat <- sample_data(biom)", br(), br(), 
                              "You can check if the features are matched and identical across feature table, taxonomic table and 
                              phylogenetic tree, and the subjects are matched and identical between feature table and metadata/sample information 
                              using following code.", br(), br(), 
                              " > identical(rownames(otu.tab), rownames(tax.tab))", br(), 
                              " > identical(rownames(otu.tab), tree$tip.label)", br(), 
                              " > identical(colnames(otu.tab), rownames(sam.dat))", style = "font-size:11pt")
  INPUT_INDIVIDUAL_DATA_COMMENT = p("Description:", br(), br(), 
                                    "1) The feature table (.txt or .csv) should contain counts, where rows are features (OTUs or ASVs) and columns are subjects 
                                    (row names are feature IDs and column names are subject IDs). Alternatively, you can upload .biom file processed by QIIME.", br(), 
                                    "2) The taxonomic table (.txt) should contain taxonomic names, where rows are features and columns are seven taxonomic ranks 
                                    (row names are feature IDs and column names are 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species' or 
                                    'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'). Alternatively, you can upload .tsv file processed by QIIME.", br(), 
                                    "3) The metadata/sample information (.txt or .csv) should contain variables for the subjects about host phenotypes, medical interventions, 
                                    disease status or environmental/behavioral factors, where rows are subjects and columns are variables (row names are subject IDs, and 
                                    column names are variable names).", br(), 
                                    "4) The phylogenetic tree (.tre or .nwk) should be a rooted tree. Otherwise, MiCloud automatically roots the tree through midpoint 
                                    rooting (phangorn::midpoint). The tip labels of the phylogenetic tree are feature IDs.", br(), br(), 
                                    "* The features should be matched and identical across feature table, taxonomic table and phylogenetic tree. 
                                    The subjects should be matched and identical between feature table and metadata/sample information. 
                                    MiCloud will analyze only the matched features and subjects.", style = "font-size:11pt")
  INPUT_INDIVIDUAL_DATA_COMMENT2 = p("You can download example microbiome data 'biom.zip'. This zip file contains four necessary data, feature table (otu.tab.txt), 
                                     taxonomic table (tax.tab.txt), metadata/sample information (sam.dat.txt), and phylogenetic tree (tree.tre).", br(), br(),
                                     "> setwd('/yourdatadirectory/')", br(), br(), 
                                     "> otu.tab <- read.table(file = 'otu.tab.txt', check.names = FALSE) ", br(), 
                                     "> tax.tab <- read.table(file = 'tax.tab.txt', check.names = FALSE)", br(), 
                                     " > sam.dat <- read.table(file = 'sam.dat.txt', check.names = FALSE) ", br(),
                                     "> tree <- read.tree(file = 'tree.tre')", br(), br(), 
                                     "You can check if the features are matched and identical across feature table, taxonomic table and phylogenetic tree, 
                                     and the subjects are matched and identical between feature table and metadata/sample information using following code.", br(), br(), 
                                     " > identical(rownames(otu.tab), rownames(tax.tab))", br(), 
                                     " > identical(rownames(otu.tab), tree$tip.label)", br(), 
                                     " > identical(colnames(otu.tab), rownames(sam.dat))", style = "font-size:11pt")
  
  EXTERNAL_RESOURCE_COMMENT = p("MiCloud does not take raw sequence data. For the raw sequence data processing and microbiome profiling, we recommend following popular and well-established bioinformatic pipelines.", br(), br(),
                              "For web platforms:", p(" ", style = "margin-bottom: 10px;"),"Nephele (https://nephele.niaid.nih.gov), Qiita (https://qiita.ucsd.edu), QIIME2 (q2studio) (https://qiime2.org) and PUMAA (https://sites.google.com/g.ucla.edu/pumaa)", br(), br(),
                              "For command line interfaces:", p(" ", style = "margin-bottom: 10px;"), "QIIME (http://qiime.org), QIIME2 (q2cli) (https://qiime2.org), MG-RAST (https://www.mg-rast.org), Mothur (https://mothur.org), MEGAN (http://ab.inf.uni-tuebingen.de/software/megan6) and MetaPhlAn (https://huttenhower.sph.harvard.edu/metaphlan)", style = "font-size:11pt")
  
  QC_KINGDOM_COMMENT = p("A microbial kingdom to be analyzed. Default is 'Bacteria' for 16S data. Alternatively, you can type 'Fungi' for ITS data 
                         or any other kingdom of interest for shotgun metagenomic data.", style = "font-size:11pt")
  QC_LIBRARY_SIZE_COMMENT1 = p("Remove subjects that have low library sizes (total read counts). Default is 3,000.", style = "font-size:11pt")
  QC_LIBRARY_SIZE_COMMENT2 = p("Library size: The total read count per subject.", style = "font-size:11pt")
  QC_MEAN_PROP_COMMENT1 = p("Remove features (OTUs or ASVs) that have low mean relative abundances (Unit: %). Default is 0.002%.",style = "font-size:11pt")
  QC_MEAN_PROP_COMMENT2 = p("Mean proportion: The average of relative abundances (i.e., proportions) per feature.", style = "font-size:11pt")
  QC_TAXA_NAME_COMMENT1 = p('Remove taxonomic names in the taxonomic table that are completely matched with the specified character strings. 
                            Multiple character strings should be separated by a comma. Default is "", "metagenome", "gut metagenome", "mouse gut metagenome".',
                            style = "font-size:11pt")
  QC_TAXA_NAME_COMMENT2 = p('Remove taxonomic names in the taxonomic table that are partially matched with the specified character strings (i.e., taxonomic names that contain 
                            the specified character strings). Multiple character strings should be separated by a comma. Default is "uncultured", "incertae", "Incertae",
                            "unidentified", "unclassified", "unknown".',style = "font-size:11pt")
  
  ALPHA_COMMENT = p("Calculate alpha-diversity indices: Richness (Observed), Shannon (Shannon, 1948), Simpson (Simpson, 1949), Inverse Simpson (Simpson, 1949), 
                    Fisher (Fisher et al., 1943), Chao1 (Chao, 1984), ACE (Chao and Lee, 1992), ICE (Lee and Chao, 1994), PD (Faith, 1992).")
  ALPHA_REFERENCES = p("1. Chao A, Lee S. Estimating the number of classes via sample coverage. J Am Stat Assoc. 1992:87:210-217.", br(),
                       "2. Chao A. Non-parametric estimation of the number of classes in a population. Scand J Stat. 1984:11:265-270.", br(),
                       "3. Faith DP. Conservation evaluation and phylogenetic diversity. Biol Conserv. 1992:61:1-10.", br(),
                       "4. Fisher RA, Corbet AS, Williams CB. The relation between the number of species and the number of individuals 
                       in a random sample of an animal population. J Anim Ecol. 1943:12:42-58.", br(),
                       "5. Lee S, Chao A. Estimating population size via sample coverage for closed capture-recapture models. Biometrics. 1994:50:1:88-97.", br(),
                       "6. Shannon CE. A mathematical theory of communication. Bell Syst Tech J. 1948:27:379-423 & 623-656.", br(),
                       "7. Simpson EH. Measurement of diversity. Nature 1949:163:688.", br())
  BETA_COMMENT = p("Calculate beta-diversity indices: Jaccard dissimilarity (Jaccard, 1912), Bray-Curtis dissimilarity (Bray and Curtis, 1957), Unweighted UniFrac distance 
                   (Lozupone and Knight, 2005), Generalized UniFrac distance (Chen et al., 2012), Weighted UniFrac distance (Lozupone et al., 2007).")
  BETA_REFERENCES = p("1. Bray JR, Curtis JT. An ordination of the upland forest communities of Southern Wisconsin. Ecol Monogr. 1957;27(32549).", br(),
                      "2. Chen J, Bittinger K, Charlson ES, Hoffmann C, Lewis J, Wu GD., et al. Associating microbiome composition with environmental 
                      covariates using generalized UniFrac distances. Bioinformatics. 2012;28(16):2106-13.", br(),
                      "3. Jaccard P. The distribution of the flora in the alpine zone. New Phytol. 1912;11(2):37-50.", br(),
                      "4. Lozupone CA, Hamady M, Kelley ST, Knight R. Quantitative and qualitative 棺-diversity measures lead to 
                      different insights into factors that structure microbial communities. Appl Environ Microbiol. 2007;73(5):1576-85.", br(),
                      "5. Lozupone CA, Knight R. UniFrac: A new phylogenetic method for comparing microbial communities. Appl Environ Microbiol. 2005;71(12):8228-35.")
  DATA_TRANSFORM_COMMENT = p("Transform the data into four different formats 1) count, 2) count (rarefied), 3) proportion, 
                             4) CLR (centered log ratio) (Aitchison, 1982) for each taxonomic rank (phylum, class, order, familiy, genus, species).")
  DATA_TRANSFORM_REFERENCE = p("1. Aitchison J. The statistical analysis of compositional data. J R Statist Soc B. 1982;44:2:139-77")
}

# UI
{
  ui = dashboardPage(
    title = "MiCloud",
    dashboardHeader(title = span(TITLE, style = "float:left;font-size: 20px"), titleWidth = "100%"),
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
                           menuSubItem("Comparison/Association", tabName = "taxaAnalysis", icon = icon("align-left"))))),
    dashboardBody(
      tags$head(tags$style(HTML(".content { padding-top: 2px;}"))),
      tags$script(src = "fileInput_text.js"),
      useShinyjs(),
      shinyDashboardThemes(theme = "onenote"),
      uiOutput("themes"),
      tabItems(
        
        ##### HOME ####
        tabItem(tabName = "home",
                div(id = "homepage", br(), HOME_COMMENT,
                    tags$ol(
                      tags$li(HOME_COMMENT1), tags$li(HOME_COMMENT2), tags$li(HOME_COMMENT3), tags$li(HOME_COMMENT4), tags$li(HOME_COMMENT5),
                      style = "font-size:13pt"),
                    HOME_COMMENT6, 
                    fluidRow(
                      column(3,
                             selectInput("selectTheme", strong("Select Theme", style = "font-size:14pt"), 
                                         c("Choose one" = "", "Flat Red", "Gray Dark", "Gray Light",
                                           "Onenote (Default)", "Poor Mans Flatly", "Purple Gradient")))))),
        
        ##### DATA INPUT ####
        tabItem(tabName = "step1", br(),
                column(width = 6, style='padding-left:0px',
                       box(
                         width = NULL, status = "primary", solidHeader = TRUE,
                         title = strong("Data Input", style = "color:black"),
                         selectInput("inputOption", h4(strong("Data Type?")), c("Choose one" = "", "Phyloseq", "Individual Data"), width = '30%'),
                         div(id = "optionsInfo", tags$p("You can choose phyloseq or individual data.", style = "font-size:11pt"), style = "margin-top: -15px"),
                         uiOutput("moreOptions"))),
                column(width = 6, style='padding-left:0px', uiOutput("addDownloadinfo"))),
        
        ##### QC ####
        tabItem(tabName = "step2", br(), 
                sidebarLayout(
                  position = "left",
                  sidebarPanel(width = 3,
                               textInput("kingdom", h4(strong("Bacteria?")), value = "Bacteria"),
                               QC_KINGDOM_COMMENT,
                               tags$style(type = 'text/css', '#slider1 .irs-grid-text {font-size: 1px}'),
                               tags$style(type = 'text/css', '#slider2 .irs-grid-text {font-size: 1px}'), 
                               
                               sliderInput("slider1", h4(strong("Library Size?")), min=0, max=10000, value = 3000, step = 1000),
                               QC_LIBRARY_SIZE_COMMENT1,
                               QC_LIBRARY_SIZE_COMMENT2,
                               
                               sliderInput("slider2", h4(strong("Mean Proportion?")), min = 0, max = 0.1, value = 0.002, step = 0.001,  post  = " %"),
                               QC_MEAN_PROP_COMMENT1,
                               QC_MEAN_PROP_COMMENT2,
                               
                               br(),
                               p(" ", style = "margin-bottom: -20px;"),
                               
                               h4(strong("Erroneous Taxonomic Names?")),
                               textInput("rem.str", label = "Complete Match", value = ""),
                               QC_TAXA_NAME_COMMENT1,
                               
                               textInput("part.rem.str", label = "Partial Match", value = ""),
                               QC_TAXA_NAME_COMMENT2,
                               
                               actionButton("run", (strong("Run!")), class = "btn-info"), br(), br(),
                               uiOutput("moreControls")),
                  mainPanel(width = 9,
                            fluidRow(width = 12,
                                     status = "primary", solidHeader = TRUE, 
                                     valueBoxOutput("sample_Size", width = 3),
                                     valueBoxOutput("OTUs_Size", width = 3),
                                     valueBoxOutput("phyla", width = 3),
                                     valueBoxOutput("classes", width = 3)),
                            fluidRow(width = 12, 
                                     status = "primary", solidHeader = TRUE,
                                     valueBoxOutput("orders", width = 3),
                                     valueBoxOutput("families", width = 3),
                                     valueBoxOutput("genera", width = 3),
                                     valueBoxOutput("species", width = 3)),
                            fluidRow(style = "position:relative",
                                     tabBox(width = 6, title = strong("Library Size", style = "color:black"), 
                                            tabPanel("Histogram",
                                                     plotlyOutput("hist"),
                                                     sliderInput("binwidth", "# of Bins:",min = 0, max = 100, value = 50, width = "100%"),
                                                     chooseSliderSkin("Round", color = "#112446")),
                                            tabPanel("Box Plot", 
                                                     plotlyOutput("boxplot"))),
                                     tabBox(width = 6, title = strong("Mean Proportion", style = "color:black"), 
                                            tabPanel("Histogram",
                                                     plotlyOutput("hist2"),
                                                     sliderInput("binwidth2", "# of Bins:",min = 0, max = 100, value = 50, width = "100%"),
                                                     chooseSliderSkin("Round", color = "#112446")),
                                            tabPanel("Box Plot",
                                                     plotlyOutput("boxplot2"))))))),
        
        ##### DIVERSITY Calculation ####
        tabItem(tabName = "divCalculation", br(),
                column(width = 6, style = 'padding-left:0px',
                       box(title = strong("Diversity Calculation", style = "color:black"), width = NULL, status = "primary", solidHeader = TRUE,
                           ALPHA_COMMENT, 
                           BETA_COMMENT, 
                           actionButton("divCalcRun", (strong("Run!")), class = "btn-info")),
                       uiOutput("divCalcDownload")),
                column(width = 6, style='padding-left:0px',
                       box(title = strong("References", style = "color:black"), width = NULL, status = "primary", solidHeader = TRUE,
                           p("Alpha Diversity", style = "font-size:12pt"),
                           ALPHA_REFERENCES,
                           p("Beta Diversity", style = "font-size:12pt"),
                           BETA_REFERENCES))),
        
        ##### ALPHA DIVERSITY ####
        tabItem(tabName = "alphaDivanalysis", br(),
                fluidRow(
                  tabBox(width = 12,
                         tabPanel(
                           title ="Cross-Sectional",
                           sidebarLayout(
                             position = "left",
                             sidebarPanel(width = 3,
                                          uiOutput("primvars"),
                                          uiOutput("prim_vars_types"),
                                          uiOutput("covariates"), br(), 
                                          uiOutput("alpha_downloadTable"),
                                          uiOutput("alpha_references")),
                             mainPanel(width = 9,
                                       fluidRow(width = 12, 
                                                uiOutput("alpha_display_results"))))),
                         tabPanel(
                           title ="Longitudinal", 
                           sidebarLayout(
                             position = "left",
                             sidebarPanel(width = 3, 
                                          uiOutput("primvars_long"),
                                          uiOutput("prim_vars_types_long"),
                                          uiOutput("covariates_long"), br(), 
                                          uiOutput("alpha_downloadTablelong"),
                                          uiOutput("alpha_references_long")),
                             mainPanel(width = 9,
                                       fluidRow(width = 12, 
                                                uiOutput("alpha_display_resultslong")))))))),
        
        ##### BETA DIVERSITY ####
        tabItem(tabName = "betaDivanalysis", br(),
                fluidRow(
                  tabBox(width = 12,
                         tabPanel(
                           title ="Cross-Sectional",
                           sidebarLayout(
                             position = "left",
                             sidebarPanel(width = 3,
                                          uiOutput("beta_primvar_cross"),
                                          uiOutput("beta_prim_vars_types_cross"),
                                          uiOutput("beta_covariates_cross"), br(), 
                                          uiOutput("beta_downloadTable"),
                                          uiOutput("beta_references")),
                             mainPanel(width = 9,
                                       fluidRow(width = 12, 
                                                uiOutput("beta_display_results_cross"))))),
                         tabPanel(
                           title ="Longitudinal",
                           sidebarLayout(
                             position = "left",
                             sidebarPanel(width = 3,
                                          uiOutput("beta_primvars_long"),
                                          uiOutput("beta_prim_vars_types_long"),
                                          uiOutput("beta_covariates_long"), br(), 
                                          uiOutput("beta_downloadTablelong"),
                                          uiOutput("beta_references_long")),
                             mainPanel(width = 9,
                                       fluidRow(width = 12, 
                                                uiOutput("beta_display_resultslong")))))))),
        
        ##### Data Transformation ####
        tabItem(tabName = "dataTransform", br(),
                column(width = 6, style='padding-left:0px',
                       box(title = strong("Data Transformation", style = "color:black"), width = NULL, status = "primary", solidHeader = TRUE,
                           DATA_TRANSFORM_COMMENT,
                           actionButton("datTransRun", (strong("Run!")), class = "btn-info") ),
                       uiOutput("datTransDownload")),
                column(width = 6, style='padding-left:0px', 
                       box(title = strong("References", style = "color:black"), width = NULL, status = "primary", solidHeader = TRUE,
                           DATA_TRANSFORM_REFERENCE))),
        
        ##### Taxa Analysis ####
        tabItem(tabName = "taxaAnalysis", br(),
                fluidRow(
                  tabBox(width = 12,
                         tabPanel(
                           title = "Cross-Sectional",
                           sidebarLayout( 
                             position = "left",
                             sidebarPanel(width = 3,
                                          uiOutput("primvars_taxa"),
                                          uiOutput("morePrimvar_opt_taxa"),
                                          uiOutput("covariates_taxa"), br(),
                                          uiOutput("downloadTable_taxa"),
                                          uiOutput("taxa_references")),
                             mainPanel(width = 9,
                                       fluidRow(width = 12, 
                                                div(style='height:800px;overflow-y: scroll;', uiOutput("taxa_display")), br(),br(),
                                                uiOutput("taxa_display_dend"))))),
                         tabPanel(
                           title = "Longitudinal",
                           sidebarLayout( 
                             position = "left",
                             sidebarPanel(width = 3,
                                          uiOutput("primvars_taxa.long"),
                                          uiOutput("morePrimvar_opt_taxa.long"),
                                          uiOutput("covariates_taxa.long"), br(),
                                          uiOutput("downloadTable_taxalong"),
                                          uiOutput("taxa_references_long")),
                             mainPanel(width = 9,
                                       fluidRow(width = 12,
                                                div(style='height:800px;overflow-y: scroll;', uiOutput("taxa_displaylong")), br(),br(),
                                                uiOutput("taxa_displaylongd"))))))))
      )
    )
  )
}
