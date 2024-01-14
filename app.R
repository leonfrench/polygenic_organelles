#TODO put on github
#TODO put on shinyapps
#TODO add python notebook to project
library(matrixTests)
library(matrixStats)
library(pheatmap)
library(readxl)
library(here)
library(shiny)
library(readr)
library(magrittr)
library(tibble)
library(tidyr)
library(dplyr)
library(purrr)
#library(shinyjs)
source("./AUCFunction.R")
source("./string_processing.R")
#deploy with:
#rsconnect::deployApp(appFileManifest = "appFileManifest.txt", account = "poly-brain") 

ui <- fluidPage(
  shinyjs::useShinyjs(),
  #tags$head(includeHTML("google-analytics.html")),
  # App title ----
  
  titlePanel("Polygenic tester for the global organelle proteome profiling dataset from the Chan Zuckerberg Biohub"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    # Sidebar panel for inputs ----
    sidebarPanel(
      textAreaInput(
        inputId = "genelist",
        label = "Input your gene list:",
        value = 'APOE\nABCA1\nABCA7\nABI3\nACE\nADAM17\nADAMTS1\nANK3\nANKH\nAPH1B\nAPP\nBCKDK\nBIN1\nBLNK\nCASS4\nCD2AP\nCLNK\nCLU\nCOX7C\nCR1\nCTSB\nCTSH\nDOC2A\nEED\nEPDR1\nEPHA1\nFERMT2\nFOXF1\nGRN\nHLA-DQA1\nHS3ST5\nICA1\nIDUA\nIGHG3\nIGHV3-65\nIL34\nINPP5D\nJAZF1\nKLF16\nLILRB2\nMAF\nMINDY2\nMME\nMS4A4A\nMYO15A\nNCK2\nPLCG2\nPLEKHA1\nPRDM7\nPRKD3\nPTK2B\nRASGEF1C\nRBCK1\nRHOH\nSCIMP\nSEC61G\nSHARPIN\nSIGLEC11\nSLC24A4\nSLC2A4RG\nSNX1\nSORL1\nSORT1\nSPDYE3\nSPI1\nSPPL2A\nTMEM106B\nTNIP1\nTPCN1\nTREM2\nTREML2\nTSPAN14\nTSPOAP1\nUMAD1\nUNC5CL\nUSP6NL\nWDR12\nWDR81\nWNT3',
        rows = 10
      ),
      textAreaInput(
        inputId = "background_genelist",
        label = "Background gene list (optional):",
        value = '',
        rows = 2
      ),
      selectInput(
        inputId = 'species',
        label = 'Species of input genes:',
        choices = c('Human', 'Mouse', 'Rhesus Macaque')
      ),
      selectInput(
        inputId = 'data_version',
        label = 'Analysis level:',
        choices = c('Marker', 'Organelle')
      ),
      actionButton(inputId = "submit",
                   label = "Submit"),
      br(),
      br(),
      downloadButton(outputId = "download_data", label = "Download results as .csv"),
      downloadButton(outputId = "download_heatmap", label = "Download heatmap"),
      hr(),
      tags$b("Data was made available by Hein, Peng, Todorova, McCarthy, Kim, Liu and collegues and is available at: "),
      br(),
      
      tags$a(href="https://organelles.czbiohub.org/", "Online web tool"),
      br(),
      tags$a(href="https://www.biorxiv.org/content/10.1101/2023.12.18.572249v1.full", "Preprint (2023)"),
      br(),
      br(),
      tags$div("Default gene list is from the Bellenguez et al. GWAS of Alzheimer's disease."),
      tags$a(href="https://www.nature.com/articles/s41588-022-01024-z", "(Nature genetics, 2022, Table S5)"),
      br()
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      div(
        id = "main",
        # Output: Verbatim text for data summary ----
        verbatimTextOutput("summary"),
        br(),
        dataTableOutput("view"),
        br(),
        #plotlyOutput("dotplot"),
        verbatimTextOutput("info")
      )
    )
  )
)


# Define server logic process and output top cortical layers/zones ----
server <- function(input, output) {
  output$summary <- renderPrint({
    cat("\nResults will load here when complete. The default gene list is Alzheimer's disease associated genes from Bellenguez et al. (Table S5).")
    cat("\n")
    #print(gc())
    #print(Sys.info()['nodename'])
  })
  
  observeEvent(input$submit, {
    start <- Sys.time()
    
    cleaned_gene_list <-
      isolate(process_input_genes(input$genelist))
    
    
    # load reference data
    if (input$data_version == "Marker") {
      cell_expression_ranks <- read_csv(here("data", "organelles.czbiohub", "processed_ranks_markers.csv"))
      cell_meta_data <- read_csv(here("data", "organelles.czbiohub", "processed_ranks_markers_metadata.csv"))
    } else {
      cell_expression_ranks <- read_csv(here("data", "organelles.czbiohub", "processed_ranks_organelle.csv"))
      cell_meta_data <- read_csv(here("data", "organelles.czbiohub", "processed_ranks_organelle_metadata.csv"))
    }
    cell_meta_data %<>% mutate(cell_type = as.character(cell_type))
    #gene_universe <- read_table('./data/gene_universe.txt', col_names = F) %>% .$X1
    #cell_expression_ranks %<>% filter(gene_symbol %in% gene_universe)
    cleaned_gene_list <- convert2human(input_genes = cleaned_gene_list, in_species = input$species)
    first_cleaned_gene_list <- cleaned_gene_list
    background_cleaned_gene_list <- isolate(process_input_genes(input$background_genelist))
    background_cleaned_gene_list <- convert2human(input_genes = background_cleaned_gene_list, in_species = input$species)
    
    if (length(background_cleaned_gene_list) == 1 && background_cleaned_gene_list == "") {
      background_cleaned_gene_list <- cell_expression_ranks$gene_symbol
    } else { #if given a background
      cell_expression_ranks %<>% filter(gene_symbol %in% background_cleaned_gene_list)
      #slow - could use matrix version, or skipped if no background loaded
      cell_expression_ranks %<>% mutate_if(is.numeric, rank)
    }
    cleaned_gene_list <- intersect(cleaned_gene_list, background_cleaned_gene_list)
    #re rank based on new background
    #for indices - use dplyr for ease
    forIndices <- as_tibble(cell_expression_ranks$gene_symbol)
    names(forIndices) <- 'gene_symbol'
    forIndices %<>% mutate(isTargetGene = gene_symbol %in% cleaned_gene_list)
    targetIndices <- forIndices$isTargetGene
    
    df <- cell_expression_ranks %>%
      select(-gene_symbol)
    
    AUROC <- auroc_analytic_ranked_profiles(df, targetIndices)
    
    get_p_from_AUC <- function(n.x, n.y, AUC, alternative = "two.sided", correct = T) {
      n.x <- as.double(n.x)
      n.y <- as.double(n.y)
      
      STATISTIC = AUC*(n.x*n.y)
      z <- STATISTIC - n.x * n.y/2
      SIGMA <- sqrt((n.x * n.y/12) * ((n.x + n.y + 1) ))
      
      if (correct) {
        CORRECTION <- switch(alternative, two.sided = sign(z) * 0.5, greater = 0.5, less = -0.5)
      } else {
        CORRECTION <- 0
      }
      z <- (z - CORRECTION)/SIGMA
      PVAL <- switch(alternative, less = pnorm(z), greater = pnorm(z, lower.tail = FALSE), two.sided = 2 * min(pnorm(z), pnorm(z, lower.tail = FALSE)))
      PVAL
    }
    
    #wilcox_tests <- col_wilcoxon_twosample(df[targetIndices, ], df[!targetIndices,])
    #wilcox_tests <- as_tibble(t(wilcox_tests[, "pvalue", drop=F]), )
    
    # group results together in a single table
    table <- gather(AUROC, key = cell_type, value = AUROC)
    
    table %<>% rowwise() %>% mutate(pValue = get_p_from_AUC(length(cleaned_gene_list),length(background_cleaned_gene_list), AUROC))
    table %<>% ungroup()
    
    # these are the values for the results table
    table %<>% mutate(pValue = signif(pValue, digits = 3), 
                      AUROC = signif(AUROC, digits = 3),
                      adjusted_P = signif(p.adjust(pValue, method = "bonferroni"), digits = 3))
    
    #add in meta data
    table <- inner_join(cell_meta_data, table)
    #table %<>% select(cell_type, )
    if (input$data_version == "Marker") {
      table %<>% select(Organelle, Marker, AUROC, pValue, adjusted_P, everything())
      table %<>% select(-cell_type)
    } else {
      table %<>% select(Organelle = cell_type, Markers, AUROC, pValue, adjusted_P, everything())
    }
    
    table %<>% arrange(-AUROC)
    
    output$summary <- renderPrint({
      #count of intersection of submitted genes with total gene list
      cat(paste("Time taken:", round(Sys.time() - start), "seconds"))
      cat(paste(
        "\nGenes found in data:",
        length(cleaned_gene_list),
        "of",
        length(first_cleaned_gene_list)
      ))
      cat(paste(
        "\nBackground genes:", length(background_cleaned_gene_list)
      ))
    })
    
    output$view <- renderDataTable(    table, escape = FALSE, options = list(
      paging =FALSE,
      pageLength =  51 
    ))
    
    #code duplication with monkey
    output$download_heatmap <-
      downloadHandler(
        filename = "polygenic_heatmap.pdf",
        content = function(file) {
          df <- as.data.frame(cell_expression_ranks %>% filter(gene_symbol %in% cleaned_gene_list))
          rownames(df) <- df$gene_symbol
          df <- df[ , (names(df) != "gene_symbol"), drop = T]
          
          plot_out <- pheatmap(df, main = paste0("Heatmap for ", length(cleaned_gene_list), " genes\ncolor scale represents specific expression rank (higher is more specific)"))
          
          if (input$data_version == "Class") {
            width = 5
          } else {
            width = 7
          }
          
          pdf(file, width=width, height=3+nrow(df) *.2)
          grid::grid.newpage()
          grid::grid.draw(plot_out$gtable)
          dev.off()
        }
      )
    
    
    
    output$download_data <-
      downloadHandler(
        filename = "polygenic_HPA_results.csv",
        content = function(file) {
          write_csv(table, file)
        }
      )
  })
}

shinyApp(ui, server)