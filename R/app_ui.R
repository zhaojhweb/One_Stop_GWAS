#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_ui <- function(request) {
  dashboardPage(
    skin = 'blue',
    dashboardHeader(title = 'rGWAS'),
    dashboardSidebar(
      sidebarMenu(
        menuItem('GWAS Structural', tabName = 'layout1', icon = icon('square-poll-vertical')),
        menuItem('GWAS Phenotypic', tabName = 'layout2', icon = icon('leaf')),
        menuItem('GWAS Run', tabName = 'layout3', icon = icon('dashboard')),
        menuItem('GWAS Visualization', tabName = 'layout4', icon = icon('images')),
        menuItem('GWAS SNP to Gene', tabName = 'layout5', icon = icon('table')),
        menuItem('GO/KEGG Enrichment', tabName = 'layout6', icon = icon('align-left')),
        menuItem('GO/KEGG Visualization', tabName = 'layout7', icon = icon('images'))
      )
    ),
    dashboardBody(
      tabItems(
        tabItem(
          tabName = 'layout1',
          fluidRow(
            box(title = "Upload hmp file",
                solidHeader = TRUE, status = "primary",
                fileInput("uploadhmp",
                          "Please upload your hmp file",
                          accept = c(".txt", ".hmp"),
                          buttonLabel = "Browse...",
                          placeholder = "No file selected",
                          capture = NULL),
                width = 12),
            box(title = "Generate a phylogenetic tree",
                solidHeader = TRUE, status = "primary",
                selectInput("Tree_type", "Please select a tree type",
                            choices = c("rectangular", "dendrogram","slanted","ellipse","roundrect","fan","circular","inward_circular","radial","equal_angle","daylight","ape")),
                actionButton("structural_analyse", 'Generate', icon = icon('play')),
                plotOutput("PhyloTree_plot"),
                numericInput("tree_width", "image_width", value = 10),
                numericInput("tree_height", "image_height", value = 10),
                radioButtons("tree_farmat", "Select the format in which to download the image", choices = c("jpeg", "pdf", "png", "svg", "tiff")),
                downloadButton("download_tree", "Download"),
                width = 12),
            box(title = "Principal component analysis",
                solidHeader = TRUE, status = "primary",
                width = 12,
                box(title = "2D type",
                    solidHeader = TRUE, status = "primary",
                    actionButton("PCA_2D", 'Generate', icon = icon('play')),
                    numericInput("X_2D", "X-axis num", value = '1', min = '1'),
                    numericInput("y_2D", "Y-axis num", value = '2', min = '1'),
                    plotOutput("PCA_plot_2D"),
                    numericInput("PCA_plot_2D_width", "image_width", value = 10),
                    numericInput("PCA_plot_2D_height", "image_height", value = 10),
                    radioButtons("PCA_plot_2D_farmat", "Select the format in which to download the image", choices = c("jpeg", "pdf", "png", "svg", "tiff")),
                    downloadButton("download_PCA_plot_2D", "Download"),
                    width = 6),
                box(title = "3D type",
                    solidHeader = TRUE, status = "primary",
                    actionButton("PCA_3D", 'Generate', icon = icon('play')),
                    numericInput("X_3D", "X-axis num", value = '1', min = '1'),
                    numericInput("y_3D", "Y-axis num", value = '2', min = '1'),
                    numericInput("Z_3D", "Z-axis num", value = '3', min = '1'),
                    plotOutput("PCA_plot_3D"),
                    numericInput("PCA_plot_3D_width", "image_width", value = 10),
                    numericInput("PCA_plot_3D_height", "image_height", value = 10),
                    radioButtons("PCA_plot_3D_farmat", "Select the format in which to download the image", choices = c("jpeg", "pdf", "png", "svg", "tiff")),
                    downloadButton("download_PCA_plot_3D", "Download"),
                    width = 6)
            )
          )
        ),
        tabItem(tabName = 'layout2',
                fluidRow(
                  box(title = "Calculate the BLUP",
                      solidHeader = TRUE, status = "primary",
                      fileInput("uploadpheno", "Please upload your phenotype file",
                                accept = c(".txt"), buttonLabel = "Browse...",
                                placeholder = "No file selected",
                                capture = TRUE),
                      actionButton("Calblup", "Calculation", icon = icon('edit')),
                      dataTableOutput("Bluptable"),
                      radioButtons("blup_table_format", "Select the format for downloading the form", choices = c("csv","tsv", "txt")),
                      downloadButton("download_blup", "Download"),
                      width = 12),
                  box(title = "HistPlot of mean",
                      solidHeader = TRUE, status = "primary",
                      plotOutput("Histplot_mean"),
                      numericInput("Histplot_width", "image_width", value = 10),
                      numericInput("Histplot_height", "image_height", value = 10),
                      radioButtons("Histplot_faormat", "Select the format in which to download the image", choices = c("jpeg", "pdf", "png", "svg", "tiff")),
                      downloadButton("download_histplot", "Download"),
                      width = 12),
                  box(title = "Phenotypic correlation analysis",
                      solidHeader = TRUE, status = "primary",
                      fileInput("Uploadtraits", "Please upload your processed phenotypic file",
                                accept = c(".txt"), buttonLabel = "Browse...",
                                placeholder = "No file selected",
                                capture = TRUE),
                      selectInput("Cor_method", "Please choice a method", choices = c("pearson", "kendall", "spearman")),
                      actionButton("Cal_corr", "Cal correlation", icon = icon('play')),
                      dataTableOutput("Cor_tab"),
                      radioButtons("Cor_format", "Select the format for downloading the form", choices = c("csv","tsv", "txt")),
                      downloadButton("download_Cor", "Download"),
                      plotOutput("Cor_heatmap"),
                      numericInput("Cor_heatmap_width", "image_width", value = 10),
                      numericInput("Cor_heatmap_height", "image_height", value = 10),
                      radioButtons("Cor_heatmap_farmat", "Select the format in which to download the image", choices = c("jpeg", "pdf", "png", "svg", "tiff")),
                      downloadButton("download_Cor_heatmap", "Download"),
                      width = 12
                  )
                )),
        tabItem(tabName = 'layout3',
                fluidRow(
                  box(title = "upload your data",
                      solidHeader = TRUE, status = "primary",
                      fileInput("uploadpheno_GWAS", "Please upload your new pheno file",
                                accept = c(".txt"), buttonLabel = "Browse...",
                                placeholder = "No file selected",
                                capture = TRUE),
                      fileInput("uploadgeno", "Please upload your new geno file",
                                accept = c(".txt"), buttonLabel = "Browse...",
                                placeholder = "No file selected",
                                capture = TRUE),
                      selectInput("GWAS_model", "GWAS model", choices = c("BLINK", "CMLM", "GLM", "MLM", "MMLM", "SUPER", "FarmCPU", "EMMAxP3D")),
                      numericInput("PCA_num", "PCA number", value = '3'),
                      actionButton("RunGWAS", "Run GWAS", icon = icon('random')),
                      width = 12),
                  box(title = "GWAS results",
                      solidHeader = TRUE, status = "primary",
                      dataTableOutput("GWASresult"),
                      radioButtons("GWAS_table_format", "Select the format for downloading the form", choices = c("csv","tsv", "txt")),
                      downloadButton("download_GWAS", "Download the result"),
                      width = 12)
                )),
        tabItem(tabName = 'layout4',
                fluidRow(
                  box(title = "QQ plot",
                      solidHeader = TRUE, status = "primary",
                      plotOutput("QQ_plot"),
                      actionButton("Plot_QQ", "Plot", icon = icon('play')),
                      numericInput("QQ_width", "image_width", value = 10),
                      numericInput("QQ_height", "image_height", value = 10),
                      radioButtons("QQ_farmat", "Select the format in which to download the image", choices = c("jpeg", "pdf", "png", "svg", "tiff")),
                      downloadButton("download_QQ", "Download"),
                      width = 12),
                  box(title = "Manhattan plot",
                      solidHeader = TRUE, status = "primary",
                      numericInput("Manhattan_y", "Select a threshold", value = '3'),
                      plotlyOutput("Manhattan_plot", width = "100%", height = "400px"),
                      actionButton("Plot_Manhattan", "Plot", icon = icon('play')),
                      numericInput("Manhattan_width", "image_width", value = 10),
                      numericInput("Manhattan_height", "image_height", value = 10),
                      radioButtons("Manhattan_farmat", "Select the format in which to download the image", choices = c("jpeg", "pdf", "png", "svg", "tiff")),
                      downloadButton("download_Manhattan", "Download"),
                      width = 12),
                  box(title = "LD heatmap",
                      solidHeader = TRUE, status = "primary",
                      textInput("LD_Chromosome", "Select the Chromosome for LD heatmap", value = "1A"),
                      numericInput("LD_left", "left_Pos", value = '1000000'),
                      numericInput("LD_right", "right_Pos", value = '5000000'),
                      plotOutput("LD_plot"),
                      actionButton("Plot_LD", "Plot", icon = icon('play')),
                      numericInput("LD_width", "image_width", value = 10),
                      numericInput("LD_height", "image_height", value = 10),
                      radioButtons("LD_farmat", "Select the format in which to download the image", choices = c("jpeg", "pdf", "png", "svg", "tiff")),
                      downloadButton("download_LD", "Download"),
                      width = 12)
                )),
        tabItem(tabName = 'layout5',
                fluidRow(
                  box(title = "SNP to Gene",
                      solidHeader = TRUE, status = "primary",
                      fileInput("uploadGXF", "Please upload your GXF file",
                                accept = c(".txt", ".gff", ".gtf", ".gff", ".gff3", ".GFF", ".GTF", ".GFF", ".GFF3"), buttonLabel = "Browse...",
                                placeholder = "No file selected",
                                capture = TRUE),
                      numericInput("up_down_stream", "Location range(bp)", value = '100000'),
                      numericInput("Significant_SNP", "The P value is selected according to the GWAS result", value = '0.01'),
                      actionButton("Extract", "Extract", icon = icon('play')),
                      width = 12),
                  box(title = "Result",
                      solidHeader = TRUE, status = "primary",
                      dataTableOutput("SNP_Gene"),
                      radioButtons("SNP_Gene_format", "Select the format for downloading the form", choices = c("csv","tsv", "txt")),
                      downloadButton("download_SNP_Gene", "Download"),
                      width = 12
                  )
                )),
        tabItem(tabName = 'layout6',
                fluidRow(
                  box(title = "Enrichment",
                      solidHeader = TRUE, status = "primary",
                      fileInput("GeneInfo", "Please upload your annatation file",
                                accept = c(".txt"), buttonLabel = "Browse...",
                                placeholder = "No file selected",
                                capture = TRUE),
                      fileInput("GeneList", "Please upload your Gene list",
                                accept = c(".txt"), buttonLabel = "Browse...",
                                placeholder = "No file selected",
                                capture = TRUE),
                      width = 12),
                  box(title = "Run GO Enrichment",
                      solidHeader = TRUE, status = "primary",
                      numericInput("GO_pvalue", "Select the P-value", value = '1'),
                      #numericInput("GO_qvalue", "select the Q-value", value = '1'),
                      actionButton("Run_GO", "Run", icon = icon('play')),
                      dataTableOutput("GO_result"),
                      radioButtons("GO_table_format", "Select the format for downloading the form", choices = c("csv","tsv", "txt")),
                      downloadButton("download_GO_table", "Download"),
                      width = 12),
                  box(title = "Run KEGG Enrichment",
                      solidHeader = TRUE, status = "primary",
                      numericInput("KEGG_pvalue", "Select the P-value", value = '1'),
                      #numericInput("KEGG_qvalue", "select the Q-value", value = '1'),
                      actionButton("Run_KEGG", "Run", icon = icon('play')),
                      dataTableOutput("KEGG_result"),
                      radioButtons("KEGG_table_format", "Select the format for downloading the form", choices = c("csv","tsv", "txt")),
                      downloadButton("download_KEGG_table", "Download"),
                      width = 12)
                )),
        tabItem(tabName = 'layout7',
                fluidRow(
                  box(title = "PLot GO Enrichment result",
                      solidHeader = TRUE, status = "primary",
                      numericInput("BP_num", "Shows the number of BP", value = '10'),
                      numericInput("CC_num", "Shows the number of CC", value = '10'),
                      numericInput("MF_num", "Shows the number of MF", value = '10'),
                      actionButton("Plot_GO", "Plot", icon = icon('play')),
                      plotlyOutput("GO_plot", width = "100%", height = "400px"),
                      numericInput("GO_width", "image_width", value = 10),
                      numericInput("GO_height", "image_height", value = 10),
                      radioButtons("GO_farmat", "Select the format in which to download the image", choices = c("jpeg", "pdf", "png", "svg", "tiff")),
                      downloadButton("download_GO", "Download"),
                      width = 12),
                  box(title = "Plot KEGG Enrichment result",
                      solidHeader = TRUE, status = "primary",
                      numericInput("KEGG_num", "Shows the number of Pathway", value = '10'),
                      actionButton("Plot_KEGG", "Plot", icon = icon('play')),
                      plotlyOutput("KEGG_plot", width = "100%", height = "400px"),
                      numericInput("KEGG_width", "image_width", value = 10),
                      numericInput("KEGG_height", "image_height", value = 10),
                      radioButtons("KEGG_farmat", "Select the format in which to download the image", choices = c("jpeg", "pdf", "png", "svg", "tiff")),
                      downloadButton("download_KEGG", "Download"),
                      width = 12)
                ))
      )
    )
  )
}

#' Add external Resources to the Application
#'
#' This function is internally used to add external
#' resources inside the Shiny application.
#'
#' @import shiny
#' @importFrom golem add_resource_path activate_js favicon bundle_resources
#' @noRd
golem_add_external_resources <- function() {
  add_resource_path(
    "www",
    app_sys("app/www")
  )

  tags$head(
    favicon(),
    bundle_resources(
      path = app_sys("app/www"),
      app_title = "rGWAS"
    )
    # Add here other external resources
    # for example, you can add shinyalert::useShinyalert()
  )
}
