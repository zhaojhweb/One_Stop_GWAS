#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd

app_server <- function(input, output, session) {
  options(shiny.maxRequestSize=3000*1024^2)
  ######################################################
  ##################structural_analyse##################
  hmpdata <- reactive({
    req(input$uploadhmp)
    data <- read.delim(input$uploadhmp$datapath, header = F)
    return(data)
  })

  observeEvent(input$structural_analyse,{
    withProgress(message = "Calculating",
                 detail = "It will take some time ^-^",
                 expr = {
                   for (i in 1:5){
                     incProgress(1/5)
                     Sys.sleep(0.5)
                   }

                   hmp_num <- reactive(GAPIT(G = hmpdata(), output.numerical = F, file.output = F))

                   treedata <- hmp_num()$GD
                   rownames(treedata) <- treedata[,1]
                   treedata_df <- dplyr::select(treedata, -1)
                   tree_plot <- nj(dist.gene(treedata_df))
                   tree_type <- input$Tree_type
                   output$PhyloTree_plot <- renderPlot({
                     ggtree(tree_plot, layout = tree_type, branch.length = "none") +
                       geom_tiplab()
                   })})

    output$download_tree <- downloadHandler(
      filename = function(){
        paste0("Phylotree.",input$tree_farmat, sep = "")
      },
      content = function(file){
        withProgress(message = "Downloading",
                     detail = "It will take some time ^-^",
                     expr = {
                       for (i in 1:3){
                         incProgress(1/3)
                         Sys.sleep(0.5)
                       }
                       if (input$tree_farmat == "jpeg"){
                         jpeg(file, height = input$tree_height, width = input$tree_width, units = "in", res = 800)
                       } else if (input$tree_farmat == "pdf"){
                         pdf(file, height = input$tree_height, width = input$tree_width)
                       } else if (input$tree_farmat == "png"){
                         png(file, height = input$tree_height, width = input$tree_width, units = "in", res = 800)
                       } else if (input$tree_farmat == "svg"){
                         svg(file, height = input$tree_height, width = input$tree_width)
                       } else if (input$tree_farmat == "tiff"){
                         tiff(file, height = input$tree_height, width = input$tree_width, units = "in", res = 800)
                       }
                       p_tree <- ggtree(tree_plot, layout = tree_type, branch.length = "none") +
                         geom_tiplab()
                       plot(p_tree)
                       dev.off()
                     }
        )})
  })
  ######################################################
  ######################################################



  ######################################################
  ########################PCA_2D########################
  hmpdata <- reactive({
    req(input$uploadhmp)
    data <- read.delim(input$uploadhmp$datapath, header = F)
    return(data)
  })
  observeEvent(input$PCA_2D,{
    withProgress(message = "Drawing",
                 detail = "It will take some time ^-^",
                 expr = {
                   for (i in 1:5){
                     incProgress(1/5)
                     Sys.sleep(0.5)
                   }
                   hmp_num <- reactive(GAPIT(G = hmpdata(), output.numerical = F, file.output = F))
                   pcx_2D <- input$X_2D
                   pcy_2D <- input$y_2D
                   pcadata <- data.frame(hmp_num()$GD, row.names = 1 )
                   mypca_2D <- stats::prcomp(pcadata)
                   output$PCA_plot_2D <- renderPlot({
                     plot(col = 'red', x = mypca_2D$x[,pcx_2D], y = mypca_2D$x[,pcy_2D],
                          pch = 20, xlab = paste0("PC", pcx_2D),
                          ylab = paste0("PC", pcy_2D))
                   })})

    output$download_PCA_plot_2D <- downloadHandler(
      filename = function(){
        paste0("PCA_2D.",input$PCA_plot_2D_farmat, sep = "")
      },
      content = function(file){
        withProgress(message = "Downloading",
                     detail = "It will take some time ^-^",
                     expr = {
                       for (i in 1:3){
                         incProgress(1/3)
                         Sys.sleep(0.5)
                       }
                       if (input$PCA_plot_2D_farmat == "jpeg"){
                         jpeg(file, height = input$PCA_plot_2D_height, width = input$PCA_plot_2D_width, units = "in", res = 800)
                       } else if (input$PCA_plot_2D_farmat == "pdf"){
                         pdf(file, height = input$PCA_plot_2D_height, width = input$PCA_plot_2D_width)
                       } else if (input$PCA_plot_2D_farmat == "png"){
                         png(file, height = input$PCA_plot_2D_height, width = input$PCA_plot_2D_width, units = "in", res = 800)
                       } else if (input$PCA_plot_2D_farmat == "svg"){
                         svg(file, height = input$PCA_plot_2D_height, width = input$PCA_plot_2D_width)
                       } else if (input$PCA_plot_2D_farmat == "tiff"){
                         tiff(file, height = input$PCA_plot_2D_height, width = input$PCA_plot_2D_width, units = "in", res = 800)
                       }
                       plot(col = 'red', x = mypca_2D$x[,pcx_2D], y = mypca_2D$x[,pcy_2D],
                            pch = 20, xlab = paste0("PC", pcx_2D),
                            ylab = paste0("PC", pcy_2D))
                       dev.off()
                     }
        )}
    )
  })
  ######################################################
  ######################################################



  ######################################################
  ##########################PCA_3D######################
  observeEvent(input$PCA_3D,{
    withProgress(message = "Drawing",
                 detail = "It will take some time ^-^",
                 expr = {
                   for (i in 1:5){
                     incProgress(1/5)
                     Sys.sleep(0.5)
                   }
                   hmp_num <- reactive(GAPIT(G = hmpdata(), output.numerical = F, file.output = F))
                   pcx_3D <- input$X_3D
                   pcy_3D <- input$y_3D
                   pcz_3D <- input$Z_3D
                   pcadata <- data.frame(hmp_num()$GD, row.names = 1 )
                   mypca_3D <- stats::prcomp(pcadata)
                   output$PCA_plot_3D <- renderPlot({
                     scatterplot3d::scatterplot3d(x = mypca_3D$x[,pcx_3D],
                                                  y = mypca_3D$x[,pcy_3D],
                                                  z = mypca_3D$x[,pcz_3D],
                                                  xlab = paste0("PC", pcx_3D),
                                                  ylab = paste0("PC", pcy_3D),
                                                  zlab = paste0("PC", pcz_3D),
                                                  color= 'red',
                                                  pch = 20)
                   })})

    output$download_PCA_plot_3D <- downloadHandler(
      filename = function(){
        paste0("PCA_3D.",input$PCA_plot_3D_farmat, sep = "")
      },
      content = function(file){
        withProgress(message = "Downloading",
                     detail = "It will take some time ^-^",
                     expr = {
                       for (i in 1:3){
                         incProgress(1/3)
                         Sys.sleep(0.5)
                       }
                       if (input$PCA_plot_3D_farmat == "jpeg"){
                         jpeg(file, height = input$PCA_plot_3D_height, width = input$PCA_plot_3D_width, units = "in", res = 800)
                       } else if (input$PCA_plot_2D_farmat == "pdf"){
                         pdf(file, height = input$PCA_plot_3D_height, width = input$PCA_plot_3D_width)
                       } else if (input$PCA_plot_2D_farmat == "png"){
                         png(file, height = input$PCA_plot_3D_height, width = input$PCA_plot_3D_width, units = "in", res = 800)
                       } else if (input$PCA_plot_2D_farmat == "svg"){
                         svg(file, height = input$PCA_plot_3D_height, width = input$PCA_plot_3D_width)
                       } else if (input$PCA_plot_2D_farmat == "tiff"){
                         tiff(file, height = input$PCA_plot_3D_height, width = input$PCA_plot_3D_width, units = "in", res = 800)
                       }
                       scatterplot3d::scatterplot3d(x = mypca_3D$x[,pcx_3D],
                                                    y = mypca_3D$x[,pcy_3D],
                                                    z = mypca_3D$x[,pcz_3D],
                                                    xlab = paste0("PC", pcx_3D),
                                                    ylab = paste0("PC", pcy_3D),
                                                    zlab = paste0("PC", pcz_3D),
                                                    color= 'red',
                                                    pch = 20)
                       dev.off()
                     }
        )}
    )
  })
  ######################################################
  ######################################################



  ######################################################
  #########################Calblup######################
  phenodata <- reactive({
    req(input$uploadpheno)
    pheno <- read.table(input$uploadpheno$datapath, header = T, sep = "\t")
    return(pheno)
  })
  observeEvent(input$Calblup,{
    withProgress(message = "Calculating",
                 detail = "It will take some time ^-^",
                 expr = {
                   for (i in 1:3){
                     incProgress(1/3)
                     Sys.sleep(0.5)
                   }

                   blup_data <- reactive({
                     data <- blup(phenodata(), sample = "Line", loc = "Loc", rep = "Rep", year = "Year", phe = "Brix", fold = 1.5)
                   })

                   output$Bluptable <- renderDT({
                     datatable(blup_data(), options = list(scrollX = TRUE))
                   })
                   output$Histplot_mean <- renderPlot({
                     histplot(blup_data()$Mean)
                   })})

    output$download_blup <- downloadHandler(
      filename = function(){
        paste0("Pheno_BLUP.", input$blup_table_format, sep = "")
      },
      content = function(file){
        withProgress(message = "Downloading",
                     detail = "It will take some time ^-^",
                     expr = {
                       for (i in 1:3){
                         incProgress(1/3)
                         Sys.sleep(0.5)
                       }
                       if (input$blup_table_format == "csv"){
                         write.csv(blup_data(), file)
                       } else if (input$blup_table_format == "tsv"){
                         write_tsv(blup_data(), file)
                       } else if (input$blup_table_format == "txt"){
                         write.table(blup_data(), file)
                       }
                     }
        )})

    output$download_histplot <- downloadHandler(
      filename = function(){
        paste0("HistPlot_Mean.",input$Histplot_faormat, sep = "")
      },
      content = function(file){
        withProgress(message = "Downloading",
                     detail = "It will take some time ^-^",
                     expr = {
                       for (i in 1:3){
                         incProgress(1/3)
                         Sys.sleep(0.5)
                       }
                       if (input$Histplot_faormat == "jpeg"){
                         jpeg(file, height = input$Histplot_height, width = input$Histplot_width, units = "in", res = 800)
                       } else if (input$Histplot_faormat == "pdf"){
                         pdf(file, height = input$Histplot_height, width = input$Histplot_width)
                       } else if (input$Histplot_faormat == "png"){
                         png(file, height = input$Histplot_height, width = input$Histplot_width, units = "in", res = 800)
                       } else if (input$Histplot_faormat == "svg"){
                         svg(file, height = input$Histplot_height, width = input$Histplot_width)
                       } else if (input$Histplot_faormat == "tiff"){
                         tiff(file, height = input$Histplot_height, width = input$Histplot_width, units = "in", res = 800)
                       }
                       histplot(blup_data()$Mean)
                       dev.off()
                     }
        )})

  })
  ######################################################
  ######################################################



  ######################################################
  ########################Cal_corr######################
  traitdata <- reactive({
    req(input$Uploadtraits)
    traitdata <- read.table(input$Uploadtraits$datapath, header = T, row.names = 1)
    traitdata_clean <- na.omit(traitdata)
    return(traitdata_clean)
  })

  observeEvent(input$Cal_corr, {
    Cormethod <- input$Cor_method
    withProgress(message = "Calculating",
                 detail = "It will take some time ^-^",
                 expr = {
                   for (i in 1:3){
                     incProgress(1/3)
                     Sys.sleep(0.5)
                   }
                   trait_cor <- reactive(cor(traitdata(), method = Cormethod))
                   output$Cor_tab <- renderDT({
                     datatable(trait_cor())
                   })
                   output$Cor_heatmap <- renderPlot({
                     ggcorrplot(trait_cor(),
                                ggtheme = ggplot2::theme_void())
                   })})
    output$download_Cor_heatmap <- downloadHandler(
      filename = function(){
        paste0("Corr_.",input$Cor_method,".", input$Cor_heatmap_farmat, sep = "")
      },
      content = function(file){
        withProgress(message = "Downloading",
                     detail = "It will take some time ^-^",
                     expr = {
                       for (i in 1:3){
                         incProgress(1/3)
                         Sys.sleep(0.5)
                       }
                       if (input$Cor_heatmap_farmat == "jpeg"){
                         jpeg(file, height = input$Cor_heatmap_height, width = input$Cor_heatmap_width, units = "in", res = 800)
                       } else if (input$Cor_heatmap_farmat == "pdf"){
                         pdf(file, height = input$Cor_heatmap_height, width = input$Cor_heatmap_width)
                       } else if (input$Cor_heatmap_farmat == "png"){
                         png(file, height = input$Cor_heatmap_height, width = input$Cor_heatmap_width, units = "in", res = 800)
                       } else if (input$Cor_heatmap_farmat == "svg"){
                         svg(file, height = input$Cor_heatmap_height, width = input$Cor_heatmap_width)
                       } else if (input$Cor_heatmap_farmat == "tiff"){
                         tiff(file, height = input$Cor_heatmap_height, width = input$Cor_heatmap_width, units = "in", res = 800)
                       }
                       p_Cor_heatmap <- ggcorrplot(trait_cor(),
                                                   ggtheme = ggplot2::theme_void())
                       plot(p_Cor_heatmap)
                       dev.off()
                     }
        )}
    )
    output$download_Cor <- downloadHandler(
      filename = function(){
        paste0("Corr_", input$Cor_method,".", input$Cor_format, sep = "")
      },
      content = function(file){
        withProgress(message = "Downloading",
                     detail = "It will take some time ^-^",
                     expr = {
                       for (i in 1:3){
                         incProgress(1/3)
                         Sys.sleep(0.5)
                       }
                       if (input$Cor_format == "csv"){
                         write.csv(trait_cor(), file)
                       } else if (input$Cor_format == "tsv"){
                         write_tsv(trait_cor(), file)
                       } else if (input$Cor_format == "txt"){
                         write.table(trait_cor(), file)
                       }
                     }
        )}
    )
  })
  ######################################################
  ######################################################



  ######################################################
  ########################RunGWAS#######################
  GWAS_Geno <- reactive({
    req(input$uploadgeno)
    geno <- read.delim(input$uploadgeno$datapath, header = F)
    return(geno)
  })
  GWAS_pheno <-reactive({
    req(input$uploadpheno_GWAS)
    pheno_GWAS <- read.delim(input$uploadpheno_GWAS$datapath, header = T)
    return(pheno_GWAS)
  })
  observeEvent(input$RunGWAS, {

    withProgress(message = "Calculating",
                 detail = "It will take some time ^-^",
                 expr = {
                   for (i in 1:10){
                     incProgress(1/10)
                     Sys.sleep(0.5)
                   }
                   GWAS_Geno_num <- reactive(GAPIT(G = GWAS_Geno(), output.numerical = F, file.output = F))

                   PCA_num <- input$PCA_num

                   GWAS_model <- input$GWAS_model

                   GWAS_result <- GAPIT(Y = GWAS_pheno()[,c(1,2)], GD = GWAS_Geno_num()$GD, GM = GWAS_Geno_num()$GM, PCA.total = PCA_num, model = GWAS_model, file.output = F)

                   GWAS_reslut_final <- dplyr::select(GWAS_result$GWAS, !(2:3)) %>%
                     left_join(GWAS_Geno_num()$GM, by = c("SNP" = "SNP"))

                   output$GWASresult <- renderDT({
                     datatable(GWAS_reslut_final, options = list(scrollX = TRUE))
                   })})

    output$download_GWAS <- downloadHandler(
      filename = function(){
        paste0("GWAS_result.", input$download_GWAS, sep = "")
      },
      content = function(file){
        withProgress(message = "Downloading",
                     detail = "It will take some time ^-^",
                     expr = {
                       for (i in 1:3){
                         incProgress(1/3)
                         Sys.sleep(0.5)
                       }
                       if (input$GWAS_table_format == "csv"){
                         write.csv(GWAS_reslut_final, file)
                       } else if (input$GWAS_table_format == "tsv"){
                         write_tsv(GWAS_reslut_final, file)
                       } else if (input$GWAS_table_format == "txt"){
                         write.table(GWAS_reslut_final, file)
                       }
                     }
        )}
    )


    observeEvent(input$Plot_QQ, {
      withProgress(message = "Drawing",
                   detail = "It will take some time ^-^",
                   expr = {
                     for (i in 1:3){
                       incProgress(1/3)
                       Sys.sleep(0.5)
                     }

                     GWAS_reslut_final <- dplyr::select(GWAS_result$GWAS, !(2:3)) %>%
                       left_join(GWAS_Geno_num()$GM, by = c("SNP" = "SNP"))
                     output$QQ_plot <- renderPlot(

                       qqman::qq(GWAS_reslut_final$P.value,
                                 col = 'blue',
                                 cex = 2,
                                 las = 0.5,
                                 pch = 20)
                     )})

      output$download_QQ <- downloadHandler(
        filename = function(){
          paste0("QQ.",input$QQ_farmat, sep = "")
        },
        content = function(file){
          withProgress(message = "Downloading",
                       detail = "It will take some time ^-^",
                       expr = {
                         for (i in 1:3){
                           incProgress(1/3)
                           Sys.sleep(0.5)
                         }
                         if (input$QQ_farmat == "jpeg"){
                           jpeg(file, height = input$QQ_height, width = input$QQ_width, units = "in", res = 800)
                         } else if (input$QQ_farmat == "pdf"){
                           pdf(file, height = input$QQ_height, width = input$QQ_width)
                         } else if (input$QQ_farmat == "png"){
                           png(file, height = input$QQ_height, width = input$QQ_width, units = "in", res = 800)
                         } else if (input$QQ_farmat == "svg"){
                           svg(file, height = input$QQ_height, width = input$QQ_width)
                         } else if (input$QQ_farmat == "tiff"){
                           tiff(file, height = input$QQ_height, width = input$QQ_width, units = "in", res = 800)
                         }
                         qqman::qq(GWAS_reslut_final$P.value,
                                   col = 'blue',
                                   cex = 2,
                                   las = 0.5,
                                   pch = 20)
                         dev.off()
                       }
          )}
      )
    })


    observeEvent(input$Plot_Manhattan, {
      withProgress(message = "Drawing",
                   detail = "It will take some time ^-^",
                   expr = {
                     for (i in 1:3){
                       incProgress(1/3)
                       Sys.sleep(0.5)
                     }

                     GWAS_reslut_final <- dplyr::select(GWAS_result$GWAS, !(2:3)) %>%
                       left_join(GWAS_Geno_num()$GM, by = c("SNP" = "SNP"))
                     yline = input$Manhattan_y
                     output$Manhattan_plot <- renderPlotly({
                       p_GWAS_reslut_final <- ggplot(GWAS_reslut_final, aes(x = Chromosome, y = -log10(P.value)))+
                         geom_jitter(aes(color=Chromosome), alpha = 0.8, size = 2)+
                         theme_bw()+
                         theme(legend.position = "none",
                               axis.text.x = element_text(angle=45,hjust=1),
                               panel.border = element_blank(),
                               axis.line.y = element_line(),
                               axis.line.x = element_line(),
                               panel.grid.major.x = element_blank(),
                               panel.grid.minor.x = element_blank())+
                         labs(x=NULL,y="-log10(P)")+
                         scale_y_continuous(expand = c(0,0)) +
                         geom_hline(yintercept = yline,lty="dashed")
                       plotly::ggplotly(p_GWAS_reslut_final)
                     })})

      output$download_Manhattan <- downloadHandler(
        filename = function(){
          paste0("Manhattan.",input$Manhattan_farmat, sep = "")
        },
        content = function(file){
          withProgress(message = "Downloading",
                       detail = "It will take some time ^-^",
                       expr = {
                         for (i in 1:3){
                           incProgress(1/3)
                           Sys.sleep(0.5)
                         }
                         if (input$Manhattan_farmat == "jpeg"){
                           jpeg(file, height = input$Manhattan_height, width = input$Manhattan_width, units = "in", res = 800)
                         } else if (input$Manhattan_farmat == "pdf"){
                           pdf(file, height = input$Manhattan_height, width = input$Manhattan_width)
                         } else if (input$Manhattan_farmat == "png"){
                           png(file, height = input$Manhattan_height, width = input$Manhattan_width, units = "in", res = 800)
                         } else if (input$Manhattan_farmat == "svg"){
                           svg(file, height = input$Manhattan_height, width = input$Manhattan_width)
                         } else if (input$Manhattan_farmat == "tiff"){
                           tiff(file, height = input$Manhattan_height, width = input$Manhattan_width, units = "in", res = 800)
                         }
                         p_Manhattan <- ggplot(GWAS_reslut_final, aes(x = Chromosome, y = -log10(P.value)))+
                           geom_jitter(aes(color=Chromosome), alpha = 0.8, size = 2)+
                           theme_bw()+
                           theme(legend.position = "none",
                                 axis.text.x = element_text(angle=45,hjust=1),
                                 panel.border = element_blank(),
                                 axis.line.y = element_line(),
                                 axis.line.x = element_line(),
                                 panel.grid.major.x = element_blank(),
                                 panel.grid.minor.x = element_blank())+
                           labs(x=NULL,y="-log10(P)")+
                           scale_y_continuous(expand = c(0,0)) +
                           geom_hline(yintercept = yline,lty="dashed")
                         plot(p_Manhattan)
                         dev.off()
                       }
          )}
      )
    })

    observeEvent(input$Extract, {
      withProgress(message = "Extracting",
                   detail = "It will take some time ^-^",
                   expr = {
                     for (i in 1:5){
                       incProgress(1/5)
                       Sys.sleep(0.5)
                     }
                     mygxf <- reactive({
                       req(input$uploadGXF)
                       gxf <- rtracklayer::import(input$uploadGXF$datapath)
                       return(gxf)
                     })
                     mygxf_df <- dplyr::filter(as.data.frame(mygxf()), type == "gene")
                     GWAS_reslut_final <- dplyr::select(GWAS_result$GWAS, !(2:3)) %>%
                       left_join(GWAS_Geno_num()$GM, by = c("SNP" = "SNP"))

                     Sig_p <- input$Significant_SNP
                     Sig_SNP <- dplyr::filter(GWAS_reslut_final, P.value <= Sig_p)
                     colnames(mygxf_df)[1] <- "Chromosome"
                     num_choice <- input$up_down_stream
                     Sig_SNP_range <- dplyr::mutate(Sig_SNP,
                                                    upstream = as.numeric(Sig_SNP$Position) - as.numeric(num_choice),
                                                    downstream = as.numeric(Sig_SNP$Position) + as.numeric(num_choice))
                     Gene_extr <- mygxf_df %>%
                       right_join(Sig_SNP_range, by = "Chromosome") %>%
                       dplyr::filter(start <= upstream | start <= downstream,
                                     end >= downstream | end >= upstream) %>%
                       distinct(ID, start, end, Chromosome, SNP, Position, .keep_all = TRUE) %>%
                       dplyr::select(ID, start, end, Chromosome, SNP, Position)
                     output$SNP_Gene <- renderDT(
                       datatable(Gene_extr, options = list(scrollX = TRUE))
                     )})
      output$download_SNP_Gene <- downloadHandler(
        filename = function(){
          paste0("SNP to Gene.", input$SNP_Gene_format, sep = "")
        },
        content = function(file){
          withProgress(message = "Downloading",
                       detail = "It will take some time ^-^",
                       expr = {
                         for (i in 1:3){
                           incProgress(1/3)
                           Sys.sleep(0.5)
                         }
                         if (input$SNP_Gene_format == "csv"){
                           write.csv(Gene_extr, file)
                         } else if (input$SNP_Gene_format == "tsv"){
                           write_tsv(Gene_extr, file)
                         } else if (input$SNP_Gene_format == "txt"){
                           write.table(Gene_extr, file)
                         }
                       }
          )}
      )
    })
  })
  ######################################################
  ######################################################



  ######################################################
  #######################Plot_LD########################
  observeEvent(input$Plot_LD, {
    withProgress(message = "Drawing",
                 detail = "It will take some time ^-^",
                 expr = {
                   for (i in 1:3){
                     incProgress(1/3)
                     Sys.sleep(0.5)
                   }
                   LD_hmp <- reactive({
                     req(input$uploadgeno)
                     LD_hmp_data <- read_delim(input$uploadgeno$datapath,
                                              col_types = cols(chrom = col_character()))
                     return(LD_hmp_data)
                   })
                   userchoice_chrom <- input$LD_Chromosome
                   leftpos <- input$LD_left
                   rightpos <- input$LD_right
                   LDSNPdata <- dplyr::filter(LD_hmp(), chrom == userchoice_chrom & pos >leftpos & pos < rightpos)
                   SNPdata <- dplyr::select(LDSNPdata, !(2:11))
                   SNPdata_t <- data.frame(t(SNPdata), row.names = NULL)
                   colnames(SNPdata_t) <-  SNPdata_t[1,]
                   SNPdata_df <- SNPdata_t[-1,]
                   num <- ncol(SNPdata_df)
                   for(i in 1:num){
                     SNPdata_df[,i] <- as.genotype(SNPdata_df[,i])
                   }
                   SNPpos <- filter(LD_hmp(), chrom == userchoice_chrom & pos >leftpos & pos < rightpos) %>%
                     pull(pos)
                   color.rgb <- colorRampPalette(rev(c("orange", "red")), space = "rgb")
                   output$LD_plot <- renderPlot({
                     LDheatmap(SNPdata_df, SNPpos, color = color.rgb(20), flip = T, SNP.name = colnames(SNPdata_df), geneMapLocation = 0.1)
                   })})

    output$download_LD <- downloadHandler(
      filename = function(){
        paste0("LD.",input$LD_farmat, sep = "")
      },
      content = function(file){
        withProgress(message = "Downloading",
                     detail = "It will take some time ^-^",
                     expr = {
                       for (i in 1:3){
                         incProgress(1/3)
                         Sys.sleep(0.5)
                       }
                       if (input$LD_farmat == "jpeg"){
                         jpeg(file, height = input$LD_height, width = input$LD_width, units = "in", res = 800)
                       } else if (input$LD_farmat == "pdf"){
                         pdf(file, height = input$LD_height, width = input$LD_width)
                       } else if (input$LD_farmat == "png"){
                         png(file, height = input$LD_height, width = input$LD_width, units = "in", res = 800)
                       } else if (input$LD_farmat == "svg"){
                         svg(file, height = input$LD_height, width = input$LD_width)
                       } else if (input$LD_farmat == "tiff"){
                         tiff(file, height = input$LD_height, width = input$LD_width, units = "in", res = 800)
                       }
                       LDheatmap(SNPdata_df, SNPpos, color = color.rgb(20), flip = T, SNP.name = colnames(SNPdata_df), geneMapLocation = 0.1)
                       dev.off()
                     }
        )}
    )
  })
  ######################################################
  ######################################################



  ######################################################
  ########################Run_GO########################
  observeEvent(input$Run_GO, {
    withProgress(message = "Calculating",
                 detail = "It will take some time ^-^",
                 expr = {
                   for (i in 1:5){
                     incProgress(1/5)
                     Sys.sleep(0.5)
                   }
                   GeneInfo <- reactive({
                     req(input$GeneInfo)
                     Info <- read.delim(input$GeneInfo$datapath)
                     return(Info)
                   })

                   gene_list <- reactive({
                     req(input$GeneList)
                     genelist <- read.delim(input$GeneList$datapath)
                     return(genelist)
                   })

                   Gene_list <- gene_list()$GID

                   Gene_info <- dplyr::select(GeneInfo(),
                                              GID = GID,
                                              GO = GO,
                                              KEGG = KEGG,
                                              Pathway = Pathway,
                                              Gene_Name = GID)
                   gene2go <- dplyr::select(Gene_info, GID, GO) %>%
                     separate_rows(GO, sep = ',', convert = F) %>%
                     dplyr::filter(!is.na(GO))

                   GOTerms <- dplyr::select(as.data.frame(GOTERM),
                                            GO = "go_id",
                                            Term = "Term",
                                            Ontology = "Ontology") %>%
                     right_join(gene2go, by = "GO") %>%
                     na.omit()

                   Go_class <- distinct(as.data.frame(GOTERM), go_id, Ontology, .keep_all =F)

                   GO_p <- input$GO_pvalue
                   #GO_q <- input$GO_qvalue

                   go_rich <- enricher(
                     gene = Gene_list,
                     TERM2GENE = GOTerms[c('GO', 'GID')],
                     TERM2NAME = GOTerms[c('GO', 'Term', 'Ontology')],
                     pvalueCutoff = GO_p,
                     #qvalueCutoff = GO_q,
                     minGSSize = 1,
                     maxGSSize = 100000000
                   )

                   GO_result_ori1 <- go_rich@result %>%
                     left_join(Go_class, by = c("ID" = "go_id"))
                   GO_result_ori2 <- dplyr::select(GO_result_ori1, -c(p.adjust, qvalue))
                   GO_result <- dplyr::filter(GO_result_ori2, pvalue <= GO_p)

                   output$GO_result <- renderDT({
                     datatable(GO_result, options = list(scrollX = TRUE))
                   })
                 })
    output$download_GO_table <- downloadHandler(
      filename = function(){
        paste0("GO_enrichment_result.", input$GO_table_format, sep = "")
      },
      content = function(file){
        withProgress(message = "Downloading",
                     detail = "It will take some time ^-^",
                     expr = {
                       for (i in 1:3){
                         incProgress(1/3)
                         Sys.sleep(0.5)
                       }
                       if (input$GO_table_format == "csv"){
                         write.csv(GO_result, file)
                       } else if (input$GO_table_format == "tsv"){
                         write_tsv(GO_result, file)
                       } else if (input$GO_table_format == "txt"){
                         write.table(GO_result, file)
                       }
                     }
        )}
    )

    observeEvent(input$Plot_GO, {
      withProgress(message = "Drawing",
                   detail = "It will take some time ^-^",
                   expr = {
                     for (i in 1:3){
                       incProgress(1/3)
                       Sys.sleep(0.5)
                     }
                     GO_result_BP <- dplyr::filter(GO_result, Ontology == "BP")
                     GO_result_CC <- dplyr::filter(GO_result, Ontology == "CC")
                     GO_result_MF <- dplyr::filter(GO_result, Ontology == "MF")

                     GO_result_BP_num = nrow(GO_result_BP)
                     GO_result_CC_num = nrow(GO_result_CC)
                     GO_result_MF_num = nrow(GO_result_MF)

                     user_BP_num = input$BP_num
                     user_CC_num = input$CC_num
                     user_MF_num = input$MF_num

                     if (user_BP_num <= GO_result_BP_num){
                       BP_num =  user_BP_num
                     }else{
                       BP_num = GO_result_BP_num
                     }

                     if (user_CC_num <= GO_result_CC_num){
                       CC_num =  user_CC_num
                     }else{
                       CC_num = GO_result_CC_num
                     }

                     if (user_MF_num <= GO_result_MF_num){
                       MF_num =  user_MF_num
                     }else{
                       MF_num = GO_result_MF_num
                     }

                     GO_BP <- subset(GO_result, subset = (Ontology == "BP"))[1:BP_num,]
                     GO_CC <- subset(GO_result, subset = (Ontology == "CC"))[1:CC_num,]
                     GO_MF <- subset(GO_result, subset = (Ontology == "MF"))[1:MF_num,]

                     GO_DF <- rbind(GO_BP, GO_CC, GO_MF)

                     output$GO_plot <- renderPlotly({
                       p_GO <- ggplot(GO_DF, aes(GeneRatio, Description)) +
                         geom_point(aes(size = Count, color = -log10(pvalue)))+
                         scale_color_gradient(low = "green", high = "red") +
                         facet_wrap(~ Ontology, ncol = 1, scales = "free") +
                         theme_bw()
                       plotly::ggplotly(p_GO)
                     })
                   })
      output$download_GO <- downloadHandler(
        filename = function(){
          paste0("GO_enrichment.",input$GO_farmat, sep = "")
        },
        content = function(file){
          withProgress(message = "Downloading",
                       detail = "It will take some time ^-^",
                       expr = {
                         for (i in 1:3){
                           incProgress(1/3)
                           Sys.sleep(0.5)
                         }
                         if (input$GO_farmat == "jpeg"){
                           jpeg(file, height = input$GO_height, width = input$GO_width, units = "in", res = 800)
                         } else if (input$GO_farmat == "pdf"){
                           pdf(file, height = input$GO_height, width = input$GO_width)
                         } else if (input$GO_farmat == "png"){
                           png(file, height = input$GO_height, width = input$GO_width, units = "in", res = 800)
                         } else if (input$GO_farmat == "svg"){
                           svg(file, height = input$GO_height, width = input$GO_width)
                         } else if (input$GO_farmat == "tiff"){
                           tiff(file, height = input$GO_height, width = input$GO_width, units = "in", res = 800)
                         }
                         p_GO <- ggplot(GO_DF, aes(GeneRatio, Description)) +
                           geom_point(aes(size = Count, color = -log10(pvalue)))+
                           scale_color_gradient(low = "green", high = "red") +
                           facet_wrap(~ Ontology, ncol = 1, scales = "free") +
                           theme_bw()
                         plot(p_GO)
                         dev.off()
                       })
        })
    })
  })
  ######################################################
  ######################################################



  ######################################################
  #######################Run_KEGG#######################
  observeEvent(input$Run_KEGG, {
    withProgress(message = "Calculating",
                 detail = "It will take some time ^-^",
                 expr = {
                   for (i in 1:5){
                     incProgress(1/5)
                     Sys.sleep(0.5)
                   }
                   GeneInfo <- reactive({
                     req(input$GeneInfo)
                     Info <- read.delim(input$GeneInfo$datapath)
                     return(Info)
                   })

                   gene_list <- reactive({
                     req(input$GeneList)
                     genelist <- read.delim(input$GeneList$datapath)
                     return(genelist)
                   })

                   Gene_list <- gene_list()$GID

                   Gene_info <- dplyr::select(GeneInfo(),
                                              GID = GID,
                                              GO = GO,
                                              KEGG = KEGG,
                                              Pathway = Pathway,
                                              Gene_Name = GID)
                   pathway2gene <- dplyr::select(Gene_info, Pathway, GID) %>%
                     separate_rows(Pathway, sep = ',', convert = F) %>%
                     filter(str_detect(Pathway, 'ko')) %>%
                     mutate(Pathway = str_remove(Pathway, 'ko'))

                   get_path2name <- function(){
                     keggpathid2name.df <- clusterProfiler:::kegg_list("pathway")
                     keggpathid2name.df[,1] %<>% gsub("path:map", "", .)
                     colnames(keggpathid2name.df) <- c("path_id", "path_name")
                     return(keggpathid2name.df)}

                   pathway2name <- get_path2name() %>%
                     mutate(path_id = str_remove(path_id, 'map'))

                   KEGG_p <- input$KEGG_pvalue
                   KEGG_q <- input$KEGG_qvalue
                   kegg_rich <- enricher(
                     gene = Gene_list,
                     TERM2GENE = pathway2gene,
                     TERM2NAME = pathway2name,
                     pvalueCutoff = KEGG_p,
                     #qvalueCutoff = KEGG_q,
                     minGSSize = 1,
                     maxGSSize = 100000000
                   )
                   kegg_result_ori1 <- dplyr::select(kegg_rich@result, -c(p.adjust, qvalue))
                   kegg_result <- dplyr::filter(kegg_result_ori1, pvalue <= KEGG_p)

                   output$KEGG_result <- renderDT(
                     datatable(kegg_result, options = list(scrollX = TRUE))
                   )})

    output$download_KEGG_table <- downloadHandler(
      filename = function(){
        paste0("KEGG_enrichment_result.",input$KEGG_table_format, sep = "")
      },
      content = function(file){
        withProgress(message = "Downloading",
                     detail = "It will take some time ^-^",
                     expr = {
                       for (i in 1:3){
                         incProgress(1/3)
                         Sys.sleep(0.5)
                       }
                       if (input$blup_table_format == "csv"){
                         write.csv(kegg_result, file)
                       } else if (input$blup_table_format == "tsv"){
                         write_tsv(kegg_result, file)
                       } else if (input$blup_table_format == "txt"){
                         write.table(kegg_result, file)
                       }
                     })
      })

    observeEvent(input$Plot_KEGG, {
      withProgress(message = "Drawing",
                   detail = "It will take some time ^-^",
                   expr = {
                     for (i in 1:5){
                       incProgress(1/5)
                       Sys.sleep(0.5)
                     }
                     user_kegg_num = input$KEGG_num
                     kegg_num = nrow(kegg_rich@result)
                     if (user_kegg_num <= kegg_num){
                       KEGG_Num =  user_kegg_num
                     }else{
                       KEGG_Num = kegg_num
                     }
                     kegg_DF <- subset(kegg_rich@result)[1:KEGG_Num,]

                     output$KEGG_plot <- renderPlotly({
                       p_KEGG <- ggplot(kegg_DF, aes(GeneRatio, Description)) +
                         geom_point(aes(size = Count, color = -log10(pvalue)))+
                         scale_color_gradient(low = "green", high = "red") +
                         theme_bw()
                       plotly::ggplotly(p_KEGG)
                     })
                   })
      output$download_KEGG <- downloadHandler(
        filename = function(){
          paste0("KEGG_enrichment.",input$KEGG_farmat, sep = "")
        },
        content = function(file){
          withProgress(message = "Downloading",
                       detail = "It will take some time ^-^",
                       expr = {
                         for (i in 1:3){
                           incProgress(1/3)
                           Sys.sleep(0.5)
                         }
                         if (input$KEGG_farmat == "jpeg"){
                           jpeg(file, height = input$KEGG_height, width = input$KEGG_width, units = "in", res = 800)
                         } else if (input$KEGG_farmat == "pdf"){
                           pdf(file, height = input$KEGG_height, width = input$KEGG_width)
                         } else if (input$KEGG_farmat == "png"){
                           png(file, height = input$KEGG_height, width = input$KEGG_width, units = "in", res = 800)
                         } else if (input$KEGG_farmat == "svg"){
                           svg(file, height = input$KEGG_height, width = input$KEGG_width)
                         } else if (input$KEGG_farmat == "tiff"){
                           tiff(file, height = input$KEGG_height, width = input$KEGG_width, units = "in", res = 800)
                         }
                         p_KEGG <- ggplot(kegg_DF, aes(GeneRatio, Description)) +
                           geom_point(aes(size = Count, color = -log10(pvalue)))+
                           scale_color_gradient(low = "green", high = "red") +
                           theme_bw()
                         plot(p_KEGG)
                         dev.off()
                       })})

    })


  })
  ######################################################
  ########################end###########################
}
