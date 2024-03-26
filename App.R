library(shiny)
options(shiny.maxRequestSize = 10 * 1024^3)
library(colourpicker)
library(bslib)
library(readxl)
library(shinydashboard)
library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(scales)
library(RColorBrewer)
library(tibble)
library(purrr)
library(DT)
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DESeq2")
library(DESeq2)
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("fgsea")
library(fgsea)

ui <- fluidPage(
tabsetPanel(
  tabPanel("Samples",
      sidebarLayout(
        sidebarPanel(
              fileInput("file", "Upload a CSV File", accept = c(".xlsx")),
              radioButtons("select_x_axis", "Select Column for x axis:", 
                          choices = c("Treatment", "Sex", "Timepoint", "Lifestage"),
                          selected = "Sex")),
        mainPanel(
              tabsetPanel(
              tabPanel("Summary", tableOutput("summary")),
              tabPanel("Table", dataTableOutput("table")),
              tabPanel("Plot", plotOutput("histogram")))))),
  tabPanel("Counts",
      sidebarLayout(
        sidebarPanel(
              fileInput("file1", "Upload Normalized Counts File", accept = c(".csv")),
              sliderInput("variance_slider", "Variance Percentile", min = 0, 
                          max = 100, value = 50),
              sliderInput("nonzero_slider", "Non-zero Samples Threshold", min = 1, 
                              max = 80, value = 10)),
        mainPanel(
              tabsetPanel(
              tabPanel("Summary", h2("Summary of Normalized Counts Data"), 
                       textOutput("Countstext1"), textOutput("Countstext2"),
                       textOutput("Countstext3"), textOutput("Countstext4"),
                       textOutput("Countstext5"), textOutput("Countstext6")),
              tabPanel("Scatter Plots", 
                       plotOutput("scatterA"), plotOutput("scatterB")),
              tabPanel("Heatmap", plotOutput("heatmap")),
              tabPanel("PCA Plot", plotOutput("pca"), 
                       radioButtons("pc_x", "Select X-axis PC: ", 
                                    choices = c("PC1", "PC2", "PC3", "PC4"),
                                    selected = "PC1"),
                       radioButtons("pc_y", "Select Y-axis PC: ",
                                    choices = c("PC1", "PC2", "PC3", "PC4"),
                                    selected = "PC2"))
              )))),
  tabPanel("Differential Expression",
     sidebarLayout(
       sidebarPanel(
             fileInput("file2", "Upload Differential Expression File", accept = c(".csv")),
             radioButtons("select_x_axis_DE", "Select Column for x axis:", 
                          choices = c("baseMean", "log2FoldChange", "lfcSE", 
                               "stat", "pvalue", "padj"),
                          selected = "log2FoldChange"),
             radioButtons("select_y_axis_DE", "Select Column for y axis:", 
                          choices = c("baseMean", "log2FoldChange", "lfcSE", 
                               "stat", "pvalue", "padj"),
                          selected = "padj"),
             sliderInput("p_adjust_slider", "Adjust p-Value Threshold:", min = -300, 
                          max = -1, value = -100),
             colourInput("color_above", "Base point color", value = "red"),
             colourInput("color_below", "Highlight point color", value = "blue")),
        mainPanel(
             tabsetPanel(
             tabPanel("Differential Expression Table", dataTableOutput("DEtable")),
             tabPanel("Volcano Plot", plotOutput("DEvolcano")))))),
  tabPanel( "GSEA",
    tabsetPanel(
    tabPanel("Table",
       sidebarLayout(
          sidebarPanel(
             fileInput("file3", "Upload FGSEA File", accept = c(".csv")),
             sliderInput("p_adjust_sliderfgsea", "Adjust p-Value:", 
                        min = 0, max = 1, value = c(0, 0.00003)),
             radioButtons("nes_selection", "Select NES:",
                         choices = c("All", "Positive", "Negative"), selected = "All"),
             downloadButton("downloadTable", "Download Filtered Table")),
           mainPanel(dataTableOutput("tablefgsea")))),
     tabPanel("Barplot",
        sidebarLayout(
           sidebarPanel(
              fileInput("file4", "Upload FGSEA File", accept = c(".csv")),
              sliderInput("p_adjust_sliderfgsea_barplot", "Adjust p-Value:", 
                         min = 0, max = 1, value = c(0, 0.05))),
           mainPanel(plotOutput("barplotfgsea")))),
     tabPanel("Scatterplot",
        sidebarLayout(
           sidebarPanel(
              fileInput("file5", "Upload FGSEA File", accept = c(".csv")),
              sliderInput("threshold_sliderfgsea_scatterplot", "Adjust p-Value:", 
                              min = 0, max = 1, value = c(0, 0.05))),
           mainPanel(plotOutput("scatterplot_fgsea"))))))))

server <- function(input, output, session){

##############################  SAMPLES FUNCTION  ##############################
#Samples Data
samples_data <- reactive({
  req(input$file)
  samples_dataframe <- read_excel(input$file$datapath)
  return(samples_dataframe)})
  
#Samples Function - Summary
summary_samples_data <- reactive({
  data <- samples_data()[,-1:-2]
  summary_data <- data.frame(
      Column_Name = names(data),
      Type = sapply(data, function(x) if (is.character(x)) "Factor" else if (is.numeric(x)) "Numeric" else "Other"),
      Distinct_Values = sapply(data, function(x) if (is.character(x)) toString(unique(x)) else ""),
      stringsAsFactors = FALSE)
  return(summary_data)})

#Summary Tab Output
output$summary <- renderTable({summary_samples_data()})
 
       

#Samples Function - Table
table_samples_data <- reactive({
  return(samples_data())})

#Table Tab Output
output$table <- renderDataTable({
  datatable(table_samples_data(),options = list(ordering = TRUE))})
        


#Samples Function - Plot
plot_samples_data <- reactive({
  x_name <- input$select_x_axis
  hist_samples_data <- ggplot(samples_data(), aes(x = .data[[x_name]])) + 
                       geom_bar() +
                       labs(title = "Histogram for Samples", x = x_name)
  return(hist_samples_data)})
  
#Plot Tab
output$histogram <- renderPlot({plot_samples_data()})

################################################################################
      
          
          
              
#########################  NORMALIZED COUNTS FUNCTION  #########################       
#Normalized Counts Data   
load_data_norm_counts <- reactive({
  req(input$file1)
  data_counts <- read.csv(input$file1$datapath, header = TRUE)
  return(data_counts)})

#Filtered Counts Data 
filtered_data <- reactive({
  data <- load_data_norm_counts()
  data <- data[-1, -1]
  variance_threshold <- quantile(apply(data, 1, var), input$variance_slider/100)
  data <- data[rowSums(data^2) >= variance_threshold, ]
  data <- data[rowSums((data != 0)) >= input$nonzero_slider, ]
  return(data)})

#Text Output 1
output$Countstext1 <- renderText({
  data <- load_data_norm_counts()[,-1]
  filtered <- filtered_data()
  paste("Number of Samples: ", ncol(data))})

#Text Output 2
output$Countstext2 <- renderText({
  data <- load_data_norm_counts()[-1,]
  filtered <- filtered_data()
  paste("Total Number of genes: ", nrow(data))})

#Text Output 3
output$Countstext3 <- renderText({
  data <- load_data_norm_counts()
  filtered <- filtered_data()
  paste("Number of Genes Passing the Current Filter: ", nrow(filtered))})

#Text Output 4
output$Countstext4 <- renderText({
  data <- load_data_norm_counts()[-1,]
  filtered <- filtered_data()
  paste("Percentage of Genes Passing the Filter: ", percent(nrow(filtered)/nrow(data)))})

#Text Output 5
output$Countstext5 <- renderText({
  data <- load_data_norm_counts()[-1,]
  filtered <- filtered_data()
  paste("Number of Genes Not Passing the Current Filter: ", nrow(data) - nrow(filtered))})

#Text Output 6
output$Countstext6 <- renderText({
  data <- load_data_norm_counts()[-1,]
  filtered <- filtered_data()
  paste("Percentage of Genes Not Passing the Filter: ", 
        percent((nrow(data) - nrow(filtered))/nrow(data)))})



# Counts Function - Scatter Plot (Median vs Variance)
scatterplot_medianvsvar <- function(data, variance_slider){
  data <- data[-1, -1] 
  median_count <- apply(data, 1, median)
  variance <- apply(data, 1, var)
  variance_threshold <- quantile(apply(data, 1, var), variance_slider/100)
  plot_data <- data.frame(MedianCount = median_count, Variance = variance,
               GeneFilter = ifelse(rowSums(data^2) >= variance_threshold, 'Pass', 'Fail'))
  plot1 <- ggplot(plot_data, aes(x = log10(MedianCount), y = log10(Variance), color = GeneFilter)) +
          geom_point() +
          scale_color_manual(values = c('Fail' = 'lightgrey', 'Pass' = 'darkblue')) +
          labs(title = 'Median Count vs Variance (Log10 Scale)', 
               x = 'Median Count', y = 'Variance') +
          theme_minimal()
  return(plot1)}

# Counts Function - Scatter Plot (Median vs Zero Genes)
scatterplot_medianvszero <- function(data, nonzero_slider){
  data <- data[-1, -1] 
  median_count <- apply(data, 1, median)
  num_zeros <- apply(data == 0, 1, sum)
  plot_data <- data.frame(MedianCount = median_count, NumZeros = num_zeros, 
                    GeneFilter = ifelse(rowSums((data != 0)) >= nonzero_slider, 'Pass', 'Fail'))
  plot2 <- ggplot(plot_data, aes(x = log10(MedianCount), y = log10(NumZeros), color = GeneFilter)) +
    geom_point() +
    scale_color_manual(values = c('Fail' = 'lightpink', 'Pass' = 'darkred')) +
    labs(title = 'Median Count vs Number of Zeros (Log10 Scale)', 
         x = 'Median Count', y = 'Number of Zeros') +
    theme_minimal()
  return(plot2)}

#Plot Output
output$scatterA <- renderPlot({
  scatterplot_medianvsvar(load_data_norm_counts(), input$variance_slider)})
     
#Plot Output
output$scatterB <- renderPlot({
  scatterplot_medianvszero(load_data_norm_counts(), input$nonzero_slider)})



#Counts Function - Heatmap
plot_heatmap <- function(filtered_data, num_colors, palette){
  filtered_data <- as.matrix(filtered_data)
  log_data <- log10(filtered_data)
  heatmap(filtered_data, 
          col = colorRampPalette(brewer.pal(num_colors, palette))(num_colors),
          xlab = 'Samples', ylab = 'Genes', scale = 'row', 
          main = "\nHeatmap of Filtered Genes\n", Rowv = NA, Colv = NA)
  legend(x = "bottomright", legend = c("low", "medium", "high"), 
         cex = 0.8, fill = colorRampPalette(brewer.pal(num_colors, palette))(3))}

#Heatmap Output
output$heatmap <- renderPlot({plot_heatmap(filtered_data(), 9, 'PuBuGn')})



#Counts Function - PCA
pca_objects <- reactive({
  # Keep the selected columns
  data <- load_data_norm_counts()[-1,-1]
  scaled_df <- scale(data)
  pca_output <- prcomp(na.omit(scaled_df))
  # data.frame of PCs
  pcs_df <- data.frame(scaled_df, pca_output$x)
  return(list(pca_output = pca_output, 
              pcs_df = pcs_df))})

#PCA Plot
pca_biplot <- reactive({
  pcs_df <- pca_objects()$pcs_df
  pca_output <-  pca_objects()$pca_output
  var_expl_x <- round(100 * pca_output$sdev[as.numeric(gsub("[^0-9]", "", input$pc_x))]^2/sum(pca_output$sdev^2), 1)
  var_expl_y <- round(100 * pca_output$sdev[as.numeric(gsub("[^0-9]", "", input$pc_y))]^2/sum(pca_output$sdev^2), 1)
  PCx <- input$pc_x
  PCy <- input$pc_y
  pc_plot  <- ggplot(pcs_df, aes(x = get(PCx), y = get(PCy))) +
      geom_point() + theme_bw(base_size = 14) + coord_equal() +
      xlab(paste0(input$pc_x, " (", var_expl_x, "% explained variance)")) +
      ylab(paste0(input$pc_y, " (", var_expl_y, "% explained variance)")) 
    return(pc_plot)})

#PCA Plot Output
output$pca <- renderPlot({pca_biplot()})
################################################################################
          
          
          
          
######################  DIFFERENTIAL EXPRESSION FUNCTION  ######################
#Differential Expression Data    
load_data_de <- reactive({
  req(input$file2)
  data <- read.csv(input$file2$datapath, header = TRUE)
  return(data)})

#Table Tab Output
output$DEtable <- renderDataTable({
  datatable(load_data_de(), options = list(ordering = TRUE))})
          
          
          
#Differential Expression - Volcano Plot
volcano_plot <- function(dataf, x_name_de, y_name_de, slider, color1, color2) {
  ggplot(dataf, aes(x = !!sym(x_name_de), y = -log10(!!sym(y_name_de)), 
                    col = padj < 1* 10^slider)) +
                    scale_color_manual(values = c(color1, color2)) +
                    geom_point() +
                    labs(title = "Volcano Plot", x = x_name_de, y = '-log10padj')}

#Volcano Plot Output
output$DEvolcano <- renderPlot({
  volcano_plot(load_data_de(), input$select_x_axis_DE, input$select_y_axis_DE, 
               input$p_adjust_slider, input$color_above, input$color_below)})
################################################################################



###############################  FGSEA FUNCTION  ###############################
#FGSEA Data     
load_data_fgsea <- reactive({
  req(input$file3)
  data <- read_csv(input$file3$datapath)
  filtered_data <- subset(data, padj >= input$p_adjust_sliderfgsea[1] & padj <= input$p_adjust_sliderfgsea[2])
  if (input$nes_selection == "Positive") {
    filtered_data <- subset(filtered_data, NES >= 0)
  } else if (input$nes_selection == "Negative") {
    filtered_data <- subset(filtered_data, NES < 0)
  }
  return(filtered_data)})

#Table Tab Output
output$tablefgsea <- renderDataTable({datatable(load_data_fgsea(),
                                      options = list(ordering = TRUE))})

#Download Table Output
output$downloadTable <- downloadHandler(
  filename = function() {
    paste("filtered_table_", Sys.Date(), ".csv", sep = "")},
  content = function(file) {
    write.csv(load_data_fgsea(), file, row.names = FALSE)})



#FGSEA Function - Bar Plot Data
load_data_barplot <- reactive({
  req(input$file4)
  data <- read_csv(input$file4$datapath)
  # Filter based on p-adjusted value
  filtered_data <- subset(data, padj >= input$p_adjust_sliderfgsea_barplot[1] & padj <= input$p_adjust_sliderfgsea_barplot[2])
  # You can add additional filtering logic specific to the barplot here
  return(filtered_data)
})

#FGSEA Function - Bar Plot
top_pathways <- function(fgsea_results, padj_slider) {
  # Order the data by padj in ascending order and select the top 10 rows
  top_nes <- fgsea_results[order(fgsea_results$padj), ]
  # Create the barplot using ggplot2
  barplot_fgsea <- ggplot(top_nes, aes(x = padj, y = reorder(pathway, padj < padj_slider))) +
    geom_bar(stat = "identity") +
    labs(title = "Top Pathways by Adjusted p-value", x = "p Adjusted", y = "Pathway") + 
    theme_minimal()
  return(barplot_fgsea)}

#Bar Plot Output
output$barplotfgsea <- renderPlot({
  top_pathways(load_data_barplot(), input$p_adjust_sliderfgsea_barplot)
})



#FGSEA Function - Scatterplot Data
load_data_scatterplot <- reactive({
  req(input$file5)
  threshold <- input$threshold_sliderfgsea_scatterplot
  data <- read_csv(input$file5$datapath)
  data$color <- ifelse(data$padj < threshold, "grey", "red")
  return(data)
})

#FGSEA Function - Scatterplot
scatter_plot_fgsea <- function(data){
  scatter_plot <- ggplot(data, aes(x = NES, y = -log10(padj), 
                          color = color)) +
    geom_point() +
    labs(title = "Scatter Plot of NES vs. -log10 Adjusted p-value",
         x = "NES",
         y = "-log10 Adjusted p-value") +
    theme_minimal() +
  scale_color_manual(values = c("grey" = "grey", "red" = "red"))
  return(scatter_plot)
}

#Scatter Plot Output
output$scatterplot_fgsea <- renderPlot({
  scatter_plot_fgsea(load_data_scatterplot())})}
################################################################################

shinyApp(ui = ui, server = server)

