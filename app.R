library(shiny)
library(scales)
library(DT)
library(bslib)
library(shinyjs)
library(dplyr)
library(colourpicker)
library(tidyr)
library(htmltools)
library(ggplot2)
library(gplots)
library(DESeq2)
library(gridExtra)
library(plotly)

################################################################################

#Making a clean sample metadata file
clean_metadata <- function(data) {
  data$Treatment <- sub('treatment: ', '', data$Treatment)
  data$Sex <- sub('Sex: ', '', data$Sex)
  data$Timepoint <- sub('timepoint: ', '', data$Timepoint)
  data$Lifestage <- sub('lifestage: ', '', data$Lifestage)
  return(data)
}
###################################################################

#Making a custom summary
custom_summary <- function(metadata) {
  # Number of rows and columns
  num_rows <- nrow(metadata)
  num_cols <- ncol(metadata)

  # Create an empty data frame to store the results
  result_df <- data.frame(Column_Names = character(0), Type = character(0), Distinct_Values = character(0), stringsAsFactors = FALSE)

  # Print column summary
  for (col_name in colnames(metadata)) {
    col_type <- class(metadata[[col_name]])

    if (col_type == "numeric") {
      mean_val <- mean(metadata[[col_name]])
      sd_val <- sd(metadata[[col_name]])
      values <- paste0(format(mean_val, digits = 3), " (+/- ", format(sd_val, digits = 3), ")")
    } else {
      values <- paste(unique(metadata[[col_name]]), collapse = ", ")
    }

    # Append the results to the data frame
    result_df <- rbind(result_df, data.frame(Column_Names = col_name, Type = col_type, Distinct_Values = values, stringsAsFactors = FALSE))
  }

  # Remove the underscores in the column names
  colnames(result_df) <- gsub("_", " ", colnames(result_df))

  # Return the final data frame, excluding the first row (column names)
  return(result_df[-1, ])
}

###############################################################################

#Making a plot of the no of characteristics
create_grouped_bar_plot <- function(data) {
  # Reshape the data to long format
  data_long <- data %>%
    pivot_longer(cols = c(Treatment, Sex, Timepoint, Lifestage), names_to = "Variable")

  # Create a summary data frame for labels
  label_data <- data_long %>%
    group_by(Variable, value) %>%
    summarise(Count = n())

  # Create a grouped bar plot with larger bars and labels inside the bars
  plot_result <- ggplot(data_long, aes(x = Variable, fill = value)) +
    geom_bar(position = "dodge", width = 0.7) +  # Adjust the width as needed
    geom_text(
      data = label_data,
      aes(label = value, y = Count, group = value),
      position = position_dodge(width = 0.7),
      vjust = 0.5,
      hjust = -0.3,  # Center the labels inside the bars
      size = 5,
      angle = 270
    ) +
    geom_text(
      data = label_data,
      aes(label = Count, y = Count, group = value),
      position = position_dodge(width = 0.7),
      vjust = -0.5,  # Adjust the vertical justification for counts on top
      size = 5,
      hjust = 0.5  # Center the counts on top horizontally
    ) +
    labs(x = "Variables", y = "Count", fill = "") +
    theme_bw() +
    ggtitle("Sample Attributes") +
    theme(plot.title = element_text(size = 18, hjust = 0.5, face = "bold"))

  return(plot_result)
}
##################################################################################
#Generating a Heatmap
generate_heatmap <- function(data, title = "Clustered Heatmap", log_transform = TRUE) {
  if (log_transform) {
    data <- log2(data + 1)  # Log transformation with pseudocount
  }

  heatmap.2(
    as.matrix(data),
    Rowv = TRUE,
    Colv = TRUE,
    dendrogram = "both",
    trace = "none",
    col = colorRampPalette(c("yellow", "#22577A"))(100),
    margins = c(5, 10),
    main = title,
    xlab = "Samples",
    ylab = "Genes",
    key = TRUE,
    keysize = 1.5,
    density.info = "none",
    cexCol = 0.8,
    cexRow = 0.8
  )
}
################################################################################
max_pcs <- function(data) {
  num_samples <- nrow(data)
  num_features <- ncol(data)
  return(min(num_samples, num_features))
}
################################################################################
perform_pca <- function(data_v, num_pcs = max_pcs(data_v)) {
  # Ensure that data is a numeric matrix or data frame
  if (!is.numeric(data_v)) {
    stop("Input data must be a numeric matrix or data frame.")
  }

  # Perform PCA
  pca_result <- prcomp(data_v, scale. = TRUE)

  # Extract PC scores
  pc_scores <- pca_result$x[, 1:num_pcs]

  # Extract percentage of variance explained by each PC
  variance_explained <- (pca_result$sdev^2) / sum(pca_result$sdev^2) * 100

  # Create data frame for PCA scores
  pca_data <- data.frame(pc_scores)

  # Return list with PCA data and percentage of variance explained
  return(list(pca_data = pca_data, variance_explained = variance_explained))
}

plot_pca <- function(pca_data, variance_explained, pc_x, pc_y) {
  # Print information about the input data
  print(paste("Number of rows in pca_data:", nrow(pca_data)))
  print(paste("Number of columns in pca_data:", ncol(pca_data)))

  pc_x <- as.numeric(pc_x)
  pc_y <- as.numeric(pc_y)

  plot_data <- data.frame(
    Sample = rownames(pca_data),
    PC1 = pca_data[, pc_x],
    PC2 = pca_data[, pc_y]
  )

  # Print information about the plot_data
  print(head(plot_data))  # Print the first few rows of plot_data
  print(round(variance_explained[pc_x], digits = 2))
  print(round(variance_explained[pc_y], digits = 2))

  # Plot PCA results
  ggplot(plot_data, aes(x = PC1, y = PC2, color = Sample)) +
    geom_point() + ggtitle("Principal Component Analysis") +
    theme(plot.title = element_text(size = 18, hjust = 0.5, face = "bold"))+
    labs(
      x = paste("PC", pc_x, ": ", round(variance_explained[pc_x], 2), "%", sep = ""),
      y = paste("PC", pc_y, ": ", round(variance_explained[pc_y], 2), "%", sep = "")
      )
}
################################################################################
# Filter zero variance genes

filter_zero_var_genes <- function(counts_matrix) {
  b <- apply(counts_matrix, 1, var)
  c <- counts_matrix[b > 0, ]
  return(c)
}
################################################################################

ui <- fluidPage(
  tags$head(
    tags$style(
      HTML("
        body {
          background-color: AliceBlue ; /* Replace with your desired color code */
          margin: 0; /* Remove default margins */
          height: 80vh; /* 100% of the viewport height */
        }

        .center-title {
          text-align: center;
          width: 100%;
        }
        .subtitle {
          text-align: center;
          font-style: italic;
          color: #666;
        }
        .main-panel {
          width: 95%; /* Set the width to 95% */
        }

      ")
    )
  ),
  div(class = "center-title", titlePanel("DrosaGenXplorer", "title")),
  div(class = "subtitle", "Genomic Exploration of Drosophila melanogaster in Dynamic Thermal Environments (DOI: 10.1111/mec.16463)"),

  mainPanel(
    tabsetPanel(
      tabPanel("Samples",
               sidebarLayout(
                 sidebarPanel(width = 3,
                              helpText("Use this tab if you want to explore your sample metadata."),
                              fileInput("upload_samples", 'Load your sample file', accept = c(".csv", ".tsv"))
                 ),
                 mainPanel(width = 9,
                   tabsetPanel(
                     tabPanel('Summary', tableOutput("summary_output"), verbatimTextOutput("summary_info")),
                     tabPanel('Table', dataTableOutput("table_tab1")),
                     tabPanel('Plot', plotOutput("plot_tab1", height = "550px"))
                   )
                 )
               )
      ),
      tabPanel("Counts",
               sidebarLayout(
                 sidebarPanel(width = 3,
                              helpText("Use this tab if you want to explore your normalized counts matrix."),
                              HTML('<p><strong style="color: #EA1751 ;">Note:</strong> If you think your counts data is not normalized, use the De-Mystify tab.</p>'),
                              fileInput("upload_counts", 'Load your normalized counts matrix', accept = c(".csv")),
                              sliderInput('slider_v', "Percentile of variance", min = 1, max = 100, value =50),
                              sliderInput('slider_nz', "Non-zero samples", min = 0, max = 80, value = 40)
                              ),
                 mainPanel(width = 9,
                   tabsetPanel(
                     tabPanel('Summary', tableOutput("summary_table")),
                     tabPanel('ScatterPlot', plotOutput("scatter_plots", height = "550px")),
                     tabPanel('Heatmap', plotOutput('heatmap', height = "550px")),
                     tabPanel('PCA Projections',
                              sidebarLayout(
                                sidebarPanel( width = 2,
                                  # Dropdown menus for selecting principal components
                                  selectInput("pc_x", "X-axis (PC):", choices = 1:80),
                                  selectInput("pc_y", "Y-axis (PC):", choices = 1:80)
                                ),
                                mainPanel( width = 10,
                                  plotOutput("pca_plot", height = "550px")
                   )
                 )
               )
                   )
                 )
               )
               ),

      tabPanel("Differential Expression",
               sidebarLayout(
                 sidebarPanel(width = 3,
                              helpText("Use this tab if you want to explore your differential expression results."),
                              HTML('<p><strong style="color:red;">Note:</strong> If you don\'t have DE results, use the De-Mystify tab.</p>'),
                              fileInput("upload_de_result", 'Load your differential expression analysis results', accept = c(".csv"))
                 ),
                 mainPanel(width = 9,
                   tabsetPanel(
                     tabPanel("Differential Expression Results",
                              dataTableOutput('DE_data')),
                     tabPanel("Visualizing DE Results",
                      sidebarLayout(
                          sidebarPanel(width = 3,
                                      helpText(HTML('A volcano plot can be generated with "<strong>log2 fold-change</strong>" on the x-axis and "<strong>p-adjusted</strong>" on the y-axis.')),
                                      radioButtons('x_name', 'Choose the column for X-axis', choices = c('baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj'), selected = 'log2FoldChange'),
                                      radioButtons('y_name', 'Choose the column for Y-axis', choices = c('baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj'), selected = 'padj'),
                                      colourInput('color1', 'Base point color', value = '#22577A'),
                                      colourInput('color2', 'Highlight point color', value = '#FFCF56'),
                                      sliderInput('slider_p', "Magnitude of p-adjusted coloring", min = -300, max = 0, value = -150, step = 1),
                                      actionButton(inputId = 'Plot_v', label = 'Plot', class = "btn-block", style = "color: black; background-color: darkorange;")
                                             ),

                                           mainPanel(width = 9,
                                             tabsetPanel(
                                               tabPanel("Volcano Plot",
                                                        plotOutput('volcano')),
                                               tabPanel("Table",
                                                        DT::dataTableOutput('table_output'))
                                             )
                                           )
                                         )
                     )
                              )
                     )
                   )
                 ),
      tabPanel("Gene Set Enrichment Analysis",
               sidebarLayout(
                  sidebarPanel(width = 3,
                               helpText("Use this tab if you want to explore your GSEA results."),
                               HTML('<p><strong style="color: #EA1751 ;">Note:</strong> If you don\'t have GSEA results, use the De-Mystify tab.</p>'),
                               fileInput("upload_gsea", 'Load your GSEA results', accept = c(".csv", ".tsv")),
                              actionButton(inputId = 'Upload_Results', label = 'Upload', class = "btn-block", style = "color: black; background-color: darkorange;")
                              ),
                            mainPanel(width = 9,
                              tabsetPanel(
                                tabPanel("Barplots",
                                  sidebarLayout(
                                    sidebarPanel(width = 3,
                                        sliderInput('p_adj', label = 'p-adj filter for Top pathways', min = 0, max = 1, value = 0.5, step = 0.01)
                                        ),
                                    mainPanel(width = 9,
                                              plotlyOutput("GSEA_Barplots", height = "550px"),
                                              dataTableOutput("selected_pathway_table")
                                              )
                                              )
                                            ),
                                  tabPanel("Table",
                                           sidebarLayout(
                                             sidebarPanel(width = 3,
                                               sliderInput('p_adj_2', 'Filter table based on p-adj', min = 0, max = 1, value = 0.5, step = 0.01),
                                               radioButtons('nes_selection', 'Select NES Pathways', choices = c('All', 'Positive', 'Negative'), selected = 'All'),
                                               downloadButton("download_filtered_data", "Download Data")
                                             ),
                                             mainPanel(width = 9,
                                               dataTableOutput("Filtered_GSEA_Table")
                                             )
                                           )
                                  ),
                                  tabPanel("Scatter-Plot",
                                           sidebarLayout(
                                             sidebarPanel(width = 3,
                                               sliderInput('p_adj_3', 'p-adj filter for NES plot', min = 0, max = 1, value = 0.5, step = 0.01),
                                             ),
                                             mainPanel(width = 9,
                                               plotOutput("NES_scatter_plot", height = "550px")
                                             )
                                           )
                                  )
                                )
                              )
                            )
               ),
      tabPanel("De-Mystify your data",
               sidebarLayout(
                 sidebarPanel(width = 3,
                   fileInput("upload_sample_file", 'Load your sample file', accept = c(".csv", ".tsv")),
                   fileInput("upload_counts_matrix", 'Load your counts matrix file', accept = c(".csv", ".tsv")),
                   actionButton(inputId = 'Normalize', label = 'Run DeSeq2', class = "btn-block", style = "color: black; background-color: darkorange;")
                 ),
                 mainPanel(width = 9,
                   tabsetPanel(
                     tabPanel('Normalized Counts data',
                              helpText("Use this tab to get your counts data normalized."),
                              tableOutput("normalized_counts_data"),
                              downloadButton("download_normalized_data", "Download Normalized Data")),
                     tabPanel('Differential Expression Results',
                              helpText("Use this tab to get the Differential Expression Results."),
                              tableOutput("DE_results"),
                              downloadButton("download_de_data", "Download DE Results")),
                     tabPanel('GSEA Results',
                              sidebarLayout(
                                sidebarPanel(width = 5,
                                  helpText("Use this tab to get the Gene Set Enrichment Results from your DE Results"),
                                  fileInput("Pathway file", "Provide the gene sets to compare pathways"),
                                  helpText("You can download the gene sets from GSEA site.Eg: c2.cp.v2023.Hs.symbols.gmt"),
                                  radioButtons("Ranking Metric", "Select Ranking Metric", choices = c("p-adj", "Base Mean", "log fold change"), selected = "log fold change"),
                                  radioButtons("min_size", "Select Minimum Reliability Size", choices = c("5", "10", "15"), selected = "15"),
                                  radioButtons("max_size", "Select Maximum Reliability Size", choices = c("50","100", "500"), selected = "100"),
                                  radioButtons("file_format", "Choose File Format:", choices = c("CSV", "TSV"), selected = "TSV"),
                                  downloadButton("download_fgsea_results","Download GSEA Results")
                                ),
                                mainPanel(width = 7,
                                          tableOutput("GSEA_Results_table")

                                )
                              )
                   )
                 )
               )
               )
  )
  ),class = "main-panel"  # Add the main-panel class
  )
)
###############################################################################
                    #SERVER FUNCTION#
###############################################################################
server <- function(input, output, session) {
  options(shiny.maxRequestSize = 20 * 1024^2)

  # Reactive function to read data when a file is uploaded
  cleaned_data <- reactive({
    req(input$upload_samples)
    data_filter <- read.csv(input$upload_samples$datapath, sep = ',', header = TRUE, stringsAsFactors = FALSE)
    data_filter <- as.data.frame(data_filter)
    clean_metadata(data_filter)
  })

  observeEvent(input$upload_samples, {
    req(cleaned_data())
    metadata <- cleaned_data()

    # Create an HTML string for the summary info
    summary_info <- paste(
      "Number of rows:", nrow(metadata),
      "\nNumber of columns:", ncol(metadata)
    )

    # Update the summary_info output
    output$summary_info <- renderPrint({
      cat(summary_info)
    })

    # Update the summary_output table
    output$summary_output <- renderTable({
      custom_summary(metadata)
    })
  })

  # Display cleaned data in the "Table" tab
  output$table_tab1 <- renderDataTable({
    req(cleaned_data())
    cleaned_data()
  })

  # Render the grouped bar plot in the "Plot" tab
  output$plot_tab1 <- renderPlot({
    req(cleaned_data())
    create_grouped_bar_plot(cleaned_data())
  })

################################################################################
  counts_data <- reactive({
    req(input$upload_counts)

  normalized_counts <- read.csv(input$upload_counts$datapath, sep = ',', header = TRUE, row.names = 1)
  variance_counts <- filter_zero_var_genes(normalized_counts)
  # Filter based on sliders
  variance_threshold <- quantile(apply(variance_counts, 1, var), probs = input$slider_v / 100)
  non_zero_threshold <- input$slider_nz

  filtered_genes <- rownames(variance_counts)[apply(variance_counts, 1, function(x) var(x) > variance_threshold && sum(x > 0) > non_zero_threshold)]

  # Create a summary data frame
  summary_data <- data.frame(
    "Total Genes" = nrow(normalized_counts),
    "Number of Samples" = ncol(normalized_counts),
    "Genes_Passing_Filter" = length(filtered_genes),
    "Genes_Not_Passing_Filter" = nrow(normalized_counts) - length(filtered_genes),
    "Percent Passing Filter" = sprintf("%.2f%%", 100 * length(filtered_genes) / nrow(normalized_counts)),
    "Percent Not Passing Filter" = sprintf("%.2f%%", 100 * (1 - length(filtered_genes) / nrow(normalized_counts)))
  )

  names(summary_data) <- c(
    "Total Genes\n",
    "Number of Samples\n",
    "Genes Passing Filter\n",
    "Genes Not Passing Filter\n",
    "Percent Passing Filter\n",
    "Percent Not Passing Filter\n"
  )

  summary_transpose <- as.data.frame(t(summary_data))
  summary_transpose$Attributes <- colnames(summary_data)
  colnames(summary_transpose) <- c("Counts", "Attributes")
  # Reorder the columns
  summary_transpose <- select(summary_transpose, Attributes, Counts)

  scatter_data <- data.frame(
    Gene = rownames(variance_counts),
    Median_Count = apply(variance_counts, 1, median),
    Variance = apply(variance_counts, 1, var),
    Num_Zeros = apply(variance_counts, 1, function(x) sum(x == 0)),
    Filtered = ifelse(rownames(variance_counts) %in% filtered_genes, "Passing Filter", "Not Passing Filter")
  )

  return(list(filtered_genes = filtered_genes, variance_counts = variance_counts, summary_transpose = summary_transpose , summary_data = summary_data, scatter_data = scatter_data, normalized_counts = normalized_counts))
})

# Display summary table
output$summary_table <- renderTable({
  counts_data()$summary_transpose
  })

output$scatter_plots <- renderPlot({
  # Create scatter plots with color differentiation
  p1 <- ggplot(counts_data()$scatter_data, aes(y = Variance, x = log(Median_Count), color = Filtered)) +
    geom_point() + geom_smooth(method='gam',se=FALSE, color='red')+ scale_y_continuous(trans = 'log10') +
    scale_color_manual(values = c("Passing Filter" = "black", "Not Passing Filter" = "lightgrey")) +
    labs(y = "Log(Variance)", x = "Log(Median Count)") + theme_bw() +
    ggtitle("Median Count vs Variance") + theme(plot.title = element_text(size = 18, hjust = 0.5, face = "bold"))

  p2 <- ggplot(counts_data()$scatter_data, aes(y = Num_Zeros, x = log(Median_Count), color = Filtered)) +
    geom_point() + geom_smooth(method='gam',se=FALSE, color='red') +
    scale_color_manual(values = c("Passing Filter" = "black", "Not Passing Filter" = "lightgrey")) +
    labs(y = "Number of Zeros", x = "Log(Median Count)") + theme_bw() +
    ggtitle("Median Count vs Number of Zeros") + theme(plot.title = element_text(size = 18, hjust = 0.5, face = "bold"))

  # Display plots side by side
  grid.arrange(p1, p2, ncol = 2)
  })

output$heatmap <- renderPlot({
  heatmap_data <- counts_data()$normalized_counts[rownames(counts_data()$normalized_counts) %in% counts_data()$filtered_genes, ]

  generate_heatmap(heatmap_data, title = "Clustered Heatmap of Counts", log_transform = TRUE)
  })
################################################################################

#PCA Plots

pca_result <- reactive({
  req(counts_data()$normalized_counts)  # Ensure normalized_counts is available

  pca_input <- counts_data()$normalized_counts[rownames(counts_data()$normalized_counts) %in% counts_data()$filtered_genes, ]

  transposed_data <- t(pca_input)
  perform_pca(transposed_data, num_pcs = max_pcs(transposed_data))
  result <- perform_pca(transposed_data, num_pcs = max_pcs(transposed_data))
  print(result$variance_explained)  # Add this line to print the pca_data for debugging
  return(result)
})

# Render PCA plot
output$pca_plot <- renderPlot({
  pca_data_plot <- pca_result()$pca_data
  pca_variance_exp <- pca_result()$variance_explained
  plot_pca(pca_data_plot,pca_variance_exp, pc_x = input$pc_x, pc_y = input$pc_y)
})
################################################################################

#Differential Expression Analysis

load_data <- reactive({
  req(input$upload_de_result)
  dataf <- read.csv(input$upload_de_result$datapath)
  return(dataf)
})

output$DE_data <- renderDataTable({
  req(input$upload_de_result)
  load_data()
}
)


volcano_plot <- function(dataf, x_name, y_name, slider, color1, color2) {
  ggplot(data = dataf, aes(x = !!sym(x_name), y = -log10(!!sym(y_name)), color = !!sym(y_name) < 10^slider)) +
    geom_point() +
    labs(x = "log2(fold-change)", y = "-log10(p adjusted)") +
    scale_color_manual(
      name = paste(ensym(y_name), "< 1 x 10^", slider),
      values = c("FALSE" = color1, "TRUE" = color2)) + theme_bw() + theme(legend.position = 'bottom')
}

draw_table <- function(dataf, slider_p) {
  print(str(dataf))
  print(head(dataf))

  dataf <- na.omit(dataf)

  # Filter based on the slider value
  filtered_data_de <- dataf[dataf$padj < 10^slider_p, ]

  # Format the p-value and p-adjusted columns to display more digits
  filtered_data_de$pvalue <- formatC(filtered_data_de$pvalue, format = "e", digits = 5)
  filtered_data_de$padj <- formatC(filtered_data_de$padj, format = "e", digits = 5)
  names(filtered_data_de)[1] <- 'gene'

  return(filtered_data_de)
}

plotClicked <- reactiveVal(FALSE)

selectedColor1 <- reactiveVal('#22577A')
selectedColor2 <- reactiveVal('#FFCF56')
selectedSlider <- reactiveVal(-150)
selectedXName <- reactiveVal('log2FoldChange')
selectedYName <- reactiveVal('padj')

tableChanged <- reactiveVal(FALSE)

observe({
  req(input$upload)
  plotClicked(TRUE)
})

observeEvent(input$Plot_v, {
  plotClicked(TRUE)
})

observe({
  req(input$upload)
  tableChanged(TRUE)
})

observeEvent(input$Plot_v, {
  tableChanged(TRUE)
})

# Reset the reactive value when the slider changes
observeEvent(c(input$slider_p, input$color1, input$color2, input$x_name, input$y_name), {
  selectedSlider(input$slider_p)
  selectedColor1(input$color1)
  selectedColor2(input$color2)
  selectedXName(input$x_name)
  selectedYName(input$y_name)
  plotClicked(FALSE)
})

observeEvent(c(input$color1, input$color2, input$x_name, input$y_name), {
  selectedSlider(input$slider_p)
  selectedColor1(input$color1)
  selectedColor2(input$color2)
  selectedXName(input$x_name)
  selectedYName(input$y_name)
  tableChanged(TRUE)
})

observeEvent(c(input$slider_p), {
  tableChanged(FALSE)
})


#' These outputs aren't really functions, so they don't get a full skeleton,
#' but use the renderPlot() and renderTabel() functions to return() a plot
#' or table object, and those will be displayed in your application.
output$volcano <- renderPlot({
  req(plotClicked())
  volcano_plot(load_data(), selectedXName(), selectedYName(), selectedSlider(), selectedColor1(), selectedColor2())
},width = 'auto', height = 650)


# Same here, just return the table as you want to see it in the web page
output$table_output <- DT::renderDataTable({
  req(tableChanged())
  draw_table(load_data(), selectedSlider())
}
)
################################################################################

#Gene Set Enrichment Analysis
# Gene Set Enrichment Analysis
uploaded_gsea <- eventReactive(input$Upload_Results,{
  req(input$upload_gsea)
  read.table(input$upload_gsea$datapath, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
})

# Reactive expression for filtered results based on adjusted p-value
filtered_gsea <- reactive({
  req(uploaded_gsea())
  uploaded_gsea_filtered <- uploaded_gsea()
  uploaded_gsea_filtered[uploaded_gsea_filtered$padj <= input$p_adj, ]
})

output$selected_pathway_table <- renderDataTable({
  req(filtered_gsea())
  filtered_gsea()
})

# Render the barplot for the "Barplots" tab using plot_ly for interactivity
selected_pathway <- reactiveVal(NULL)

# Render the barplot for the "Barplots" tab using plot_ly for interactivity
output$GSEA_Barplots <- renderPlotly({
  req(filtered_gsea())

  # Check if filtered_gsea() is not NULL and contains data
  if (!is.null(filtered_gsea()) && nrow(filtered_gsea()) > 0) {
    # Filter based on adjusted p-value
    filtered_gsea_display <- filtered_gsea()[filtered_gsea()$padj <= input$p_adj_2, ]

    # Order levels of the "pathway" factor based on NES
    filtered_gsea_display$pathway <- factor(filtered_gsea_display$pathway,
                                            levels = filtered_gsea_display$pathway[order(filtered_gsea_display$NES)])

    # Create barplot with positive and negative colors
    p <- plot_ly(data = filtered_gsea_display, x = ~NES, y = ~pathway, type = "bar",
                 color = ~factor(sign(NES), levels = c(-1, 1), labels = c("Negative", "Positive")),
                 colors = c("darkcyan", "indianred")) %>%
      layout(title = "Pathway Analysis", xaxis = list(title = "Normalized Enrichment Score"), yaxis = list(title = "Pathways"))

    # Add click event handler
    p <- htmlwidgets::onRender(p, '
      function(el, x) {
        el.on("plotly_click", function(event) {
          var pathway = event.points[0].y;
          Shiny.setInputValue("selected_pathway", pathway);
        });
      }
    ')

    p
  }
})

# Observe the selected pathway and update the reactive value
observeEvent(input$selected_pathway, {
  selected_pathway(input$selected_pathway)
})

# Render the filtered and sortable data table for the "Table" tab
output$selected_pathway_table <- DT::renderDataTable({
  req(filtered_gsea(), selected_pathway())

  # Filter the table based on the selected pathway
  filtered_table <- filtered_gsea()
  if (!is.null(selected_pathway())) {
    filtered_table <- filtered_table[filtered_table$pathway == selected_pathway(), ]
  }

  DT::datatable(filtered_table, options = list(order = list(2, 'desc')))
})

# Render the filtered and sortable data table for the "Table" tab
output$Filtered_GSEA_Table <- DT::renderDataTable({
  req(uploaded_gsea())

  # Filter based on adjusted p-value
  uploaded_gsea_filtered <- uploaded_gsea()
  uploaded_gsea_filtered <- uploaded_gsea_filtered[uploaded_gsea_filtered$padj <= input$p_adj_2, ]

  # Filter based on NES selection (positive, negative, or all)
  if (input$nes_selection == 'Positive') {
    uploaded_gsea_filtered <- uploaded_gsea_filtered[uploaded_gsea_filtered$NES > 0, ]
  } else if (input$nes_selection == 'Negative') {
    uploaded_gsea_filtered <- uploaded_gsea_filtered[uploaded_gsea_filtered$NES < 0, ]
  }

  DT::datatable(uploaded_gsea_filtered, options = list(order = list(2, 'desc')))
})

# Download filtered data
output$download_filtered_data <- downloadHandler(
  filename = function() {
    paste("filtered_data_", Sys.Date(), ".csv", sep = "\t")
  },
  content = function(file) {
    write.csv(filtered_table(), file, row.names = FALSE)
  }
)

# NES scatter data
NES_scatter_data <- reactive({
  req(uploaded_gsea())
  uploaded_gsea_filtered <- uploaded_gsea()

  # Filter based on adjusted p-value
  uploaded_gsea_filtered
})

# Render the Scatter Plot
output$NES_scatter_plot <- renderPlot({
  req(NES_scatter_data())

  # Scatter plot with grey color for gene sets below threshold
  ggplot(NES_scatter_data(), aes(x = NES, y = -log10(padj))) +
    geom_point(aes(color = ifelse(padj <= input$p_adj_3, "Above Threshold", "Below Threshold"))) +
    scale_color_manual(values = c("Below Threshold" = "grey", "Above Threshold" = "#FFCF56")) +
    labs(x = "NES", y = "-log10(Adjusted p-value)", color = "Threshold") +
    theme_minimal() + ggtitle("Pathway Enrichment Landscape: Exploring NES vs. Significance") +
    theme(plot.title = element_text(size = 18, hjust = 0.5, face = "bold"))
})

################################################################################
#Demystify tab

deSeqResults <- reactiveVal(NULL)

sample_data <- reactive({
  req(input$upload_sample_file)
  filter_data <- read.csv(input$upload_sample_file$datapath, sep = ',', header = TRUE, row.names = 1, stringsAsFactors = FALSE)
  filter_data <- as.data.frame(filter_data)
  clean_metadata(filter_data)
})

counts_metadata <- reactive({
  req(input$upload_counts_matrix)
  counts <- read.csv(input$upload_counts_matrix$datapath, sep = ',', header = TRUE, row.names = 1)
  counts <- as.data.frame(counts)
})

runDeSeq2 <- function(sample_data, counts_metadata) {
  dds <- DESeqDataSetFromMatrix(countData = counts_metadata, colData = sample_data, design = ~ Treatment + Timepoint)
  dds <- DESeq(dds)
  results <- results(dds)
  return(list(dds_normalized = dds, de_results = results))
}

observeEvent(input$Normalize, {
  # Run DeSeq2 when the button is clicked
  deSeqResultsVal <- runDeSeq2(sample_data(), counts_metadata())

  # Update the reactive value with the results
  deSeqResults(deSeqResultsVal)
})

# Download handler for normalized data
output$download_normalized_data <- downloadHandler(
  filename = function() {
    "normalized_data.csv"
  },
  content = function(file) {
    # Access the reactive value to get deSeqResults
    write.csv(assay(deSeqResults()$dds_normalized), file, row.names = FALSE)
  }
)

# Download handler for DE results
output$download_de_data <- downloadHandler(
  filename = function() {
    "DE_results.csv"
  },
  content = function(file) {
    # Access the reactive value to get deSeqResults
    write.csv(deSeqResults()$de_results, file, row.names = FALSE)
  }
)
}


shinyApp(ui = ui, server = server)
