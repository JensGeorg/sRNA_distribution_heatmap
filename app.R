# Required Libraries
library(shiny)
library(shinydashboard)
library(ComplexHeatmap)
library(grid)
library(colourpicker)
library(circlize)

# Load the data (update the path as necessary)
load("heatdata.Rdata")
tax_table <- heatdata[[1]]
coor_all <- heatdata[[2]]
anno2 <- heatdata[[3]]

# Default organism file
default_organism_file <- "salmonella_all.txt"

# UI for the Shiny app using shinydashboard
ui <- dashboardPage(
  dashboardHeader(title = "sRNA Heatmap"),
  
  dashboardSidebar(
    width = 380,
    div(style = "height: 90vh; overflow-y: auto; padding: 10px;",
        
        # Help Section
        box(
          title = "Information", status = "info", solidHeader = TRUE,
          width = 12, collapsible = TRUE, collapsed = TRUE,
          
          h4("Application Overview"),
          p("This application generates interactive heatmaps showing the phylogentic distribution of sRNA homologs. 
          It processes genomic data to visualize the presence and abundance of sRNA homologs in various taxonomic groups.
            sRNA homologs are detected using the GLASSgo tool. Taxonomic information was linked using a custom R script (github) and TaxonKit (ref)"),
          
          
        ),
        
        box(
          title = "sRNA selection", 
          status = "primary", 
          solidHeader = TRUE,
          width = 12, 
          collapsible = TRUE,

          textInput(
            "include", 
            "Included sRNA clusters (comma-separated):", 
            value = paste(names(coor_all), collapse = ",")
          ),
          helpText("Leave empty to include all clusters."),
          hr(), # Adds a visual separator line
          #--- 4. Identity Threshold ---#
          numericInput(
            "identity", 
            "Identity threshold (%):", 
            min = 0, 
            max = 100, 
            value = 57, 
            step = 1
          ),
          helpText("Minimum sequence identity to reference homolog for inclusion."),
          
          hr(), # Adds a visual separator line
          
          div(
            tags$label("Use organism/group file:", class = "control-label", style = "display: block; font-weight: bold; margin-bottom: 5px;"),
            checkboxInput("use_cluster_file", "Use sRNA cluster file:", value = FALSE),
            style = "margin-bottom: 15px;"
          ),
          
            # Item 3: The Help Text
            
          #),
          #--- 2. Conditional File Upload ---#
          conditionalPanel(
            condition = "input.use_cluster_file",
            helpText(
              # style = "margin-left: 20px; margin-bottom: -10;",
              HTML("Textfile with one line per sRNA e.g.:<br>cluster_3<br>RyhB")
            ),
            fileInput(
              "cluster_file", 
              "Upload cluster file:",
              accept = ".txt",
              buttonLabel = "Browse...",
              placeholder = "No file selected"
            )
            
          )
        ),
        
        # Organism Selection Section
        box(
          p(HTML("Selections from <b> taxonomic level </b>,<b> hierarchical taxonomic Filters </b> and <b> custom organism file </b> are additive"), 
            style = "font-size: 12px; color: #777;"),
          hr(), # Adds a visual separator line
          title = "Organism Selection", status = "success", solidHeader = TRUE,
          width = 12, collapsible = TRUE,
          
          div(
            tags$label("Use organism/group file:", class = "control-label", style = "display: block; font-weight: bold; margin-bottom: 5px;"),
            checkboxInput("use_organism_file", NULL, value = FALSE),
            style = "margin-bottom: 15px;"
          ),
          
          conditionalPanel(
            condition = "input.use_organism_file",
            div(
              tags$label("Upload organism file:", class = "control-label", style = "display: block; font-weight: bold; margin-bottom: 5px;"),
              fileInput("organism_file", NULL,
                        accept = c(".txt"),
                        buttonLabel = "Browse...",
                        placeholder = "No file selected"),
              helpText(HTML("One organism ID or taxonomic name per line, e.g.:<br>
                  Escherichia coli<br>
                  Salmonella enterica")),
              style = "margin-bottom: 15px;"
            )
          ),
          hr(), # Adds a visual separator line
          div(
            tags$label("Show taxonomic filters:", class = "control-label", style = "display: block; font-weight: bold; margin-bottom: 5px;"),
            checkboxInput("expand_organism_options", NULL, value = TRUE),
            style = "margin-bottom: 15px;"
          ),
          
          conditionalPanel(
            condition = "input.expand_organism_options",
            hr(), # Adds a visual separator line
            div(
              tags$label("Taxonomic level:", class = "control-label", style = "display: block; font-weight: bold; margin-bottom: 5px;"),
              selectInput("tax_level", NULL,
                          choices = c("None", colnames(tax_table)),
                          selected = "family"),
              style = "margin-bottom: 15px;"
            ),
            p("sRNA numbers are cumulated for the selected taxonomic level", 
              style = "font-size: 12px; color: #777;"),
            hr(), # Adds a visual separator line
            h5("Hierarchical Taxonomic Filters:", style = "color: #555; font-weight: bold;"),
            p("Enter comma-separated values. When selecting a taxnomic group all relevant members of the next lower taxonomic level belonging to this group are displayed.
              (e.g. when selecting the genus 'Salmonella' all Species belonging to this genus are displayed)", 
              style = "font-size: 12px; color: #777;"),
            
            div(
              tags$label("Phylum:", class = "control-label", style = "display: block; font-weight: bold; margin-bottom: 5px;"),
              textInput("individual_phylum", NULL, placeholder = "e.g., Proteobacteria, Bacillota, Firmicutes"),
              style = "margin-bottom: 10px;"
            ),
            
            div(
              tags$label("Class:", class = "control-label", style = "display: block; font-weight: bold; margin-bottom: 5px;"),
              textInput("individual_class", NULL, placeholder = "e.g., Gammaproteobacteria, Bacilli"),
              style = "margin-bottom: 10px;"
            ),
            
            div(
              tags$label("Order:", class = "control-label", style = "display: block; font-weight: bold; margin-bottom: 5px;"),
              textInput("individual_order", NULL, placeholder = "e.g., Enterobacterales, Bacillales"),
              style = "margin-bottom: 10px;"
            ),
            
            div(
              tags$label("Family:", class = "control-label", style = "display: block; font-weight: bold; margin-bottom: 5px;"),
              textInput("individual_family", NULL,  placeholder = "e.g., Enterobacteriaceae, Staphylococcaceae"),
              style = "margin-bottom: 10px;"
            ),
            
            div(
              tags$label("Genus:", class = "control-label", style = "display: block; font-weight: bold; margin-bottom: 5px;"),
              textInput("individual_genus", NULL, placeholder = "e.g., Salmonella, Escherichia"),
              style = "margin-bottom: 10px;"
            ),
            
            div(
              tags$label("Species:", class = "control-label", style = "display: block; font-weight: bold; margin-bottom: 5px;"),
              textInput("individual_species", NULL,  
                        placeholder = "e.g., Salmonella enterica"),
              style = "margin-bottom: 10px;"
            )
          )
        ),
        
        # Heatmap Options Section
        box(
          title = "Heatmap Display Options", status = "warning", solidHeader = TRUE,
          width = 12, collapsible = TRUE,
          
          fluidRow(
            column(6,
                   div(
                     tags$label("Show row names:", class = "control-label", style = "display: block; font-weight: bold; margin-bottom: 5px;"),
                     checkboxInput("show_row", NULL, value = FALSE),
                     style = "margin-bottom: 15px;"
                   ),
                   div(
                     tags$label("Cluster rows:", class = "control-label", style = "display: block; font-weight: bold; margin-bottom: 5px;"),
                     checkboxInput("cluster_row", NULL, value = TRUE),
                     style = "margin-bottom: 15px;"
                   ),
                   div(
                     tags$label("Presence/Absence mode:", class = "control-label", style = "display: block; font-weight: bold; margin-bottom: 5px;"),
                     checkboxInput("presence_absence", NULL, value = FALSE),
                     helpText("Convert counts to binary (0/1)"),
                     style = "margin-bottom: 15px;"
                   )
            ),
            column(6,
                   div(
                     tags$label("Show column names:", class = "control-label", style = "display: block; font-weight: bold; margin-bottom: 5px;"),
                     checkboxInput("show_col", NULL, value = TRUE),
                     style = "margin-bottom: 15px;"
                   ),
                   div(
                     tags$label("Cluster columns:", class = "control-label", style = "display: block; font-weight: bold; margin-bottom: 5px;"),
                     checkboxInput("cluster_col", NULL, value = TRUE),
                     style = "margin-bottom: 15px;"
                   )
            )
          ),
          
          div(
            tags$label("Minimum column sum for inclusion:", class = "control-label", style = "display: block; font-weight: bold; margin-bottom: 5px;"),
            numericInput("min_col_sum", NULL, value = 1, min = 0, step = 1),
            helpText("Exclude organisms with total counts below this threshold"),
            style = "margin-bottom: 15px;"
          ),
          
          h5("Label Sizes:", style = "color: #555; font-weight: bold;"),
          fluidRow(
            column(6, 
                   div(
                     tags$label("Row label size:", class = "control-label", style = "display: block; font-weight: bold; margin-bottom: 5px;"),
                     numericInput("row_label_size", NULL, value = 10, min = 6, max = 20, step = 1),
                     style = "margin-bottom: 10px;"
                   )
            ),
            column(6, 
                   div(
                     tags$label("Column label size:", class = "control-label", style = "display: block; font-weight: bold; margin-bottom: 5px;"),
                     numericInput("col_label_size", NULL, value = 10, min = 6, max = 20, step = 1),
                     style = "margin-bottom: 10px;"
                   )
            )
          )
        ),
        
        # Color Scale Section
        box(
          title = "Color Scale Configuration", status = "danger", solidHeader = TRUE,
          width = 12, collapsible = TRUE, collapsed = FALSE,
          
          helpText("Define 5 split points and corresponding colors for the heatmap scale. Values will be interpolated between these points."),
          
          div(style = "background-color: #f9f9f9; padding: 10px; border-radius: 5px; margin-bottom: 15px;",
              fluidRow(
                column(6, h5("Split Points", style = "text-align: center; color: #333;")),
                column(6, h5("Colors", style = "text-align: center; color: #333;"))
              )
          ),
          
          fluidRow(
            column(6, 
                   div(
                     tags$label("Point 1:", class = "control-label", style = "display: block; font-weight: bold; margin-bottom: 5px;"),
                     numericInput("split1", NULL, value = 0, step = 1),
                     style = "margin-bottom: 10px;"
                   )
            ),
            column(6, 
                   div(
                     tags$label("Color 1:", class = "control-label", style = "display: block; font-weight: bold; margin-bottom: 5px;"),
                     colourInput("color1", NULL, value = "#FFFFFF", 
                                 palette = "square", showColour = "both"),
                                      style = "margin-bottom: 10px;"
                   )
            )
          ),
          fluidRow(
            column(6, 
                   div(
                     tags$label("Point 2:", class = "control-label", style = "display: block; font-weight: bold; margin-bottom: 5px;"),
                     numericInput("split2", NULL, value = 1, step = 1),
                     style = "margin-bottom: 10px;"
                   )
            ),
            column(6, 
                   div(
                     tags$label("Color 2:", class = "control-label", style = "display: block; font-weight: bold; margin-bottom: 5px;"),
                     colourInput("color2", NULL, value = "#FFA500"),
                     style = "margin-bottom: 10px;"
                   )
            )
          ),
          fluidRow(
            column(6, 
                   div(
                     tags$label("Point 3:", class = "control-label", style = "display: block; font-weight: bold; margin-bottom: 5px;"),
                     numericInput("split3", NULL, value = 5, step = 1),
                     style = "margin-bottom: 10px;"
                   )
            ),
            column(6, 
                   div(
                     tags$label("Color 3:", class = "control-label", style = "display: block; font-weight: bold; margin-bottom: 5px;"),
                     colourInput("color3", NULL, value = "#FFFF00"),
                     style = "margin-bottom: 10px;"
                   )
            )
          ),
          fluidRow(
            column(6, 
                   div(
                     tags$label("Point 4:", class = "control-label", style = "display: block; font-weight: bold; margin-bottom: 5px;"),
                     numericInput("split4", NULL, value = 25, step = 1),
                     style = "margin-bottom: 10px;"
                   )
            ),
            column(6, 
                   div(
                     tags$label("Color 4:", class = "control-label", style = "display: block; font-weight: bold; margin-bottom: 5px;"),
                     colourInput("color4", NULL, value = "#008000"),
                     style = "margin-bottom: 10px;"
                   )
            )
          ),
          fluidRow(
            column(6, 
                   div(
                     tags$label("Point 5:", class = "control-label", style = "display: block; font-weight: bold; margin-bottom: 5px;"),
                     numericInput("split5", NULL, value = 50, step = 1),
                     style = "margin-bottom: 10px;"
                   )
            ),
            column(6, 
                   div(
                     tags$label("Color 5:", class = "control-label", style = "display: block; font-weight: bold; margin-bottom: 5px;"),
                     colourInput("color5", NULL, value = "#0000FF"),
                     style = "margin-bottom: 10px;"
                   )
            )
          )
        ),
        
        
        
        # Action and Export Section
        box(
          title = "Generate & Export", status = "info", solidHeader = TRUE,
          width = 12,
          
          div(
            actionButton("generate", "Generate Heatmap",
                         class = "btn-primary btn-lg",
                         style = "width: 80%; margin-bottom: 20px; padding: 12px;"),
            style = "text-align: center;"
          ),
          
          h5("PDF Export Settings:", style = "color: #555; font-weight: bold; margin-bottom: 15px;"),
          fluidRow(
            column(6, 
                   div(
                     tags$label("Width (inches):", class = "control-label", style = "display: block; font-weight: bold; margin-bottom: 5px;"),
                     numericInput("pdf_width", NULL, value = 10, min = 5, max = 20, step = 1),
                     style = "margin-bottom: 10px;"
                   )
            ),
            column(6, 
                   div(
                     tags$label("Height (inches):", class = "control-label", style = "display: block; font-weight: bold; margin-bottom: 5px;"),
                     numericInput("pdf_height", NULL, value = 8, min = 5, max = 20, step = 1),
                     style = "margin-bottom: 10px;"
                   )
            )
          ),
          
          div(
            downloadButton("download_pdf", "Download PDF",
                           class = "btn-success",
                           style = "width: 100%; margin-bottom: 10px; padding: 8px;"),
            style = "text-align: center;"
          ),
          
          div(
            downloadButton("download_matrix", "Download Data Matrix (csv)",
                           class = "btn-warning",
                           style = "width: 100%; padding: 8px;"),
            style = "text-align: center;"
          )
        )
    )
  ),
  
  dashboardBody(
    tags$head(
      tags$style(HTML("
        .box p, .box li {
         color: #000 !important;
        }
        .content-wrapper, .right-side {
          background-color: #f4f4f4;
        }
        .box {
          margin-bottom: 10px;
        }
        .control-label {
          font-weight: bold !important;
          color: #333 !important;
          font-size: 13px !important;
        }
        .form-group {
          margin-bottom: 15px;
        }
        h5 {
          color: #555;
          font-weight: bold;
          margin-top: 15px;
        }
        .help-block {
          font-size: 11px;
          color: #777;
          margin-top: 2px;
        }
        .btn {
          border-radius: 4px;
        }
        .box-title {
          font-size: 16px;
          font-weight: bold;
        }
        .shiny-input-container {
          margin-bottom: 5px;
        }
      "))
    ),
    
    fluidRow(
      # Main heatmap display
      box(
        title = "Heatmap Visualization", status = "primary", solidHeader = TRUE,
        width = 12, height = "800px",
        
        div(style = "text-align: center; padding: 20px;",
            conditionalPanel(
              condition = "input.generate == 0",
              div(style = "margin-top: 200px; color: #999;",
                  h3("Click 'Generate Heatmap' to display visualization"),
                  br(),
                  icon("chart-bar", class = "fa-3x"),
                  br(), br(),
                  p("Configure your data filters and display options in the sidebar, then generate your heatmap.",
                    style = "color: #666; max-width: 400px; margin: 0 auto;")
              )
            ),
            conditionalPanel(
              condition = "input.generate > 0",
              plotOutput("heatmapPlot", height = "700px")
            )
        )
      )
    ),
    
    fluidRow(
      # Data summary
      box(
        title = "📊 Data Summary", status = "info", solidHeader = TRUE,
        width = 6, collapsible = TRUE,
        
        verbatimTextOutput("data_summary")
      ),
      
      # Filter status
      box(
        title = "⚙️ Current Filter Settings", status = "success", solidHeader = TRUE,
        width = 6, collapsible = TRUE,
        
        verbatimTextOutput("filter_status")
      )
    )
  )
)

# Server logic (same as before)
server <- function(input, output, session) {
  
  filtered_data <- eventReactive(input$generate, {
    tmp_orgs <- c()
    
    # Load clusters from file if provided and selected
    if (input$use_cluster_file && !is.null(input$cluster_file)) {
      clusters <- readLines(input$cluster_file$datapath)
      inc <- clusters
    } else {
      if (length(input$include) == 0 || input$include == "") {
        inc <- names(coor_all)
      } else {
        inc <- strsplit(input$include, ",")[[1]]
        inc <- trimws(inc)  # Remove whitespace
      }
    }
    
    # Load organisms from file if provided and selected
    if (input$use_organism_file && !is.null(input$organism_file)) {
      tmp_orgs <- readLines(input$organism_file$datapath)
    } else if (input$use_organism_file && file.exists(default_organism_file)) {
      tmp_orgs <- readLines(default_organism_file)
    } else if (isTRUE(input$expand_organism_options)) {
      if (input$tax_level != "None") {
        tmp <- tax_table[, input$tax_level]
        tmp_orgs <- c(tmp_orgs, unique(tmp))
      }
      
      if (input$individual_phylum != "") {
        phylum <- strsplit(input$individual_phylum, ",")[[1]]
        phylum <- trimws(phylum)
        tmp <- which(tax_table[, "phylum"] %in% phylum)
        tmp <- tax_table[tmp, "class"]
        tmp_orgs <- c(tmp_orgs, unique(tmp))
      }
      
      if (input$individual_class != "") {
        class <- strsplit(input$individual_class, ",")[[1]]
        class <- trimws(class)
        tmp <- which(tax_table[, "class"] %in% class)
        tmp <- tax_table[tmp, "order"]
        tmp_orgs <- c(tmp_orgs, unique(tmp))
      }
      
      if (input$individual_order != "") {
        order <- strsplit(input$individual_order, ",")[[1]]
        order <- trimws(order)
        tmp <- which(tax_table[, "order"] %in% order)
        tmp <- tax_table[tmp, "family"]
        tmp_orgs <- c(tmp_orgs, unique(tmp))
      }
      
      if (input$individual_family != "") {
        families <- strsplit(input$individual_family, ",")[[1]]
        families <- trimws(families)
        tmp <- which(tax_table[, "family"] %in% families)
        tmp <- tax_table[tmp, "genus"]
        tmp_orgs <- c(tmp_orgs, unique(tmp))
      }
      
      if (input$individual_genus != "") {
        genus <- strsplit(input$individual_genus, ",")[[1]]
        genus <- trimws(genus)
        tmp <- which(tax_table[, "genus"] %in% genus)
        tmp <- tax_table[tmp, "species"]
        tmp_orgs <- c(tmp_orgs, unique(tmp))
      }
      
      if (input$individual_species != "") {
        species <- strsplit(input$individual_species, ",")[[1]]
        species <- trimws(species)
        tmp <- which(tax_table[, "species"] %in% species)
        tmp <- tax_table[tmp, "subspecies"]
        tmp_orgs <- c(tmp_orgs, unique(tmp))
      }
    }
    
    tmp_orgs <- na.omit(tmp_orgs)
    tmp_orgs <- unique(tmp_orgs)
    
    if (length(tmp_orgs) == 0) {
      showNotification("No organisms match the current filters. Please adjust your settings.", 
                       type = "warning", duration = 5)
      return(NULL)
    }
    
    tmp_coor <- coor_all[inc]
    heat_data <- matrix(0, length(tmp_coor), length(tmp_orgs))
    rownames(heat_data) <- names(tmp_coor)
    colnames(heat_data) <- tmp_orgs
    
    # Build heatmap data
    for (i in 1:length(tmp_coor)) {
      tmp <- tmp_coor[[i]]
      iden <- which(as.numeric(tmp[, "identity"]) >= input$identity)
      
      if (length(iden) > 0) {
        tmp <- tmp[iden, 10:19]
        pos <- table(tmp)
        pos2 <- match(names(pos), colnames(heat_data))
        if (length(na.omit(pos2)) > 0) {
          pos3 <- which(!is.na(pos2))
          heat_data[i, na.omit(pos2)] <- heat_data[i, na.omit(pos2)] + pos[pos3]
        }
      }
    }
    
    # Convert to presence/absence if selected
    if (input$presence_absence) {
      heat_data[heat_data > 1] <- 1
    }
    
    # Exclude columns with sum below the minimum threshold
    col_sums <- colSums(heat_data)
    heat_data <- heat_data[, col_sums >= input$min_col_sum, drop = FALSE]
    
    if (ncol(heat_data) == 0) {
      showNotification("No columns remain after applying the minimum column sum filter. Please reduce the threshold.", 
                       type = "warning", duration = 5)
      return(NULL)
    }
    
    heat_data
  })
  
  # Render heatmap
  output$heatmapPlot <- renderPlot({
    req(filtered_data())
    heatmap_data <- filtered_data()
    
    if (is.null(heatmap_data)) return(NULL)
    
    # Generate custom color scale
    split_points <- c(input$split1, input$split2, input$split3, input$split4, input$split5)
    colors <- c(input$color1, input$color2, input$color3, input$color4, input$color5)
    
    color_mapping <- colorRamp2(split_points, colors)
    
    ht <- Heatmap(
      heatmap_data,
      col = color_mapping,
      show_row_names = input$show_row,
      show_column_names = input$show_col,
      cluster_rows = input$cluster_row,
      cluster_columns = input$cluster_col,
      row_names_gp = gpar(fontsize = input$row_label_size),
      column_names_gp = gpar(fontsize = input$col_label_size),
      heatmap_legend_param = list(title = "Count")
    )
    
    draw(ht)
    ht
  })
  
  # Data summary output
  output$data_summary <- renderText({
    data <- filtered_data()
    
    if (is.null(data)) {
      return("Click 'Generate Heatmap' to see data summary\n\nThis panel will show:\n• Matrix dimensions\n• Value statistics\n• Data completeness")
    }
    
    paste0(
      "Matrix dimensions: ", nrow(data), " rows × ", ncol(data), " columns\n",
      "Total data points: ", length(data), "\n",
      "Non-zero values: ", sum(data > 0), " (", round(100 * sum(data > 0) / length(data), 1), "%)\n",
      "Value range: ", min(data), " - ", max(data), "\n"
     # "Mean value: ", round(mean(data), 2), "\n",
     # "Median value: ", round(median(data), 2), "\n",
     # "Standard deviation: ", round(sd(data), 2)
    )
  })
  
  # Filter status output
  output$filter_status <- renderText({
    paste0(
      "Identity threshold: ", input$identity, "%\n",
      "Min column sum: ", input$min_col_sum, "\n",
      "Presence/Absence mode: ", ifelse(input$presence_absence, "Enabled", "Disabled"), "\n",
      "Row clustering: ", ifelse(input$cluster_row, "Enabled", "Disabled"), "\n",
      "Column clustering: ", ifelse(input$cluster_col, "Enabled", "Disabled"), "\n",
      "Show row names: ", ifelse(input$show_row, "Yes", "No"), "\n",
      "Show column names: ", ifelse(input$show_col, "Yes", "No"), "\n",
      "Using cluster file: ", ifelse(input$use_cluster_file, "Yes", "No"), "\n",
      "Using organism file: ", ifelse(input$use_organism_file, "Yes", "No")
    )
  })
  
  # Download handlers
  output$download_pdf <- downloadHandler(
    filename = function() { 
      paste("heatmap_identity_", input$identity, "_", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      req(filtered_data())
      
      pdf(file, width = input$pdf_width, height = input$pdf_height)
      
      suppressMessages({
        heatmap_data <- filtered_data()
        
        if (!is.null(heatmap_data)) {
          split_points <- c(input$split1, input$split2, input$split3, input$split4, input$split5)
          colors <- c(input$color1, input$color2, input$color3, input$color4, input$color5)
          
          color_mapping <- colorRamp2(split_points, colors)
          
          ht <- Heatmap(
            heatmap_data,
            col = color_mapping,
            show_row_names = input$show_row,
            show_column_names = input$show_col,
            cluster_rows = input$cluster_row,
            cluster_columns = input$cluster_col,
            row_names_gp = gpar(fontsize = input$row_label_size),
            column_names_gp = gpar(fontsize = input$col_label_size),
            heatmap_legend_param = list(title = "Count")
          )
          
          draw(ht)
        }
      })
      
      dev.off()
    }
  )
  
  output$download_matrix <- downloadHandler(
    filename = function() { 
      paste("heatmap_data_identity_", input$identity, "_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      req(filtered_data())
      heatmap_data <- filtered_data()
      
      if (!is.null(heatmap_data)) {
        ht <- Heatmap(
          heatmap_data,
          cluster_rows = input$cluster_row,
          cluster_columns = input$cluster_col
        )
        ht_draw <- draw(ht)
        row_order <- row_order(ht_draw)
        col_order <- column_order(ht_draw)
        
        heatmap_data <- heatmap_data[row_order, col_order, drop = FALSE]
        write.csv(heatmap_data, file)
      }
    }
  )
}

# Run the app
shinyApp(ui = ui, server = server)