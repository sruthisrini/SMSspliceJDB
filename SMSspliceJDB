library(shiny)
library(GenomicRanges)
library(ggplot2)
library(trackViewer)
library(maftools)
library(shinythemes)

# Read the files 
file1 <- read.table("genes_from_gencode.txt", sep = "\t", header = TRUE) #This file contains genes from GENCODE consortium which includes Chr_id, Start, End, Interval,Strand, Transcript and Gene name
file2 <- read.delim("variant.txt", header = TRUE) #This file contains the variant information from COSMIC database associated with 6 splice site positions (3 from acceptors and 3 from donors)

domain <- read.table("domain_new.txt", sep = "\t", header = TRUE) #This file contains information about protein domains of the genes
exp_val <- read.delim("exp_val_mutations.txt", header = TRUE) #This file contains experimentally validated mutations associated with 6 splice site positions (3 from acceptors and 3 from donors)

splice_site_colors <- c("A" = "lightblue", "C" = "lightgreen", "G" = "lightcoral", "T" = "lightsalmon")
mutations_colors <- c(
  "A->C" = "#a6cee3",   # Light Blue
  "A->G" = "#1f78b4",   # Dark Blue
  "A->T" = "#b2df8a",   # Light Green
  "C->A" = "#33a02c",   # Dark Green
  "C->G" = "#fb9a99",   # Light Coral
  "C->T" = "#e31a1c",   # Red
  "G->A" = "#fdbf6f",   # Light Orange
  "G->C" = "#ff7f00",   # Orange
  "G->T" = "#cab2d6",   # Light Purple
  "T->A" = "#6a3d9a",   # Dark Purple
  "T->C" = "#a6bddb",   # Lavender
  "T->G" = "#1A5276"    # Dark Blue
)

ui <- fluidPage(theme = shinytheme("sandstone"),

  navbarPage("SMSspliceJDB",
             tabPanel("Potential mutations affecting splice site",
                      sidebarLayout(
                        sidebarPanel(width = 3,
                          textInput("geneInput", "Enter Gene:",width = '50%'),
                          actionButton("plotButton", "Submit"),
                        ),
                        mainPanel(align="center",
                          plotOutput("lollipopPlot", width = "100%"),
                          plotOutput("Protein_lollipop"),
                          plotOutput("barPlot3"),
                          div(
                            
                            tableOutput("filteredTable"),  # Table output
                            downloadButton("downloadData", "Download Data") )
                        )
                      )
             ),
             tabPanel("Experimentally validated mutations",
                      sidebarLayout(
                        sidebarPanel(width = 3,
                          textInput("geneInputExp", "Enter Gene:",width = '50%'),  
                          actionButton("plotButtonExp", "Submit")
                        ),
                        mainPanel(align="center",
                          plotOutput("lollipopPlot_EXP"),
                          plotOutput("Protein_lollipop_EXP"),
                          plotOutput("barPlot_EXP"),
                          div(
                            tableOutput("filteredTable_EXP"),  # Table output
                            downloadButton("downloadData_EXP", "Download Data") )
                        )
                      )
             )
  )
)



server <- function(input, output, session) {
  geneDataFiltered_file2 <- reactive({
    req(input$plotButton)
    isolate({
      userGene <- toupper(input$geneInput)
      filtered_data <- file2[file2$Hugo_Symbol == userGene, ]
    })
  })
  geneDataFiltered_file1 <- reactive({
    req(input$plotButton)
    isolate({
      userGene <- toupper(input$geneInput)
      filtered_data <- file1[file1$Hugo_Symbol == userGene, ]
    })
  })
  geneDataFiltered_domain <- reactive({
    req(input$plotButton)
    isolate({
      userGene <- toupper(input$geneInput)
      filtered_data <- domain[domain$Gene == userGene, ]
    })
  })
  geneDataFiltered_exp <- reactive({
    req(input$plotButtonExp)
    isolate({
      userGene <- toupper(input$geneInputExp)
      filtered_data <- exp_val[exp_val$Gene == userGene, ]
    })
  })
  observeEvent(input$plotButton, {
    output$filteredTable <- renderTable({
      geneData2 <- geneDataFiltered_file2()
      
      return(geneData2)
    })
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("filtered_data_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      geneData2 <- geneDataFiltered_file2()
      write.csv(geneData2, file)
    }
  )
  

  observeEvent(input$plotButton, {
    output$lollipopPlot <- renderPlot({
      geneData1 <- geneDataFiltered_file1()
      geneData2 <- geneDataFiltered_file2()
      if (nrow(geneData2) > 0) {
        snp <- c(geneData2$Start_Position)
        values <- c(geneData2$c.DNA_change)
        protein <- c(geneData2$Amino_acid_change)
        site <- c(geneData2$Splice.site)

        Mutations <- GRanges(seqnames = geneData2$Chromosome, IRanges(snp, width = 1, names = paste0(values, ":", protein, "(", site, ")")))
        Mutations$score <- sample.int(15, length(Mutations), replace = TRUE)
        Mutations$color <- sample.int(6, length(snp), replace = TRUE)
        
        features <- GRanges(seqnames = geneData1$Chromosome, IRanges(c(geneData1$Start_Position), width = c(geneData1$Interval)))
        feature_colors <- c("red", "blue", "green", "orange", "purple", "pink")
        
        features$fill <- sample(feature_colors, length(geneData1$Interval), replace = TRUE)
        lolliplot(Mutations, features, lollipop_style_switch_limit = 45,yaxis = FALSE)
      } else {
        plot(NULL, main = "No data available for the specified gene")
      }
    })
  }) 

  observeEvent(input$plotButton, {
    output$Protein_lollipop <- renderPlot({
      geneDatadomain <- geneDataFiltered_domain()
      geneData2 <- geneDataFiltered_file2()
      if (nrow(geneData2) > 0) {
        geneData2$Protein_length <- as.numeric(geneData2$Protein_length)
        geneDatadomain$DS <- as.numeric(geneDatadomain$DS)
        snp <- c(geneData2$Protein_length)
        label <- c(geneData2$c.DNA_change)
        protein_cosmic <- c(geneData2$Amino_acid_change)
        site <- c(geneData2$Splice.site)
        
        
        values <- c(geneDatadomain$DN)
        protein <- c(geneDatadomain$Protein_int)
        
        
        Mutations <- GRanges(geneData2$Chromosome, IRanges(snp, width = 1, names = paste0(label, ":", protein_cosmic, "(", site, ")")))
        Mutations$score <- sample.int(15, length(Mutations), replace = TRUE)
        Mutations$color <- sample.int(6, length(snp), replace = TRUE)
        
        features <- GRanges(geneDatadomain$Chromosome, IRanges(geneDatadomain$DS, width = geneDatadomain$Domain_int,names = paste0(values)))
        feature_colors <- c("red", "blue", "green", "orange", "purple", "pink")
        
        features$fill <- sample(feature_colors, length(protein), replace = TRUE)
        lolliplot(Mutations, features, lollipop_style_switch_limit = 45, yaxis = FALSE)
      } else {
        plot(NULL, main = "No data available for the specified gene")
      }
    })
  })
  
  
  
  observeEvent(input$plotButton, {
    output$barPlot3 <- renderPlot({
      geneData2 <- geneDataFiltered_file2()
      if (nrow(geneData2) > 0) {
        geneData2$Mutations <- factor(geneData2$Mutations)
        site_counts <- table(geneData2$Mutations, geneData2$Splice.site)
        
        # Create a bar plot with assigned colors
        barplot(
          site_counts,
          main = "Distribution of Mutations by Splice Site",
          xlab = "Mutations",
          ylab = "Count",
          col = mutations_colors[unique(geneData2$Mutations)],
          legend.text = TRUE
        )
      } else {
        plot(NULL, main = "No data available for the specified gene")
      }
    })
  })
  
 
  observeEvent(input$plotButtonExp, {
    output$filteredTable_EXP <- renderTable({
      geneData_exp <- geneDataFiltered_exp()
      #head(geneData2, n = input$numEntries)
      return(geneData_exp)
    })
  })
  
  output$downloadData_EXP <- downloadHandler(
    filename = function() {
      paste("filtered_data_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      geneData_exp <- geneDataFiltered_exp()
      write.csv(geneData_exp, file)
    }
  )
  
  observeEvent(input$plotButtonExp, {
    output$lollipopPlot_EXP <- renderPlot({
      geneData1 <- geneDataFiltered_file1()
      geneData_exp <- geneDataFiltered_exp()
      
      if (nrow(geneData_exp) > 0) {
        snp <- c(geneData_exp$Genomic_position)
        values <- c(geneData_exp$c.DNA_change)
        protein <- c(geneData_exp$Amino_acid_change)
        site <- c(geneData_exp$Splice_site)
        
        Mutations <- GRanges(geneData_exp$Chromosome, IRanges(snp, width = 1, names = paste0(values, ":", protein, "(", site, ")")))
        Mutations$score <- sample.int(15, length(Mutations), replace = TRUE)
        Mutations$color <- sample.int(6, length(snp), replace = TRUE)
        
        features <- GRanges(geneData1$Chromosome, IRanges(c(geneData1$Start_Position), width = c(geneData1$Interval)))
        feature_colors <- c("red", "blue", "green", "orange", "purple", "pink")
        
        features$fill <- sample(feature_colors, length(geneData1$Interval), replace = TRUE)
        trackViewer::lolliplot(Mutations, features, lollipop_style_switch_limit = 45, yaxis = FALSE, lollipop_spacing_factor = 0.005)
      } 
      else 
      {
        plot(NULL, main = "No data available for the specified gene")
      }
    })
  })
  
  
  observeEvent(input$plotButtonExp, {
    output$Protein_lollipop_EXP <- renderPlot({
      geneDatadomain <- geneDataFiltered_domain()
      geneData_exp <- geneDataFiltered_exp()
      if (nrow(geneData_exp) > 0) {
        geneData_exp$Protein_length <- as.numeric(geneData_exp$Protein_length)
        geneDatadomain$DS <- as.numeric(geneDatadomain$DS)
        snp <- c(geneData_exp$Protein_length)
        label <- c(geneData_exp$c.DNA_change)
        protein_cosmic <- c(geneData_exp$Amino_acid_change)
        site <- c(geneData_exp$Splice_site)
        
        values <- c(geneDatadomain$DN)
        protein <- c(geneDatadomain$Protein_int)
        
        Mutations <- GRanges(geneData_exp$Chromosome, IRanges(snp, width = rep(1, length(snp)), names = paste0(label, ":", protein_cosmic, "(", site, ")")))
        
        Mutations$score <- sample.int(15, length(Mutations), replace = TRUE)
        Mutations$color <- sample.int(6, length(snp), replace = TRUE)
        
        features <- GRanges(geneDatadomain$Chromosome, IRanges(geneDatadomain$DS, width = geneDatadomain$Domain_int,names = paste0(values)))
        feature_colors <- c("red", "blue", "green", "orange", "purple", "pink")
        
        features$fill <- sample(feature_colors, length(protein), replace = TRUE)
        trackViewer::lolliplot(Mutations, features, lollipop_style_switch_limit = 45, yaxis = FALSE)
      } else {
        plot(NULL, main = "No data available for the specified gene")
      }
    })
  })
  
  
  observeEvent(input$plotButtonExp, {
    output$barPlot_EXP <- renderPlot({
      geneData_exp <- geneDataFiltered_exp()
      if (nrow(geneData_exp) > 0) {
        geneData_exp$Mutation <- factor(geneData_exp$Mutation)
        site_counts <- table(geneData_exp$Mutation, geneData_exp$Splice_site)
        
        barplot(
          site_counts,
          main = "Distribution of Mutations by Splice Site",
          xlab = "Mutations",
          ylab = "Count",
          col = mutations_colors[unique(geneData_exp$Mutation)],
          legend.text = TRUE
        )
      } else {
        plot(NULL, main = "No data available for the specified gene")
      }
    })
  })
  
}

# Run the app
shinyApp(ui, server)

