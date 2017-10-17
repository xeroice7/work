#install.packages("shiny")
library(shiny)

ui <- fluidPage(
  titlePanel("Surface Proteome Commonalities"),
                
  sidebarLayout(
    sidebarPanel(
      #helptext("Choose which cell line you want to examine: "),
        checkboxGroupInput("checkGroup", label = h3("Cell Lines:"),
                            choices = list("A549" = 1, 
                                          "MCF7" = 2,
                                          "MDA-MB-231" = 3,
                                          "NCI-H1299" = 4,
                                          "FiPS4F5" = 5,
                                          "HUViPS4F1" = 6)
                          ),
        hr(),
        fluidRow(column(2, verbatimTextOutput("value")))
                      
        ),
    mainPanel(
      plotOutput("plot1")
    )
  )
)

server <- function(input, output) {
  output$plot1 <- renderPlot({
    if (any(input$checkGroup == 1))
      #A549
      a549 <- data.frame(read.csv(paste(file.path, "A549SummaryFDRContaminants.csv", sep=""), header = TRUE, sep = ",", stringsAsFactors = FALSE))
      splitcolumns <- str_split_fixed(a549$X, "/", 2)
      a549 <- cbind.data.frame(splitcolumns, a549)
      names(a549)[names(a549) == '1'] <- 'Line'  #Change the first column to "Line"
      names(a549)[names(a549) == '2'] <- 'YN'  #Change the second column to "YN"
      splitcolumns <- str_split_fixed(a549$Line, " ", 3)
      a549 <- cbind.data.frame(splitcolumns, a549)
      names(a549)[names(a549) == '1'] <- 'Cell'  #Change the first column to "Cell"
      names(a549)[names(a549) == '3'] <- 'Dox'  #Change the first column to "Line"
      splitcolumns <- str_split_fixed(a549$YN, " ", 2)
      splitcolumns <- splitcolumns[,1]
      a549 <- cbind.data.frame(splitcolumns, a549)
      names(a549)[names(a549) == "splitcolumns"] <- "Biotin"
      a549 <- subset(a549, select=-c(`2`, Line, YN, X))
      a549 <- filter(a549, Species == "HUMAN") #Filtering only human hits
      a549 <- filter(a549, Dox == "No Dox") #Filtering only No Dox hits
      splitcolumns <- str_split_fixed(a549$Gene.Names," ", 2)
      a549 <- cbind.data.frame(splitcolumns, a549)
      names(a549)[names(a549) == '1'] <- 'Gene'  #Change the first column to "Cell"
      names(a549)[names(a549) == '2'] <- 'Synonyms'  #Change the first column to "Line"
      a549 <- subset(a549, select=-Gene.Names)
      a549slim <- data.frame(a549$Cell, a549$Gene, a549$Synonyms, a549$Biotin, a549$Unused, a549$Total, a549$Gene.Ontology)
      colnames(a549slim) <- c("Line", "Gene", "Synonyms", "Biotin", "Unused", "Total", "Gene Ontology")
    
      a549biotin <- filter(a549slim, Biotin == "+") #Keep only +Biotin for MCF7
      colnames(a549biotin)[colnames(a549biotin) == 'Unused'] <- 'UnusedA549Biotin'
      colnames(a549biotin)[colnames(a549biotin) == 'Total'] <- 'TotalA549Biotin'
    
      a549nobiotin <- filter(a549slim, Biotin == "No") #Keep only No Biotin for FiPS
      colnames(a549nobiotin)[colnames(a549nobiotin) == 'Unused'] <- 'UnusedA549NoBiotin'
      colnames(a549nobiotin)[colnames(a549nobiotin) == 'Total'] <- 'TotalA549NoBiotin'
    
      a549 <- merge.data.frame(a549biotin, a549nobiotin, by = "Gene", all = TRUE)
      a549$UnusedA549Biotin <- as.numeric(as.character(a549$UnusedA549Biotin))
      a549$UnusedA549NoBiotin <- as.numeric(as.character(a549$UnusedA549NoBiotin))
      a549$TotalA549Biotin <- as.numeric(as.character(a549$TotalA549Biotin))
      a549$TotalA549NoBiotin <- as.numeric(as.character(a549$TotalA549NoBiotin))
    
      a549$UnusedA549Biotin <- na.zero(a549$UnusedA549Biotin)
      a549$UnusedA549NoBiotin <- na.zero(a549$UnusedA549NoBiotin)
      a549$TotalA549Biotin <- na.zero(a549$TotalA549Biotin)
      a549$TotalA549NoBiotin <- na.zero(a549$TotalA549NoBiotin)
    
      df6 <- data.frame(a549nobiotin$Gene, a549nobiotin$Line)
      colnames(df6) <- c("Gene", "Line")
      df7 <- data.frame(a549biotin$Gene, a549biotin$Line)
      colnames(df7) <- c("Gene", "Line")
    
      df8 <- setdiff(df6, df7)
      no6 <- nrow(df6)
      df9 <- setdiff(df7, df6)
      no7 <- nrow(df7)
      df10 <- semi_join(df6, df7)
      no10 <- nrow(df10)
      
      grid.newpage()
      draw.pairwise.venn(no6, no7, no10, category = c("MDA No Biotin", "MDA + Biotin"), lty = rep("solid", 2), 
                         lwd = rep(2, 2), fontface = "bold", cat.fontface = "bold", col = rep("black", 2), fill = c("yellow", "dark green"), alpha = rep(0.5, 1), 
                         cat.pos = c(-10, 0), cat.dist = rep(0.025, 2), scaled = TRUE)
  
      })
}

shinyApp(ui = ui, server = server)
