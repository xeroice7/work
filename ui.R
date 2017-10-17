#install.packages("shiny")
library(shiny)

shinyUI(fluidPage(
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
))
