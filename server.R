#install.packages("shiny")
library(shiny)

server <- function(input, output) {
  output$plot1 <- renderPlot({
      if (any(input$checkGroup == 1))
          #A549
          draw.pairwise.venn(a549no1, a549no2, a549no5, category = c("A549 No Biotin", "A549 + Biotin"), lty = rep("solid", 2), 
                       lwd = rep(2, 2), fontface = "bold", cat.fontface = "bold", col = rep("black", 2), fill = c("yellow", "dark green"), alpha = rep(0.5, 1), 
                       cat.pos = c(-10, 0), cat.dist = rep(0.025, 2), scaled = TRUE)
    
  })
}

