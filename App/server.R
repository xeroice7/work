#install.packages("shiny")
library(shiny)

server <- function(input, output) {
  output$plot1 <- renderPlot({
      if (any(input$checkGroup == 1)) #A549
          draw.pairwise.venn(a549no1, a549no2, a549no5, category = c("A549 No Biotin", "A549 + Biotin"), lty = rep("solid", 2), 
                        lwd = rep(2, 2), fontface = "bold", cat.fontface = "bold", col = rep("black", 2), fill = c("dark red", "dark blue"), alpha = rep(0.5, 1), 
                        cat.pos = c(-20, 0), cat.dist = rep(0.025, 2), scaled = TRUE)
      
      if (any(input$checkGroup == 2)) #MCF7
          draw.pairwise.venn(mcf7no1, mcf7no2, mcf7no5, category = c("MCF7 No Biotin", "MCF7 + Biotin"), lty = rep("solid", 2), 
                        lwd = rep(2, 2), fontface = "bold", cat.fontface = "bold", col = rep("black", 2), fill = c("deepskyblue", "firebrick1"), alpha = rep(0.5, 1), 
                        cat.pos = c(-10, 0), cat.dist = rep(0.025, 2), scaled = TRUE)
  
      if (any(input$checkGroup == 3)) #MDA-MB-231
          draw.pairwise.venn(mdano1, mdano2, mdano5, category = c("MDA-MB-231 No Biotin", "MDA-MB-231 + Biotin"), lty = rep("solid", 2), 
                        lwd = rep(2, 2), fontface = "bold", cat.fontface = "bold", col = rep("black", 2), fill = c("light green", "magenta"), alpha = rep(0.5, 1), 
                        cat.pos = c(-10, 0), cat.dist = rep(0.025, 2), scaled = TRUE)
      
      if (any(input$checkGroup == 4)) #NCI-H1299
          draw.pairwise.venn(ncino1, ncino2, ncino5, category = c("NCI-H1299 No Biotin", "NCI-H1299 + Biotin"), lty = rep("solid", 2), 
                        lwd = rep(2, 2), fontface = "bold", cat.fontface = "bold", col = rep("black", 2), fill = c("turquoise", "orange"), alpha = rep(0.5, 1), 
                        cat.pos = c(-20, 0), cat.dist = rep(0.025, 2), scaled = TRUE)
    
      if (any(input$checkGroup == 5)) #FiPS4F5
          draw.pairwise.venn(fipsno1, fipsno2, fipsno5, category = c("FiPS4F5 No Biotin", "FiPS4F5 + Biotin"), lty = rep("solid", 2), 
                       lwd = rep(2, 2), fontface = "bold", cat.fontface = "bold", col = rep("black", 2), fill = c("purple", "pink"), alpha = rep(0.5, 1), 
                       cat.pos = c(-20, 0), cat.dist = rep(0.025, 2), scaled = TRUE)
      
      if (any(input$checkGroup == 6)) #HUViPS4F1
          draw.pairwise.venn(huvno1, huvno2, huvno5, category = c("HUViPS4F5 No Biotin", "HUViPS4F5 + Biotin"), lty = rep("solid", 2), 
                       lwd = rep(2, 2), fontface = "bold", cat.fontface = "bold", col = rep("black", 2), fill = c("yellow", "dark green"), alpha = rep(0.5, 1), 
                       cat.pos = c(-20, 0), cat.dist = rep(0.025, 2), scaled = TRUE)
  })
}

