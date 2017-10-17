install.packages("gplots")
library(gplots)
t <- seq(0, 2 * pi, length = 20)
df <- data.frame(x = sin(t), y = cos(t))
df %>% 
  ggvis() %>% 
  layer_paths(~x, ~y, fill:= "blue", fillOpacity := 0.2) %>%
  layer_paths(~x + 1, ~y, fill:= "red", fillOpacity := 0.2) %>%
  layer_text(~-1.5, ~0.8, text := "HUViPS4F1", fontSize := 20) %>%
  layer_text(~1.65, ~0.8, text := "FiPS4F5", fontSize := 20) %>%
  add_axis("x", properties = axis_props(grid = "none")) %>%
  add_axis("y", properties = axis_props(grid = "none")) %>%
  add_tooltip(~0.4, ~0.8, totalips,"hover")
  
  
  textplot(totalips$Gene, halign="center", valign="center", cex = 0.4, fixed.width = TRUE)
