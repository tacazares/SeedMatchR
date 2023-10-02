SeedMatchR.theme = ggplot2::theme(
  title = ggplot2::element_text(size=14),
  plot.background = ggplot2::element_rect(fill = 'white'),
  text = ggplot2::element_text(color = '#34495e', size=12),
  legend.justification = c(0,1),
  legend.box.background = ggplot2::element_rect(fill = 'white', colour = "#2c3e50"),
  panel.grid.major = ggplot2::element_line(colour = "#6b8cae", linetype = 2, linewidth=.2),
  panel.grid.minor = ggplot2::element_line(colour = "white", linetype = 2, linewidth=.2),
  legend.title = ggplot2::element_text(face= "bold", size=12),
  legend.text = ggplot2::element_text(size = 12),
  panel.background = ggplot2::element_rect(fill = 'white', colour = "#2c3e50")
)

# color palette
SeedMatchR.palette = c('black', '#3498db', '#2ecc71','#f1c40f', '#e74c3c',
                       '#9b59b6','#1abc9c', '#f39c12','#d35400')

