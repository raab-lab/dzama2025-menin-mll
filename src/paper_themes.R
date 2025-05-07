# paper specific themes/colors
paper_theme <- function() {
  theme(plot.title = element_text(family ="sans", color = 'black', 
                                  face = "bold", size = 10),
        axis.title = element_text(family = "sans", color = 'black', 
                                  size = 8),
        axis.text = element_text(size = 6, color = 'black', family = 'sans'), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linetype = 'dashed', color = 'grey90'),
        panel.border = element_rect(color = 'grey10', fill = NA),
        panel.background = element_blank() )
}
