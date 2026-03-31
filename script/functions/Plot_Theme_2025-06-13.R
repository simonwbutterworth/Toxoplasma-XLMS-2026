
simon_theme <-
  theme_pubr(
    border=T,
    base_family="Open Sans")+
  theme(
    panel.background = element_blank(),
    axis.ticks = element_line(colour="#000000", linewidth=0.25),
    panel.border=element_rect(colour="#000000", linewidth = 0.5, fill = NA))

#Set theme
lourido_theme <-
  theme_pubr(
    border=T,
    base_family="Open Sans",
    base_size=7)+
  theme(
    panel.background = element_blank(),
    axis.ticks = element_line(colour="#000000", linewidth=0.25),
    axis.ticks.length = unit(3, "pt"),
    axis.text = element_text(size=5),
    panel.border=element_rect(colour="#000000", linewidth = 0.25, fill = NA))
