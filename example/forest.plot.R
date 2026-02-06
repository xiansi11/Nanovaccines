
forest.plot = function(contrasts, min_val = NULL, max_val = NULL, title = "", nudge_x = NULL){

  df <- contrasts %>%
    # Transform data to tibble
    as_tibble() %>%
    # Add labs for y axis and define linetype according to p.value
    mutate(labs_y = paste0(contrast, ", p=", round(p.value, 4)),
           line_typ = ifelse(p.value <= 0.05, "solid", "dashed")) %>%
    # Change labs for y axis from "-" to "vs"
    mutate_at(vars(labs_y), function(x) gsub(" - ", " vs ", x))
  df$labs_y = factor(df$labs_y, levels = df$labs_y)
  df$contrast = factor(df$contrast, levels = df$contrast)
  
  if(is.null(min_val)){
    min_val = min(floor(df$estimate))
  }
  if(is.null(max_val)){
    max_val = max(ceiling(df$estimate))
  }
  if(is.null(nudge_x)){
    nudge_x = (max_val-min_val)/5
  }
  # plot lsmeans
  df %>%
    ggplot(aes(x = estimate,
               y = contrast)) +
    # Add point for estimate
    geom_point()+
    # Change y axis position and add labs
    scale_y_discrete(position = "right",
                     labels = df$labs_y) +
    # Add error bars, use line_typ to set line type
    geom_errorbar(aes(xmax = estimate + SE,
                      xmin = estimate - SE,
                      linetype = line_typ),
                  show.legend = F,
                  width = 0.3) +
    scale_linetype_manual(values = c("dashed", "solid")) +
    # geom_segment(aes(
    #   x = estimate + SE, xend =estimate + SE,
    #   y = contrast + Y_width, yend = contrast - Y_width
    # ),
    # size = 0.5, linetype = "solid"
    # ) +
    # geom_segment(aes(
    #   x = estimate - SE, xend =estimate - SE,
    #   y = contrast + Y_width, yend = contrast - Y_width
    # ),
    # size = 0.5, linetype = "solid"
    # )+
    # Add estimates values as text
    geom_text(aes(label = round(estimate, 2)), 
              # nudge position in x and y directions
              nudge_x = nudge_x,
              nudge_y = 0.1) +
    # Add dashed vertical line in 0 Difference (y-axis)
    geom_vline(xintercept=0, lty = "dashed") +
    # Make line for second axis (left side)
    geom_vline(xintercept = min_val )  +
    # Change axis titles
    labs(x = "Difference", y = "") +
    scale_x_continuous(position = "top",
                       sec.axis = sec_axis(~., name = ""),
                       limits = c(min_val,max_val),
                       expand = c(0,0)) +
    # Add cowplot theme
    theme_cowplot()+
    # Final theme adjustments
    # Remove ticks and text from bottom secondary axis
    theme(axis.ticks.x.bottom = element_blank(),
          axis.text.x.bottom = element_blank(),
          axis.line = element_line(colour = "black")) +
    ggtitle(title)
}
