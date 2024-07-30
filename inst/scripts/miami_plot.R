# miami plot code ====
miami_plot <- function(data,
                       x_axis_group,
                       exposure,
                       outcome,
                       estimate,
                       threshold_pos,
                       threshold_neg,
                       annotate,
                       colour,
                       colours = c("#00378f", "#ffc067", "#894300"),
                       ylims,
                       title = NA,
                       y_breaks = NULL,
                       transformation_from = -0.1,
                       transformation_to = 0.1,
                       transformation = 10
){

  discrete_palette <- colours

  # transformation function ====
  library(scales)
  squish_trans <- function(from, to, factor) {

    trans <- function(x) {

      if (any(is.na(x))) return(x)

      # get indices for the relevant regions
      isq <- x > from & x < to
      ito <- x >= to

      # apply transformation
      x[isq] <- from + (x[isq] - from)/factor
      x[ito] <- from + (to - from)/factor + (x[ito] - to)

      return(x)
    }

    inv <- function(x) {

      if (any(is.na(x))) return(x)

      # get indices for the relevant regions
      isq <- x > from & x < from + (to - from)/factor
      ito <- x >= from + (to - from)/factor

      # apply transformation
      x[isq] <- from + (x[isq] - from) * factor
      x[ito] <- to + (x[ito] - (from + (to - from)/factor))

      return(x)
    }

    # return the transformation
    return(trans_new("squished", trans, inv))
  }


  # format data
  data.tmp <- data
  data.tmp$outcome_label <- data.tmp[[outcome]]
  data.tmp$estimate <- data.tmp[[estimate]]
  data.tmp[[x_axis_group]] <- as.factor(data.tmp[[x_axis_group]])
  data.tmp <- data.tmp[order(data.tmp[[x_axis_group]], data.tmp[[exposure]]),]
  data.tmp$section_column <- factor(data.tmp[[x_axis_group]],
                                    labels = 1:nlevels(data.tmp[[x_axis_group]]))
  data.tmp$x <- with(data.tmp, ave(seq_along(data.tmp$section_column), data.tmp$section_column, FUN = seq_along))


  # Compute section sizes
  data.tmp <- data.tmp %>%
    group_by(section_column) %>%
    summarise(section_column_len=max(x)) %>%

    # Calculate cumulative position of each protein
    mutate(tot=cumsum(section_column_len)-section_column_len) %>%
    select(-section_column_len) %>%

    # Add this info to the initial dataset
    left_join(data.tmp, ., by=c("section_column"="section_column")) %>%

    # Add a cumulative position of each protein
    arrange(section_column, x) %>%
    mutate(xcum = x + tot)

  # get chromosome center positions for x-axis
  axisdata <- data.tmp %>% group_by(section_column) %>% summarize(center=( max(xcum) + min(xcum) ) / 2 )

  # highlights
  data.tmp <- data.tmp %>%
    mutate(is_annotate = ifelse(outcome_label %in% annotate, "yes", "no"))

  # plot ====
  ggplot(data.tmp) +
    geom_point(aes(x = xcum, y = .data[[estimate]], colour = .data[[colour]])) +
    scale_colour_manual(values = discrete_palette) +

    # custom axis
    scale_x_continuous(expand = c(0.01,0.01), label = axisdata$section_column, breaks = axisdata$center ) +
    scale_y_continuous(expand = c(0,0), limits = ylims, trans = squish_trans(transformation_from, transformation_to, transformation), breaks = y_breaks) +

    # title
    ggtitle(paste0(title)) +


    # add lines
    # geom_hline(yintercept = threshold_pos, linetype = "solid", col = "black") +
    # geom_hline(yintercept = threshold_neg, linetype = "solid", col = "black") +
    geom_hline(yintercept = 0, linetype = "solid", col = "black") +

    # Add label using ggrepel to avoid overlapping
    geom_text_repel(data = data.tmp[data.tmp$is_annotate=="yes",],
                    aes(label = outcome_label,
                        x = xcum, y = estimate),
                    min.segment.length = 0, box.padding = 1,
                    force = 1,
                    point.padding = NA, hjust = 1, # this will force the labels in the centre and line them up
                    max.iter = 10000, max.overlaps = 100,
                    colour = "black",
                    seed = 821) +

    # theme
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position="bottom",
      legend.title = element_blank(),
      legend.margin=margin(t=-100),
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank()
    ) +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank())
}
