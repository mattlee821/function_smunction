# devtools::install_github("gforge/forestplot")
library(forestplot)
library(ggplot2)
library(gridExtra)
library(dplyr)

data <- read.table("functions/scripts/forestplot/figure1_test.csv", sep = ",", header = T)
plot_data <- data
head(plot_data)

# forestplot ====
plot_data$HR <- round(exp(plot_data$loghr),2)
plot_data$lower_ci <- plot_data$HR - (1.96 * plot_data$se)
plot_data$upper_ci <- plot_data$HR + (1.96 * plot_data$se)
plot_data$label_HR <- paste0(plot_data$HR, " (", round(plot_data$lower_ci,2), ", ", round(plot_data$upper_ci,2), ")")
plot_data$label_HR <- gsub("NA \\(NA, NA\\)", " ", plot_data$label_HR) # Replace "NA (NA, NA)" with " "
plot_data$label_HR <- gsub("1 \\(1, 1\\)", "1 (reference)", plot_data$label_HR) # Replace "1 (1, 1)" with "1 (reference)"

## plot
tiff("/Users/leem/Downloads/figure1_test.tiff",
     units = "px", width = 1600, height = 800)
plot_data %>%
forestplot(
  labeltext = c(group_label, alc_label, cases, studies,
                label_HR, pwald, ptrend, phet_study, phet_sex),
  align = c("r", "l", "c", "c", "l", "c", "c", "c", "c"),
  mean = HR,
  lower = lower_ci,
  upper = upper_ci,
  xlab = "HR estimates (log scale)",
  xlog = T,
  graph.pos = 5,
  boxsize = 0.25,
  shapes_gp = fpShapesGp(default = gpar(lwd = 2)),
  txt_gp = fpTxtGp(xlab = gpar(cex = 1.1), 
                   ticks = gpar(cex = 1.1),
                   label = gpar(cex = 1.1)))  |>
  fp_add_header(alc_label = "",
                cases = "Cases" |> fp_txt_plain() |> fp_align_center(),
                studies = "Studies" |> fp_txt_plain() |> fp_align_center(),
                label_HR = "HR (95% CI)" |> fp_txt_plain() |> fp_align_center(),
                pwald = list(expression(P[wald^2])) |> fp_txt_plain() |> fp_align_center(),
                ptrend = list(expression(P[trend^3])) |> fp_txt_plain() |> fp_align_center(),
                phet_study = list(expression(P[study^4])) |> fp_txt_plain() |> fp_align_center(),
                phet_sex = list(expression(P[sex^5])) |> fp_txt_plain() |> fp_align_center()
                ) 
dev.off()
