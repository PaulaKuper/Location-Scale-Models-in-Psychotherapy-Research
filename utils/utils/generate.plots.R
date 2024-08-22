# Function to generate bubble plots plot
generate_bubble_plot <- function(res, dat, xlimMin, xlimMax, ylimMax, moderator, 
                                 xTicks, plotLabel, legendTextSize, xlab, 
                                 legend.box.size= 2, legend.position= "topright",
                                 color, ...) {

  
  # Generate sequence for plotting
  ns100 <- seq(xlimMin, xlimMax, length = xlimMax- xlimMin)
  
  # Predict values based on a model object
  preds <- predict(res,  newscale = ns100)
  
  # Calculate some values for plotting
  hi <- hatvalues(res)
  tau2i <- pmax(0, resid(res)^2 / (1 - hi) - dat$vi)
  wi <- sqrt(weights(res))
  size <- 0.3 + 3 * (wi - min(wi)) / (max(wi) - min(wi))
  
  # Initialize an empty plot
  plot(NA, NA, xlim = c(xlimMin, xlimMax), ylim = c(-0.13, ylimMax), 
       xlab = xlab, xaxt = "n", ylab = "Estimate of \u03C4\u00B2", 
       las = 1, bty = "l", main = plotLabel)
  
  # Add custom x-axis
  axis(1, xTicks, labels =xTicks)
  
  # Predict values based on the model object with specific transformations
  pred <- predict(res, newscale = ns100, transf = exp)
  
  # Add confidence intervals, regression lines, and data points to the plot
  polygon(c(ns100, rev(ns100)), c(pred$ci.lb, rev(pred$ci.ub)), 
          border = NA, col = "gray85")
  #lines(ns100, pred$ci.lb, lty = "dashed", col= "lightblue")
  #lines(ns100, pred$ci.ub, lty = "dashed")
  lines(ns100, pred$pred, lwd = 3, col= color)
  points(dat[[moderator]], tau2i, pch = 21, col = "white", bg = color, cex = size)
  
  # Add legends to the plot
  legend(legend.position, inset = 0.01, bg = "white", 
         pch = c(NA, NA, 22), col = c(NA, NA, color), pt.bg = c(NA, NA, "white"),
         lty = c("blank", NA, NA), lwd = c(NA, NA, 1), text.col = "white", 
         pt.cex = legend.box.size, seg.len = 3, cex= legendTextSize,
         legend = c("Studies", "Regression line", "95% CI"))
  
  legend(legend.position, inset = 0.01, bg = NA, 
         pch = c(21, NA, NA), col = color, pt.bg = c(color, NA, NA), 
         lty = c("blank", "solid", "dashed"), lwd = c(1, 3, 1), 
         text.col = "black", pt.cex = 1, seg.len = 3, cex= legendTextSize,
         legend = c("Studies", "Regression line", "95% CI"))
  

}





# Function to generate density plots for each moderator level
generate_density_plots <- function(res, dat, moderator, plotLabel, legendTextSize, 
                                   x_axis_length, y_axis_length, ...) {
  
  # Set up the plotting environment
  par(mar = c(5, 5, 2, 1))
  
  # Calculate some values for plotting
  hi <- hatvalues(res)
  tau2i <- pmax(0, resid(res)^2 / (1 - hi) - dat$vi)
  wi <- sqrt(weights(res))
  size <- 0.5 + 3 * (wi - min(wi)) / (max(wi) - min(wi))
  
  # Determine levels of the moderator variable
  levels_mod <- levels(dat[[moderator]])
  
  # Determine the number of levels
  num_levels <- length(levels_mod)
  
  # Initialize an empty plot
  plot(NA, NA, xlim = c(0, x_axis_length), ylim = c(0, y_axis_length), 
       xlab = "Estimate of \u03C4\u00B2", ylab = "Density", 
       las = 1, bty = "l", main = plotLabel)
  
  # Generate viridis color palette
  color_palette <- viridis(num_levels)
  
  # Loop through each level of the moderator variable
  for (i in 1:num_levels) {
    # Subset data for the current level
    subset_data <- tau2i[dat[[moderator]] == levels_mod[i]]
    # Plot density
    density_plot <- density(subset_data, ...)
    lines(density_plot, col = color_palette[i], lwd = 3, lty = 1)
  }
  
  # Add legend to the plot
  legend("topright", inset = 0.01, bg = "white", 
         legend = levels_mod,
         col = color_palette, lwd = 2, lty = 1, 
         text.col = "black", cex = legendTextSize)
}



generate_bubble_plot_gg <- function(res, dat, xlimMin, xlimMax, ylimMax, moderator, 
                                 xTicks, plotLabel, legendTextSize, xlab, 
                                 legend.box.size = 2, legend.position = "topright", ...) {
  # Generate sequence for plotting
  ns100 <- seq(xlimMin, xlimMax, length.out = 1250)
  
  # Predict values based on a model object
  preds <- predict(res, newscale = ns100)
  
  # Calculate some values for plotting
  hi <- hatvalues(res)
  tau2i <- pmax(0, resid(res)^2 / (1 - hi) - dat$vi)
  wi <- sqrt(weights(res))
  size <- 0.3 + 3 * (wi - min(wi)) / (max(wi) - min(wi))
  
  # Predict values based on the model object with specific transformations
  pred <- predict(res, newscale = ns100, transf = exp)
  
  # Create data frame for predicted values
  pred_data <- data.frame(
    ns100 = ns100,
    ci_lb = pred$ci.lb,
    ci_ub = pred$ci.ub,
    pred = pred$pred
  )
  
  # Create data frame for the original data
  plot_data <- data.frame(
    moderator = dat[[moderator]],
    tau2i = tau2i,
    size = size
  )
  
  # Create the ggplot
  p <- ggplot() +
    geom_polygon(data = pred_data, aes(x = ns100, y = pred, group = 1), fill = "gray85", color = NA) +
    geom_line(data = pred_data, aes(x = ns100, y = ci_lb), linetype = "dashed") +
    geom_line(data = pred_data, aes(x = ns100, y = ci_ub), linetype = "dashed") +
    geom_line(data = pred_data, aes(x = ns100, y = pred), size = 1) +
    geom_point(data = plot_data, aes(x = moderator, y = tau2i, size = size), shape = 21, color = "black", fill = "darkgray") +
    scale_size_identity() +
    scale_x_continuous(breaks = xTicks, labels = xTicks, limits = c(xlimMin, xlimMax)) +
    scale_y_continuous(limits = c(-0.13, ylimMax)) +
    labs(x = xlab, y = "Estimate of \u03C4\u00B2", title = plotLabel) +
    theme_classic() +
    theme(
      legend.position = legend.position,
      legend.text = element_text(size = legendTextSize)
    ) +
    guides(size = guide_legend(override.aes = list(shape = 21, fill = "darkgray")))
  
  return(p)
}