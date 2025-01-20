#### Combine all the simulation results into one dataframe
library(tidyverse)
library(ggpubr)
library(gridExtra)
library(cowplot)
library(latex2exp)

invlogit = function(y){
  # the inverse-logit function \vphi
  return( 1 / (1 + exp(-y)) )
}

d.invlogit = function(y){
  # calculate the derivative of the inverse-logit function
  # which is also \vphi(1-\vphi)
  return( 1 / (4*cosh(y/2)^2) )
}

#### Debiased inference ####

df = data.frame()

nsim = 200
d = 1000
n = 900
para_rule = c('1se','mincv','minfeas')
alpha_0 = 0.2

# read debiased inference results

for (it in 1){
  for (i in 0:3){
    for (rule in para_rule){
      for (jobid in 1:nsim){
        tryCatch({
          oneres = read.csv(paste0("./debias_res_AR/DebiasProg_AR_cov_d", d, "_n", n, 
                                   "_", jobid, "_x", i, "_rule", rule, "_intercept", it, ".csv"))
          oneres$x = i
          oneres$rule = rule
          # oneres$rule = paste0("debias:",rule)
          
          df = rbind(df, oneres)
          #df.rescale = rbind(df.rescale, oneres)
        }, error = function(e){
          jobid = jobid + 1
        })
      }
    }
  }
}


df = na.omit(df)

coverage = data.frame()
ppresults = data.frame()

for (i in 0:3){
  if (i == 0) {
    ## x0
    x = rep(0, d)
    x[max(d,100)] = 0.2
  }
  if (i == 1) {
    ## x1
    x = rep(0, d)
    x[c(1, 2, 3, 7, 8, 9, 10)] = c(0.1, 0.1, -0.05, 0.05, -0.05, 0.05, -0.1)
  }
  if (i == 2) {
    ## x2
    x = 0.5 ^ seq(1, d, 1)
    x = 0.2 * x / norm(x, "2")
  }
  if (i == 3) {
    ## x3
    x = 1 / seq(1, d, 1) ^ 2
    x = x * (-1) ^ seq(0, d-1, 1)
    x = 0.2 * x / norm(x, "2")
  }
  
  theta_0 = rep(0, d)
  theta_0[1:5] = 2
  theta_0[6:10] = -1
  
  m_true = c(x %*% theta_0) + alpha_0
  #par(mfrow=c(1,3))
  for (rule in c("1se","mincv","minfeas")){
    subdf = df %>% filter(x==i, intercept==1)
    subdf = subdf[subdf['rule']==rule,]
    prob = mean( abs(invlogit(subdf$m_deb) - invlogit(m_true)) / 
                   (d.invlogit(subdf$m_deb) * subdf$asym_sd) < qnorm(0.975) )
    
    normalized_m = (invlogit(subdf$m_deb) - invlogit(m_true)) / 
      (d.invlogit(subdf$m_deb) * subdf$asym_sd)
    cdf_norm = pnorm(sort(normalized_m, decreasing=F), 0,1)
    cdf_m = punif(seq(0,1, length.out=length(normalized_m)))
    
    coverage = rbind(coverage, list(prob=prob, x=i, rule=rule,
                                    CI.upp = min(prob + qnorm(0.975)*sqrt(prob*(1-prob)/nrow(subdf)), 1),
                                    CI.low = max(prob - qnorm(0.975)*sqrt(prob*(1-prob)/nrow(subdf)), 0),
                                    bias_avg = mean(subdf$m_deb) - m_true,
                                    bias_se = sd(subdf$m_deb) / sqrt(nrow(subdf)),
                                    true_mean = m_true,
                                    var_med = mean(subdf$asym_sd),
                                    var_se = sd(subdf$asym_sd) / sqrt(nrow(subdf))))
    #plot(cdf_m, cdf_norm, col="blue",  xlim=c(0,1), ylim=c(0,1), 
    #     main=paste0("x=",i,",theta=",k,",rule=",rule))
    #abline(a=0, b=1, col='red')
    temp_df <- data.frame(
      cdf_m = cdf_m,
      cdf_norm = cdf_norm,
      normalized_m = normalized_m,
      method = paste0("HDLR:", rule),
      x = i
    )
    ppresults <- rbind(ppresults, temp_df)
    
  }
}

for (i in 0:3) {
  plot_data <- ppresults %>% filter(x == i)
  
  p <- ggplot(plot_data, aes(x = cdf_norm, y = cdf_m, color = method)) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, color = 'red') +
    xlim(0, 1) +
    ylim(0, 1) +
    labs(
      title = paste0("x=", i),
      x = "Theoretical CDF",
      y = "Empirical CDF"
    ) +
    theme_minimal() + theme(
      legend.position = c(0.95, 0.05),
      legend.justification = c(1, 0),
      legend.background = element_rect(color = "black")) + guides(color = guide_legend(title = NULL))
  
  print(p)
  Sys.sleep(1)
}

#### Cross-fitted estimator ####

df.cf = data.frame()

for (it in 1){
  for (i in 0:3){
    for (rule in c("1se", "mincv", "minfeas")){
      for (jobid in 1:nsim){
        tryCatch({
          oneres = read.csv(paste0("./debiascf_res_AR/Debiascf_AR_cov_d", d, "_n", n, 
                                   "_", jobid, "_x", i, "_rule", rule, "_intercept", it, ".csv"))
          oneres$x = i
          df.cf = rbind(df.cf, oneres)
        }, error = function(e){
          jobid = jobid + 1
        })
      }
    }
  }
}



df.cf = na.omit(df.cf)

coverage.cf = data.frame()
ppresults.cf = data.frame()

for (i in 0:3){
  if (i == 0) {
    ## x0
    x = rep(0, d)
    x[max(d,100)] = 0.2
  }
  if (i == 1) {
    ## x1
    x = rep(0, d)
    x[c(1, 2, 3, 7, 8, 9, 10)] = c(0.1, 0.1, -0.05, 0.05, -0.05, 0.05, -0.1)
  }
  if (i == 2) {
    ## x2
    x = 0.5 ^ seq(1, d, 1)
    x = 0.2 * x / norm(x, "2")
  }
  if (i == 3) {
    ## x3
    x = 1 / seq(1, d, 1) ^ 2
    x = x * (-1) ^ seq(0, d-1, 1)
    x = 0.2 * x / norm(x, "2")
  }
  
  theta_0 = rep(0, d)
  theta_0[1:5] = 2
  theta_0[6:10] = -1
  
  m_true = c(x %*% theta_0) + alpha_0
  #par(mfrow=c(1,3))
  for (rule in c("1se","mincv","minfeas")){
    subdf = df.cf %>% filter(x==i, intercept==1)
    subdf = subdf[subdf['rule']==rule,]
    prob = mean( abs(invlogit(subdf$m_deb) - invlogit(m_true)) / 
                   (d.invlogit(subdf$m_deb) * subdf$asym_sd) < qnorm(0.975) )
    
    normalized_m = (invlogit(subdf$m_deb) - invlogit(m_true)) / 
      (d.invlogit(subdf$m_deb) * subdf$asym_sd)
    cdf_norm = pnorm(sort(normalized_m, decreasing=F), 0,1)
    cdf_m = punif(seq(0,1, length.out=length(normalized_m)))
    
    coverage.cf = rbind(coverage.cf, list(prob=prob, x=i, rule=rule,
                                    CI.upp = min(prob + qnorm(0.975)*sqrt(prob*(1-prob)/nrow(subdf)), 1),
                                    CI.low = max(prob - qnorm(0.975)*sqrt(prob*(1-prob)/nrow(subdf)), 0),
                                    bias_avg = mean(subdf$m_deb) - m_true,
                                    bias_se = sd(subdf$m_deb) / sqrt(nrow(subdf)),
                                    true_mean = m_true,
                                    var_med = mean(subdf$asym_sd),
                                    var_se = sd(subdf$asym_sd) / sqrt(nrow(subdf))))
    #plot(cdf_m, cdf_norm, col="blue",  xlim=c(0,1), ylim=c(0,1), 
    #     main=paste0("x=",i,",theta=",k,",rule=",rule))
    #abline(a=0, b=1, col='red')
    temp_df <- data.frame(
      cdf_m = cdf_m,
      cdf_norm = cdf_norm,
      normalized_m = normalized_m,
      method = paste0("HDLRcf:",rule),
      x = i
    )
    ppresults.cf <- rbind(ppresults.cf, temp_df)
    
  }
}

for (i in 0:3) {
  plot_data <- ppresults.cf %>% filter(x == i)
  
  p <- ggplot(plot_data, aes(x = cdf_norm, y = cdf_m, color = method)) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, color = 'red') +
    xlim(0, 1) +
    ylim(0, 1) +
    labs(
      title = paste0("x=", i),
      x = "Theoretical CDF",
      y = "Empirical CDF"
    ) +
    theme_minimal() + theme(
      legend.position = c(0.95, 0.05),
      legend.justification = c(1, 0),
      legend.background = element_rect(color = "black")) + guides(color = guide_legend(title = NULL))
  
  print(p)
  Sys.sleep(1)
}

#### SIHR ####

df.SIHR = data.frame()

## read SIHR results
for (jobid in 1:nsim){
  tryCatch({
    oneres = read.csv(paste0("./SIHR_res_AR/SIHR_AR_cov_d", d, "_n", n, "_", jobid, ".csv"))
    df.SIHR = rbind(df.SIHR, oneres)
  }, error = function(e){
    jobid = jobid
  })
}

df.SIHR = na.omit(df.SIHR)
coverage.SIHR = data.frame()
ppresults.SIHR = data.frame()

for (i in 0:3){
  if (i == 0) {
    ## x0
    x = rep(0, d)
    x[max(d,100)] = 0.2
  }
  if (i == 1) {
    ## x1
    x = rep(0, d)
    x[c(1, 2, 3, 7, 8, 9, 10)] = c(0.1, 0.1, -0.05, 0.05, -0.05, 0.05, -0.1)
  }
  if (i == 2) {
    ## x2
    x = 0.5 ^ seq(1, d, 1)
    x = 0.2 * x / norm(x, "2")
  }
  if (i == 3) {
    ## x3
    x = 1 / seq(1, d, 1) ^ 2
    x = x * (-1) ^ seq(0, d-1, 1)
    x = 0.2 * x / norm(x, "2")
  }
  
  theta_0 = rep(0, d)
  theta_0[1:5] = 2
  theta_0[6:10] = -1
  
  m_true = c(x %*% theta_0) + alpha_0
  subdf = df.SIHR %>% filter(x==i)
  prob = mean( abs(invlogit(subdf$m_deb) - invlogit(m_true)) / 
                 (d.invlogit(subdf$m_deb) * subdf$asym_sd) < qnorm(0.975) )
  
  normalized_m = (invlogit(subdf$m_deb) - invlogit(m_true)) / 
    (d.invlogit(subdf$m_deb) * subdf$asym_sd)
  cdf_norm = pnorm(sort(normalized_m, decreasing=F), 0,1)
  cdf_m = punif(seq(0,1, length.out=length(normalized_m)))
  
  coverage.SIHR = rbind(coverage.SIHR, list(prob=prob, x=i, 
                                            CI.upp = min(prob + qnorm(0.975)*sqrt(prob*(1-prob)/nrow(subdf)), 1),
                                            CI.low = max(prob - qnorm(0.975)*sqrt(prob*(1-prob)/nrow(subdf)), 0),
                                            bias_avg = mean(subdf$m_deb) - m_true,
                                            bias_se = sd(subdf$m_deb) / sqrt(nrow(subdf)),
                                            true_mean = m_true,
                                            var_med = mean(subdf$asym_sd),
                                            var_se = sd(subdf$asym_sd) / sqrt(nrow(subdf))))
  temp_df <- data.frame(
    cdf_m = cdf_m,
    normalized_m = normalized_m,
    cdf_norm = cdf_norm,
    method = 'SIHR',
    x = i
  )
  ppresults.SIHR <- rbind(ppresults.SIHR, temp_df)
}

for (i in 0:3) {
  plot_data <- ppresults.SIHR %>% filter(x == i)
  
  p <- ggplot(plot_data, aes(x = cdf_norm, y = cdf_m)) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, color = 'red') +
    xlim(0, 1) +
    ylim(0, 1) +
    labs(
      title = paste0("x=", i),
      x = "Theoretical CDF",
      y = "Empirical CDF"
    ) +
    theme_minimal() + theme(
      legend.position = c(0.95, 0.05),
      legend.justification = c(1, 0),
      legend.background = element_rect(color = "black"))
  
  print(p)
  Sys.sleep(1)
}

#### Refit and oracle ####

df.refit = data.frame()

## read refit results
for (jobid in 1:nsim){
  tryCatch({
    oneres = read.csv(paste0("./refit_res_AR/refit_AR_cov_d", d, "_n", n, "_", jobid, ".csv"))
    df.refit = rbind(df.refit, oneres)
  }, error = function(e){
    jobid = jobid + 1
  })
}

df.refit = na.omit(df.refit)
coverage.refit = data.frame()
ppresults.refit = data.frame()

for (method in c("oracle", "refit")){
  for (i in 0:3){
    if (i == 0) {
      ## x0
      x = rep(0, d)
      x[max(d,100)] = 0.2
    }
    if (i == 1) {
      ## x1
      x = rep(0, d)
      x[c(1, 2, 3, 7, 8, 9, 10)] = c(0.1, 0.1, -0.05, 0.05, -0.05, 0.05, -0.1)
    }
    if (i == 2) {
      ## x2
      x = 0.5 ^ seq(1, d, 1)
      x = 0.2 * x / norm(x, "2")
    }
    if (i == 3) {
      ## x3
      x = 1 / seq(1, d, 1) ^ 2
      x = x * (-1) ^ seq(0, d-1, 1)
      x = 0.2 * x / norm(x, "2")
    }
    
    theta_0 = rep(0, d)
    theta_0[1:5] = 2
    theta_0[6:10] = -1
    
    m_true = c(x %*% theta_0) + alpha_0
    subdf = df.refit %>% filter(x==i)
    subdf = subdf[subdf['method']==method,]
    subdf$asym_sd[subdf$asym_sd==0] = 1e-4
    prob = mean( abs(invlogit(subdf$m_est) - invlogit(m_true)) / 
                   (d.invlogit(subdf$m_est) * subdf$asym_sd) < qnorm(0.975) )
    
    normalized_m = (invlogit(subdf$m_est) - invlogit(m_true)) / 
      (d.invlogit(subdf$m_est) * subdf$asym_sd)
    cdf_norm = pnorm(sort(normalized_m, decreasing=F), 0,1)
    cdf_m = punif(seq(0,1, length.out=length(normalized_m)))
    
    coverage.refit = rbind(coverage.refit, list(prob=prob, x=i, method=method,
                                                CI.upp = min(prob + qnorm(0.975)*sqrt(prob*(1-prob)/nrow(subdf)), 1),
                                                CI.low = max(prob - qnorm(0.975)*sqrt(prob*(1-prob)/nrow(subdf)), 0),
                                                bias_avg = mean(subdf$m_est) - m_true,
                                                bias_se = sd(subdf$m_est) / sqrt(nrow(subdf)),
                                                true_mean = m_true,
                                                var_med = mean(subdf$asym_sd),
                                                var_se = sd(subdf$asym_sd) / sqrt(nrow(subdf))))
    plot(cdf_m, cdf_norm, col="blue",  xlim=c(0,1), ylim=c(0,1),
         main=paste0("x=",i,"method=",method))
    abline(a=0, b=1, col='red')
    temp_df <- data.frame(
      cdf_m = cdf_m,
      cdf_norm = cdf_norm,
      normalized_m = normalized_m,
      method = method,
      x = i
    )
    ppresults.refit <- rbind(ppresults.refit, temp_df)
    Sys.sleep(0.5)
  }
}



#### Plots ####

titles <- c(
  expression("Query Point" ~ x^{(0)}),
  expression("Query Point" ~ x^{(1)}),
  expression("Query Point" ~ x^{(2)}),
  expression("Query Point" ~ x^{(3)})
)

calculate_x_limits <- function(means, std_errors) {
  min_limit <- min(means - std_errors) * 1.1
  max_limit <- max(means + std_errors) * 1.1
  c(min(0, min_limit), max_limit)
}

draw_bar_plot <- function(means, std_errors, truncation_level, colors, title, show_labels = TRUE) {
  # Adjust means for truncation
  adjusted_means <- rev(pmax(pmin(means, truncation_level), -truncation_level))  # Truncate means and reverse
  adjusted_std_errors <- rev(ifelse(means > truncation_level | means < -truncation_level, 0, std_errors))  # Suppress error bars and reverse
  is_truncated <- rev(means > truncation_level | means < -truncation_level)  # Reverse truncation flags
  
  # Define bar colors and reverse for consistency with reversed datasets
  bar_colors <- rev(colors)  # Reverse the colors to match reversed datasets
  
  # Dynamically calculate x-axis limits
  x_limits <- calculate_x_limits(adjusted_means, adjusted_std_errors)
  
  # Update y-axis labels and reverse order
  dataset_labels <- rev(c("HDLR:1se", "HDLR:mincv", "HDLR:minfeas",
                          "HDLRcf:1se", "HDLRcf:mincv", "HDLRcf:minfeas",
                          "oracle", "refit", "SIHR"))  # Reverse the order of labels
  
  # Create horizontal bar plot with reversed y-axis labels and consistent colors
  bar_midpoints <- barplot(
    adjusted_means,
    horiz = TRUE,  # Horizontal bar plot
    col = bar_colors,  # Fill the bars with reversed colors
    border = "black",  # Keep borders for clarity
    xlim = x_limits,  # Dynamically calculated x-axis limits
    main = title,  # Set the plot title
    axes= FALSE,     # don't plot x-axis here, plot it later (below)
    cex.main = 2,  # Increase title size (adjust this value as needed)
    width = 0.8,  # Adjust bar width
    space = 0.2  # Adjust spacing between bars
  )
  
  # title(sub = subtitle, cex.sub = 1.5)
  
  # Draw bars with rugged edges for truncated bars pointing outward
  zigzag_width <- 0.02  # Width of zigzag
  zigzag_density <- 8    # Number of zigzag peaks along the bar height
  
  for (i in seq_along(adjusted_means)) {
    y_bottom <- bar_midpoints[i] - 0.4  # Adjust for consistent width
    y_top <- bar_midpoints[i] + 0.4
    x_right <- adjusted_means[i]
    x_left <- 0
    
    if (is_truncated[i]) {
      # Adjust rugged edges to point outward
      if (means[9 - i + 1] > 0) {
        # Positive truncated bar
        y_zigzag <- seq(y_bottom, y_top, length.out = zigzag_density)
        x_zigzag <- rep(c(truncation_level, truncation_level + zigzag_width), length.out = length(y_zigzag))  # Extend outward beyond truncation level
      } else {
        # Negative truncated bar
        y_zigzag <- seq(y_bottom, y_top, length.out = zigzag_density)
        x_zigzag <- rep(c(-truncation_level, -truncation_level - zigzag_width), length.out = length(y_zigzag))  # Extend outward below truncation level
      }
    } else {
      y_zigzag <- c(y_bottom, y_top)
      x_zigzag <- c(x_right, x_right)
    }
    
    # Generate coordinates for rugged and non-rugged bars
    x_coords <- c(x_left, x_zigzag[1], x_zigzag, x_zigzag[length(x_zigzag)], x_left)
    y_coords <- c(y_bottom, y_bottom, y_zigzag, y_zigzag[length(y_zigzag)], y_top)
    
    # Draw the bar with rugged edge
    polygon(x_coords, y_coords, col = bar_colors[i], border = "black")
    
    
    # Add numerical value with standard deviation for truncated bars inside the bar
    if (is_truncated[i]) {
      label <- paste0(
        round(means[9 - i + 1], 2),  # Correct indexing for means
        " (", 
        round(std_errors[9 - i + 1], 2),  # Correct indexing for standard deviations
        ")"
      )
      
      # Display the label inside the bar
      if (means[9 - i + 1] > 0) {
        text(truncation_level * 0.9, bar_midpoints[i], labels = label, cex = 1, col = "white", pos = 2)  # Position near truncation level
      } else {
        text(-truncation_level * 0.9, bar_midpoints[i], labels = label, cex = 1, col = "white", pos = 4)  # Position near -truncation level
      }
    }
    
  }
  
  # Add vertical line at zero
  abline(v = 0, col = "gray", lty = 2)
  
  # Add error bars only for non-truncated bars
  valid_indices <- which(!is_truncated)
  arrows(adjusted_means[valid_indices] - adjusted_std_errors[valid_indices], bar_midpoints[valid_indices], 
         adjusted_means[valid_indices] + adjusted_std_errors[valid_indices], bar_midpoints[valid_indices], 
         angle = 90, code = 3, length = 0.06, col = "black")
  
  # Draw a full box around the plot
  box(col = "black", lwd = 1)  # Set box line width to match axis
  
  # Add ticks and labels for the left side of the box without the axis line
  axis(
    side = 2,                       # Specify the left side of the plot
    at = bar_midpoints,             # Position ticks at bar midpoints
    labels = if (show_labels) dataset_labels else FALSE,  # Show labels only for leftmost plot
    las = 1,                        # Align labels horizontally
    cex.axis = 1.6,                 # Adjust label size for readability
    tck = -0.02,                    # Draw small ticks inward
    lwd.ticks = 1,                  # Ensure ticks are visible
    lwd = 0                         # Suppress the y-axis line
  )
  
  # Add ticks and labels for the x-axis adaptively
  tick_positions <- pretty(x_limits, n = 7)  # Generate 7 adaptive tick positions based on the data range
  axis(
    side = 1,                        # Specify the bottom side of the plot
    at = tick_positions,             # Adaptive tick positions
    labels = format(tick_positions, digits = 2),  # Adaptive labels, formatted for readability
    cex.axis = 1.5,                  # Increase label size for x-axis
    tck = -0.02,                     # Ensure consistent tick size with y-axis
    lwd.ticks = 1                    # Ensure ticks are visible
  )
  
}

draw_coverage_plot <- function(group_data, colors, title, show_labels = TRUE) {
  group_data <- group_data[nrow(group_data):1, ]
  
  # Calculate the lower limit adaptively based on the data
  x_lower_limit <- floor(min(group_data$CI.low) / 0.025) * 0.025  # Round down to nearest 0.025 increment
  x_ticks <- seq(x_lower_limit, 1, by = 0.025)  # Generate ticks in 0.025 increments
  
  # Update y-axis labels and reverse order
  dataset_labels <- group_data$method  # Reverse the order of labels
  
  # Define bar colors and reverse for consistency with reversed datasets
  coverage_colors <- rev(colors)  # Reverse the colors to match reversed datasets
  
  # Create the plot
  plot(
    group_data$prob, seq_along(group_data$method), 
    xlim = c(x_lower_limit, 1),
    ylim = c(0.5, 9.5),  # Set limits
    xlab = "",  # Remove x-axis label
    ylab = "",  # Remove y-axis label
    main = title,
    pch = 20,
    col = coverage_colors,
    cex.main = 2,  # Increase title size (adjust this value as needed)
    axes = FALSE  
  )
  
  # Add horizontal error bars with thicker lines
  for (j in seq_along(group_data$prob)) {
    arrows(
      x0 = group_data$CI.low[j], x1 = group_data$CI.upp[j],
      y0 = j, y1 = j,  # Adjust y-position
      angle = 90, code = 3, length = 0.1, col = coverage_colors[j], lwd = 2  # Increase line width with `lwd`
    )
  }
  
  # Add vertical dashed line at 0.95
  abline(v = 0.95, col = "gray", lty = 2)
  
  # Draw a full box around the plot
  box(col = "black", lwd = 1)  # Set box line width to match axis
  
  
  # Add ticks and labels for the left side of the box without the axis line
  axis(
    side = 2,                       # Specify the left side of the plot
    at = seq_along(group_data$method),  # Position ticks at bar midpoints
    labels = if (show_labels) dataset_labels else FALSE,  # Show labels only for leftmost plot
    las = 1,                        # Align labels horizontally
    cex.axis = 1.8,                 # Adjust label size for readability
    tck = -0.02,                    # Draw small ticks inward
    lwd.ticks = 1,                  # Ensure ticks are visible
    lwd = 0                         # Suppress the y-axis line
  )
  
  # Add ticks and labels for the x-axis adaptively
  axis(
    side = 1,                        # Specify the bottom side of the plot
    at = x_ticks,             # Adaptive tick positions
    labels = format(x_ticks, digits = 2),  # Adaptive labels, formatted for readability
    cex.axis = 1.8,                  # Increase label size for x-axis
    tck = -0.02,                     # Ensure consistent tick size with y-axis
    lwd.ticks = 1                    # Ensure ticks are visible
  )
  
  
  
  # Add axes
  #axis(1, at = x_ticks, labels = format(x_ticks, digits = 2))  # x-axis with ticks in 0.025 increments
  #if (show_labels) {
  #  axis(2, at = seq_along(group_data$method), labels = dataset_labels, las = 1)  # y-axis for leftmost plot
  #} else {
  #  axis(2, at = seq_along(group_data$method), labels = FALSE)  # Suppress y-axis labels for other plots
  #}
  
  
}

## generate QQ plots
# Update the function call to generate correct labels for each QQ plot
generate_qq_plot <- function(data, colors = NULL, title = "QQ Plot", subtitle = NULL,
                             legend_position = "topleft", line_spacing = 1) {
  # Check inputs
  #if (!is.list(data_list) || length(data_list) == 0) stop("data_list must be a non-empty list.")
  #if (length(data_list) > 3) stop("Up to 3 datasets are allowed.")
  #if (is.null(labels) || length(labels) != length(data_list)) stop("Provide labels matching data_list.")
  #if (!is.null(colors) && length(colors) != length(data_list)) stop("Provide colors matching data_list.")
  
  # Set default colors if not provided
  #if (is.null(colors)) {
  #  colors <- c("blue", "red", "green")[seq_along(data_list)]
  #}
  
  # Get the unique groups from the data
  methods <- unique(data$method)
  
  # Initialize the plot with fixed axes
  plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1),
       xlab = "", ylab = "", main = "") 
  
  # Add diagonal reference line
  abline(a = 0, b = 1, col = "black", lty = 2, lwd = 2)
  
  # Overlay QQ plots for each dataset
  for (i in rev(seq_along(methods))) {
    group_data <- data[data$method == methods[i], ] 
    n <- length(group_data$normalized_m)
    
    # Theoretical and sample probabilities
    theoretical_probs <- ppoints(n)
    sample_probs <- pnorm(sort(group_data$normalized_m))
    
    # Plot points within [0, 1]
    points(theoretical_probs, sample_probs, col = colors[i], pch = 16)
  }
  
  # Add legend
  legend(legend_position, legend = methods, col = colors, pch = 16, 
         bty = "n", box.lwd = 0.2, cex = 0.8, 
         text.width = max(strwidth(methods))*0.4)
  
  # Add title
  title(main = title, sub = subtitle, line = line_spacing, cex.main = 1.5, cex.sub = 1.2)
}

bar_colors <- c("HDLR:1se" = 'blue', "HDLR:mincv" = "deepskyblue", "HDLR:minfeas" = "skyblue",
                "HDLRcf:1se" = 'red', "HDLRcf:mincv" = "coral", "HDLRcf:minfeas" = "pink1",
                "oracle" = 'grey', "refit" = "purple", "SIHR" = "violet")  # example colors



#### QQ plots ####

par(mfrow = c(2,2), mar = c(2, 0, 3, 0), oma = c(0.5, 0.5, 0.5, 0.5), pty = "s")
for (i in 0:3){
  subdata.HDLR = ppresults %>% filter(x==i)
  generate_qq_plot(subdata.HDLR, colors = bar_colors[1:3], title = titles[i+1])
  
  subdata.cf = ppresults.cf %>% filter(x==i)
  generate_qq_plot(subdata.cf, colors = bar_colors[4:6], title = titles[i+1])
  
  subdata.SIHR = ppresults.SIHR %>% filter(x==i)
  generate_qq_plot(subdata.SIHR, colors = bar_colors[9], title = titles[i+1])
  
  subdata.refit = ppresults.refit %>% filter(x==i)
  generate_qq_plot(subdata.refit, colors = bar_colors[7:8], title = titles[i+1])
  Sys.sleep(0.5)
}

par(mfrow=c(1,1))

#### Variance ####

variances = data.frame()
var_deb = coverage[c("x","var_med","var_se")]
var_deb$method = rep(c("HDLR:1se","HDLR:mincv","HDLR:minfeas"), 4)
variances = rbind(variances, var_deb)

var_cf = coverage.cf[c("x", "var_med", "var_se")]
var_cf$method = rep(c("HDLRcf:1se","HDLRcf:mincv","HDLRcf:minfeas"), 4)
variances = rbind(variances, var_cf)

var_refit = coverage.refit[c("x","var_med","var_se","method")]
variances = rbind(variances, var_refit)

var_SIHR = coverage.SIHR[c("x","var_med","var_se")]
var_SIHR$method = "SIHR"
variances = rbind(variances, var_SIHR)

# Plot the 4 groups of datasets
par(mfrow = c(1, 4), mar = c(2, 2, 4, 2), oma = c(0, 9, 0, 0))  # Increased left margin
for (i in 1:4) {
  subvars = variances %>% filter(x==i-1)
  draw_bar_plot(
    means = subvars$var_med,
    std_errors = subvars$var_se,
    truncation_level = 0.5,
    colors = bar_colors,
    title = titles[i],
    show_labels = (i == 1)  # Only show y-axis labels on the first plot
  )
}



#### Coverage of CI ####

cover_prob = data.frame()
cover_deb = coverage[c("x", "prob","CI.upp","CI.low")]
cover_deb$method = rep(c("HDLR:1se","HDLR:mincv","HDLR:minfeas"), 4)
cover_prob = rbind(cover_prob, cover_deb)

cover_cf = coverage.cf[c("x", "prob", "CI.upp", "CI.low")]
cover_cf$method = rep(c("HDLRcf:1se","HDLRcf:mincv","HDLRcf:minfeas"), 4)
cover_prob = rbind(cover_prob, cover_cf)

cover_refit = coverage.refit[c("x", "prob", "CI.upp", "CI.low", "method")]
cover_prob = rbind(cover_prob, cover_refit)

cover_SIHR = coverage.SIHR[c("x", "prob", "CI.upp", "CI.low")]
cover_SIHR$method = "SIHR"
cover_prob = rbind(cover_prob, cover_SIHR)

par(
  mfrow = c(1, 4),          # 1 row and 4 columns of plots
  mar = c(2, 2, 4, 2),      # Increased left margin (second value in `mar`)
  oma = c(0, 10, 0, 0)       # Adjusted outer margins
)

for (i in 1:4) {
  # Get the data for the current group
  group_data = cover_prob %>% filter(x==i-1)
  
  # Call the function to draw a single plot
  draw_coverage_plot(
    group_data,
    colors = bar_colors,
    title = titles[i],
    #labels = c("HDLR:1se", "HDLR:mincv", "HDLR:minfeas",
    #           "HDLR:1se", "HDLRcf:mincv", "HDLRcf:minfeas",
    #           "oracle", "refit", "SIHR"),
    show_labels = (i == 1)  # Only show y-axis labels on the first plot
  )
}



#### Bias ####

biases = data.frame()
bias_deb = coverage[c("x", "bias_avg", "bias_se")]
bias_deb$method = rep(c("debias:1se","debias:mincv","debias:minfeas"), 4)
biases = rbind(biases, bias_deb)

bias_cf = coverage.cf[c("x", "bias_avg", "bias_se")]
bias_cf$method = rep(c("HDLRcf:1se","HDLRcf:mincv","HDLRcf:minfeas"), 4)
biases = rbind(biases, bias_cf)

bias_refit = coverage.refit[c("x", "bias_avg", "bias_se", "method")]
biases = rbind(biases, bias_refit)

bias_SIHR = coverage.SIHR[c("x", "bias_avg", "bias_se")]
bias_SIHR$method = "SIHR"
biases = rbind(biases, bias_SIHR)

par(mfrow = c(1, 4), mar = c(2, 2, 4, 2), oma = c(0, 9, 0, 0))  # Increased left margin
for (i in 1:4) {
  subbias = biases %>% filter(x==i-1)
  draw_bar_plot(
    means = subbias$bias_avg,
    std_errors = subbias$bias_se,
    truncation_level = 0.25,
    colors = bar_colors,
    title = titles[i],
    show_labels = (i == 1)  # Only show y-axis labels on the first plot
  )
}



#### Others ####

nsim = 250
d = 1000
n = 900
para_rule = c('1se','mincv','minfeas')
alpha_0 = 0.2

df.eg = data.frame()
for (i in 0:0){
  for (k in 0:2){
    for (rule in para_rule){
      for (jobid in 1:nsim){
        tryCatch({
          oneres = read.csv(paste0("./debias_res_AR/DebiasProg_AR_cov_d", d, "_n", n, 
                                   "_", jobid, "_x", i, "_theta", k, "_rule", rule, ".csv"))
          oneres$x = i
          oneres$theta = k
          oneres$rule = rule
          
          df.eg = rbind(df.eg, oneres)
        }, error = function(e){
          jobid = jobid + 1
        })
      }
    }
  }
}

df.eg = na.omit(df.eg)

# generate histogram
subdf = df.eg %>% filter(x==0, theta==2)
subdf = subdf[subdf["rule"]=="minfeas",]
subdf_long = data.frame(
  value = c(invlogit(subdf$m_cur), invlogit(subdf$m_deb)),
  group = rep(c("Lasso","Debiased"), each=nrow(subdf))
)

x = rep(0, d)
#x[c(1, 2, 3, 7, 8)] = c(1, 1/2, 1/4, 1/2, 1/8)
#x[c(1, 2, 3, 7, 8, 9, 10)] = c(0.1, 0.1, -0.05, 0.05, -0.05, 0.05, -0.1)
x[1] = 0.2
#s_theta = 25
theta_0 = rep(0, d)
#theta_0[1:s_theta] = sqrt(1)
theta_0[1:5] = sqrt(c(9,7,5,3,1)) * (-1)^(0:4)
m_true = c(x %*% theta_0) + alpha_0

ggplot(subdf_long) + geom_histogram(aes(x=value, y = after_stat(density), 
                                        fill=group), binwidth = 0.01) + 
  geom_vline(xintercept = invlogit(m_true), col="red", linetype = 'dashed') + labs(x = 'probability') + 
  theme_minimal() + theme(text = element_text(size = 15), legend.title = element_blank())

normalized_m = (invlogit(subdf$m_deb) - invlogit(m_true)) / (subdf$asym_sd * d.invlogit(subdf$m_deb))
normalized_lasso = (invlogit(subdf$m_cur) - invlogit(m_true)) / sd(invlogit(subdf$m_cur))

norm_df = data.frame(
  value = c(normalized_m, normalized_lasso),
  group = rep(c("Debiased", "Lasso"), each=nrow(subdf))
)

ggplot(norm_df) + geom_histogram(aes(x=value, y = after_stat(density), 
                                     fill=group), binwidth = 0.25, position = 'identity', alpha = 0.6) + 
  # facet_wrap(~group, ncol = 1) + 
  stat_function(fun = dnorm, args = list(mean = 0, sd = 1), color = "black", linetype = "dashed") + 
  xlim(-5,5) + ylim(0, 0.6) + theme_minimal() + labs(x = "z-score") + 
  theme(legend.position = "right", text = element_text(size = 15)) + 
  scale_fill_manual(values = c("red", "blue"))
