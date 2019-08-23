# DESCRIPTION:
#   Plot QC metrics from fMRI such as Frame Displacement (FD), or DVARS.
#
# Inputs:   data,       A node x volume/tp matrix
#
# Ouput:    p,           A ggplot list that can be print()/ggsave() or feed 
#                        into grid.arrange()
# #########################################################################
# myc & epp, UTD 2019/08/22 - Initial
# #########################################################################
plot_time_series <- function(data, min=NULL, max=NULL, colpal="Greys") {

  # Convert TP matrix to long format for plotting
  df <- data.frame(roi=1:nrow(data),data)
  names(df)[2:ncol(df)] <- 1:(ncol(df)-1)
  df_long <- gather(df, key = "vol", value = "bold", 2:(ncol(df)))
  df_long$vol <- as.numeric(df_long$vol)
  
  # Setup the Plot
  if(is.null(min) | is.null(max)){
    min <- min(df_long$bold)
    max <- max(df_long$bold)
  }
  max_diff <- c(min, max) %>% abs() %>% max()
  
  # Plot
  ggplot(df_long, aes(vol, roi, fill = bold)) +
    geom_raster() +
    scale_fill_distiller(palette = colpal, limits = c(-max_diff, max_diff),
                         guide = guide_colorbar(label = TRUE, frame.colour="black")) + # add border around colorbar)
    scale_y_reverse(expand=c(0, 0)) + 
    scale_x_continuous(expand = c(0, 0)) +
    theme_minimal() 
}

# DESCRIPTION:
#   Plot 6 motion paremters from fMRI (e.g., rp_file from SPM).
#    - X, Y, Z, P, Y, R
#
# Inputs:   data,       A dataframe with 6 motion parameters with column names 
#                       by X, Y, Z, P, Y.1, and R
#
# Ouput:    p,           A ggplot list that can be print()/ggsave() or feed 
#                        into grid.arrange()
# #########################################################################
# myc, UTD 2019/08/22 - Initial
# #########################################################################
plot_motion <- function(data){
  
  # Convert motion data frame to wide-format
  data$vol <- 1:nrow(data) # add a column for "vol"
  df <- gather(data, key = "XYZPYR", value="motion", 1:6) 
  df$XYZPYR<- factor(df$XYZPYR, c("X","Y","Z","P","Y.1","R"))
  df$vol <- as.numeric(df$vol)
  
  # Plotting
  ggplot(df, aes(x=vol, y=XYZPYR, fill=motion)) +
    geom_raster() +
    scale_fill_distiller(palette = "Greys") +
    scale_x_continuous(expand = c(0, 0)) +
    theme_minimal() +
    theme(legend.position = "none",
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
}

# DESCRIPTION:
#   Plot QC metrics from fMRI such as Frame Displacement (FD), or DVARS.
#
# Inputs:   qc,           A numeric vector of qc metric (e.g., FD, DVARS).
#           qc_thres,     A numeric value indicating the threshold to flag a volume. The # of 
#                         volumes exceeding the threshold will be printed in the figure
#           qc_name,      Name of the QC metric (Default = "QC") 
#           qc_color,     Color of the line graph (e.g.,"red","#hex value"). (Default="red")
#           miny/maxy,    The min/max of the y-axis. (Default = min/max of the qc values)
#
# Ouput:    p,           A ggplot list that can be print()/ggsave() or feed 
#                        into grid.arrange()
# #########################################################################
# myc, UTD 2019/08/22 - Initial
# #########################################################################
plot_qc <- function(qc, qc_thres=NULL, qc_name="QC", qc_color="red", miny=NULL, maxy=NULL){
  df <- data.frame(vol=1:length(qc), qc=qc)
  
  if(is.null(miny) | is.null(maxy)){
    miny=min(df$qc)
    maxy=max(df$qc)
  }
  
  p <- ggplot(data = df, aes(x=vol, y=qc)) +
    geom_line(color=qc_color) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0), limits = c(miny,maxy)) +
    ylab(label = qc_name)
  
  if(!is.null(qc_thres)){
    qc_flag_x <- which(df$qc > qc_thres)
    qc_text <- sprintf("%s > %.2f = %d vol", qc_name, qc_thres, length(qc_flag_x))
    
    p <- p + geom_hline(yintercept = qc_thres, size=0.2, linetype="dashed") +
      annotate("text", -Inf, Inf, label=qc_text, hjust = 0, vjust = 1, color=qc_color)
  }
  p <- p + theme_minimal() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  return(p)
}