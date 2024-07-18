# HARP_samplecode.R
# Version 2.2 (released July 2024)

# Author: Melanie E Roberts
# Date: June 2022
# Exemplar code for HARP Hysteresis metrics

# This project is in collaboration with Donghwan Kim, Jing Lu and David Hamilton
# If you use this code please cite this GitHub repository 
# and the associated paper Melanie E. Roberts, Donghwan Kim, Jing Lu & David P. Hamilton,
# HARP: A suite of parameters to describe the hysteresis of streamflow and water quality constituents,
# Journal of Hydrology, 2023, doi:10.1016/j.jhydrol.2023.130262.
# 
# ------------------------------------------------------------------------------
# CHANGE LOG
#  8 February 2024 - correction to plotting labels with no figure 8
#  13 March 2024 - Version 2.0 Improved interpolation
#  18 May 2024 - Version 2.1 Changed how figure eight loop area signs are calculated
#  20 July 2024 - Version 2.2 Corrections to bugs introduced with previous update
#
# ------------------------------------------------------------------------------


# 1. load required libraries
library(tidyverse)
library(viridis)
library(zoo)
library(dplyr)
library(shotGroups)
library(caTools)
library(latex2exp)
library(docstring) 
library(ggplot2)


filename = "sampleData_fig8.csv"
data_column_names =  c("elapsedTime", "flux", "concentration")  # names of the columns in the dataset in order elapsedTime, flux and concentration 

constituent_name = "TSS"
eventname = "Sample"

plot_filename = "samplePlot.pdf"  # name for the output file of the plot

HARP <- function(fileName, eventName, constituentName, plotFileName){
  #' HARP
  #'
  #' Function to produce a HARP plot for a given dataset
  #' 
  #' fileName - name of the csv file containing the data 
  
  
  # 1. Read in the observations
  
  df_obs <- read.csv(filename, header = TRUE)  # dataframe of observations: columns elapsedTime, flux, concentration
  df_obs <- df_obs %>% rename(  # rename columns
    Etime = elapsedTime,
    Q = flux,
    C = concentration
  )
  
  df_obs[] <- lapply(df_obs, as.numeric)  # ensure that is recognised as numeric
  
  # OPTIONAL plot to explore the original data
  # ggplot(data = df_obs, aes(x = Q, y = C, color = Etime)) +
  #   geom_path() +
  #   geom_point()
  
  
  # 2. Set up metric_df to store calculated metrics
  metric_df <- data.frame(matrix(NA, nrow = 1, ncol = 7))
  columnNames <- c(
    "area",  # area metric (+) diluting, (-) enriching
    "residual",  # residual metric, (+) finish higher than start, (-) finish lower than start
    "peak_C",  # proportional timing peak concentration is reached
    "peak_Q",  # proportional timing peak flux is reached
    "peaktime_Q",  # time (unscaled) of peak flux, intermediate calculation
    "peaktime_C",  # time (unscaled) of peak concentration, intermediate calculation
    "rad_minmax"  # radius of circle with equivalent area to curve, intermediate calculation
  )
  colnames(metric_df) <- columnNames
  
  # 3. Functions used for data wrangling:
   # genCircle - generates the circle in the plots for Area metric
   # interpolated - to interpolate the datasets to a finer grid
   # minmaxNormalisation - for normalising data to [0,1]
   # normalise_observations_flux - to normalise the flux observations based on the interpolated dataset truncated to a closed loop
   # normalise_observations_C - to normalise the concentration observations based on the interpolated dataset truncated to a closed loop
  
  genCircle <-   function(rad = 1.0,  
                          center = c(0.5, 0.5),
                          npoints = 100) {
    #' genCircle
    #' 
    #' used to generate the circle that represents the Area metric in plots
    #' 
    #' Default generates a circle centred at (0.5, 0.5) with radius 1.  
    
    r = rad
    tt <- seq(0, 2 * pi, length.out = npoints)
    xx <- center[1] + r * cos(tt)
    yy <- center[2] + r * sin(tt)
    return(data.frame(x = xx, y = yy))
  }
  
  interpolatedX <- function(i, x, y, ynew) {
    #' interpolatedX
    #' 
    #' a function to interpolate the data to a finer grid wrt the flux (y-value)
    #' using the method of linear interpolation between observations
    #' 
    #' i - index of the current value in the dataframe
    #' x - x-value of the observed data (typically time)
    #' y - y-value of the observed data (typically flux)
    #' ynew - the new y-value to be interpolated to (finer grid, typically flux)
    #' 
    m <- (y[i + 1] - y[i]) / (x[i + 1] - x[i])  # calculate the gradient
    c <- y[i] - m * x[i]  # calculate the intercept
    
    Xnew <- (ynew - c) / m  # interpolated x-value
    return(Xnew)
  }
  
  interpolate_data <- function(x, y, ynew) {
    #' interpolate_data
    #' 
    #' a function to work through the observed data and interpolate between subsequent points
    #' here the new values are in the response (typically flux) and we seek the x-values (typically time)
    #' x - x-values of the observed data (typically time)
    #' y - y-values of the observed data (typically flux)
    #' ynew - the new y-value to be interpolated to (finer grid, typically flux)
    #' 
    mydata <- data.frame(X = numeric(), y_new = numeric()) 
    
    for (i in 1:(length(y) - 1)) {
      # iterate through the list of y-values (typically flux) and interpolate between each pair of observations
      ystart <- y[i]
      yend <- y[i + 1]
      
      if (i == length(y) - 1) {  # if we are at the last point
        if (ystart < yend) {
          Y_inbetween_start_and_end <- which(ynew >= ystart & ynew <= yend)  # identify the points in the new y-values that are between the pair of observations
        } else {  # yend < ystart
          Y_inbetween_start_and_end <- which(ynew >= yend & ynew <= ystart)
        }
      } else { # if we are not at the last point
        if (ystart < yend) {
          Y_inbetween_start_and_end <- which(ynew >= ystart & ynew < yend)
        } else { # yend < ystart
          Y_inbetween_start_and_end <- which(ynew > yend & ynew <= ystart)
        }
      }
      
      for (j in 1:(length(Y_inbetween_start_and_end))){
        pos_in_ynew <- Y_inbetween_start_and_end[j]
        newX <- interpolatedX(i, x, y, ynew[pos_in_ynew] )
        mydata <- rbind(mydata, data.frame(X = newX, y_new = ynew[pos_in_ynew]))
      }
    }
    
    mydata <- mydata %>% drop_na() %>% arrange(X)  # remove any NA values and reorder the dataframe
    return(mydata)
  }
  
  minmaxNormalisation <- function(theData) {
    #' minmaxNormalisation
    #' 
    #' Produces the min-max normalised version of the dataset
    #' 
    #' theData - column in the dataframe that is to be normalised: df$column_name
    
    miny = min(theData, na.rm = TRUE)
    maxy = max(theData, na.rm = TRUE)
    
    return((theData - miny) / (maxy - miny))
  }
  
  normalise_observations_flux <- function(df) { 
    #' normalise_observations_flux
    #' 
    #' min-max normalisation of the observed fluxes
    #'  using information from the interpolated, truncated flux data
    #' 
    #' df - dataframe that has a column Qi, which is the interpolated flux data
    #'   and column Q, which is the observed flux data
    
    miny = min(df$Qi, na.rm = TRUE)
    maxy = max(df$Qi, na.rm = TRUE)
    mydata = df$Q
    return((mydata - miny) / (maxy - miny))
  }
  
  normalise_observations_C <- function(df) {
    #' normalise_observations_C
    #' 
    #' min max normalisation specifically on the observed concentration 
    #'  using information from the interpolated, truncated concentration data 
    #' 
    #' df - dataframe that has a column Ci, which is the interpolated concentration data 
    #'  and column C, which is the observed concentration data
    
    miny = min(df$Ci, na.rm = TRUE)
    maxy = max(df$Ci, na.rm = TRUE)
    mydata = df$C
    return((mydata - miny) / (maxy - miny))
  }
  
  # 4. Scale, interpolate and truncate the dataset
  # a. interpolate and truncate to a closed loop in Q
  maxQ = max(df_obs$Q)
  switchpoint = which(df_obs$Q == maxQ)
  rising_branch = df_obs[1:switchpoint,] # df for just the rising branch
  
  endpoint = length(df_obs$Q)  # df for just the falling branch
  falling_branch <- df_obs[(switchpoint):endpoint,]
  
  minQ_rise = min(rising_branch$Q)
  minQ_fall = min(falling_branch$Q)
  
  max_of_the_mins = max(c(minQ_rise, minQ_fall))  # flux value at which to truncate for closed loop
  
  # interpolate along the rising and falling branches
  Qlist = seq(max_of_the_mins, maxQ, length = 500) 
  
    # rising branch
  
  tmp_time <- interpolate_data(rising_branch$Etime, rising_branch$Q, Qlist)  #  (x_observations, y_observations, y_new)  Interpolate time given finer Q list for interpolation
  dfi_rise <- data.frame("Itime" = tmp_time$X, "Qi" = tmp_time$y_new)  # df for interpolated data rising branch
  
  tmp_conc <- approx(x = rising_branch$Etime, y = rising_branch$C, xout = dfi_rise$Itime, method = "linear", rule = 2) # interpolate concentration using interpolated time,. rule 2 ensures that numerical rounding errors with the x_out does not lead to NA values as it allows extrapolation
  dfi_rise$Ci <- tmp_conc$y  # add the interpolated concentration to the dataframe
  
    # falling branch
  tmp_time <- interpolate_data(falling_branch$Etime, falling_branch$Q, Qlist)  #  (x_observations, y_observations, y_new)  Interpolate time given finer Q list for interpolation
  dfi_fall <- data.frame("Itime" = tmp_time$X, "Qi" = tmp_time$y_new)  # df for interpolated data falling branch
  
  tmp_conc <- approx(x = falling_branch$Etime, y = falling_branch$C, xout = dfi_fall$Itime, method = "linear", rule = 2) #   interpolate concentration using interpolated time,  rule 2 ensures that numerical rounding errors with the x_out does not lead to NA values as it allows extrapolation
  dfi_fall$Ci <- tmp_conc$y  # add the interpolated concentration to the dataframe
  
  dfi = rbind(dfi_rise, dfi_fall)  # df for interpolated data
  dfi = unique(dfi) # remove the duplicate at the turning point
  
  # sort in time
  dfi = dfi[order(dfi$Itime),]
  row.names(dfi)<- NULL
  
  # scale the "dfi" dataframe
  dfi$Qsi = minmaxNormalisation(dfi$Qi)
  dfi$Csi = minmaxNormalisation(dfi$Ci)
  
  # observations and interpolated data
  df_all <- full_join(dfi, df_obs, by = c("Itime" = "Etime"))
  firsttime <- min(dfi_rise$Itime)  # remove observations outside closed loop
  lasttime <- max(dfi_fall$Itime)
  df_all <- df_all[df_all$Itime >= firsttime,]  
  df_all <- df_all[df_all$Itime <= lasttime,]  
  
  df_all = df_all[order(df_all$Itime),]
  row.names(df_all)<- NULL
  
  # renormalise
  df_all$Qsi <- minmaxNormalisation(df_all$Qi)
  df_all$Csi <- minmaxNormalisation(df_all$Ci)
  
  df_all$Qs <- normalise_observations_flux(df_all)
  df_all$Cs <- normalise_observations_C(df_all)
  
  # 5. Calculate metrics for all example events
  
  # for ease I need to set what some of my columns will be
  # Q - flux unscaled
  # Qi - interpolated flux unscaled
  # Qs - min/max normalised flux [NOT REQUIRED]
  # Qsi - min-max normalised flux interpolated
  
  # C - constituent concentration (e.g. TSS, NO3, TN)
  # Ci - interpolated unscaled
  # Csi - interpolated min-max normalised
  
  # Etime - elapsed time since start of event
  # Itime - interpolated time
  
  
  # df_obs - original data frame
  # dfi - data frame that incorporates all the interpolations
  # df_all - data frame that combines the observation and interpolated data with scalings etc.
  
  
  
  # Peaks
  peaktime_flux <- dfi$Itime[dfi$Qsi == 1]
  if (length(peaktime_flux > 1)) {
    peaktime_flux <- peaktime_flux[[1]]
  }
  peaktime_flux_proportion <- peaktime_flux / max(dfi$Itime) * 100
  
  metric_df[1, "peaktime_Q"] = peaktime_flux
  metric_df[1, "peak_Q"] = peaktime_flux_proportion
  
  peaktime_C <- dfi$Itime[dfi$Csi == 1]
  if (length(peaktime_C > 1)) {
    peaktime_C <- peaktime_C[[1]]
  }
  peaktime_C_proportion <- peaktime_C / max(dfi$Itime) * 100
  
  metric_df[1, "peaktime_C"] = peaktime_C
  metric_df[1, "peak_C"] = peaktime_C_proportion
  
  
  # Area
  
  # a. do we have a fig 8?
  # locating intercepts
  dfi_rise <- left_join(dfi_rise, dfi, by = c("Qi", "Ci", "Itime"))  # to get the scaled info
  dfi_fall <- left_join(dfi_fall, dfi, by = c("Qi", "Ci", "Itime"))
  
  segment_intercept_finder <- function(rising_point1_x, rising_point1_y, rising_point2_x, rising_point2_y, falling_point1_x, falling_point1_y, falling_point2_x, falling_point2_y) {
    #' This function identifies intersections between the rising and falling branches (i.e. locating figure 8s)
    #' The inputs are the start and end points of the segment in the rising and falling branch
    #' The parameteric form of an equation along each segment is used to find potential intersections
    #' The function returns the coordinates of the intersection point if it exists, otherwise NULL
    
    # Define the coordinates of the points
    R_x1 <- rising_point1_x;    R_y1 <- rising_point1_y
    R_x2 <- rising_point2_x;    R_y2 <- rising_point2_y
    F_x1 <- falling_point1_x;   F_y1 <- falling_point1_y
    F_x2 <- falling_point2_x;   F_y2 <- falling_point2_y
    
    # Define the coefficients of the linear equations
    A <- matrix(c(R_x2 - R_x1, R_y2 - R_y1, F_x1 - F_x2,  F_y1 - F_y2), nrow = 2)
    B <- matrix(c(F_x1 - R_x1, F_y1 - R_y1), nrow = 2)
    
    # Check if the determinant is zero
    if (is.na(det(A))) {
      return(NULL) # Lines do not intersect or are collinear
    } 
    else if (det(A) == 0) {
      return(NULL)  # Lines do not intersect or are collinear
    } else {
      # Solve for the parameters in the parametric form of the equation
      t_u <- solve(A, B)
      t <- t_u[1]
      u <- t_u[2]
      
      # Check if t and u lie between 0 and 1, i.e. whether the lines intersect between the end points of each segment
      if (t >= 0 && t <= 1 && u >= 0 && u <= 1) {
        # Calculate the intersection point
        intersection_x <- R_x1 + t * (R_x2 - R_x1)
        intersection_y <- R_y1 + t * (R_y2 - R_y1)
        
        return(c(intersection_x, intersection_y)) # Return the intersection point
      } else {
        return(NULL)  # The intersection point does not lie within the line segments
      }
    }
  }
  
  fig_8_intercepts <- function(rising_limb, falling_limb) {
    # Initialize a list to store intersection points
    intersections <- list()
    
    A <- data.frame(
      X = rising_limb$Q,
      Y = rising_limb$C
    )
    
    B <- data.frame(
      X = falling_limb$Q,
      Y = falling_limb$C
    )
    
    # Loop through all adjacent points in A and B
    for (i in 1:(nrow(A) - 1)) {
      for (j in 1:(nrow(B) - 1)) {
        # Get the coordinates of the adjacent points
        
        x1 <- A$X[i]
        y1 <- A$Y[i]
        x2 <- A$X[i + 1]
        y2 <- A$Y[i + 1]
        
        x3 <- B$X[j]
        y3 <- B$Y[j]
        x4 <- B$X[j + 1]
        y4 <- B$Y[j + 1]
        
        # Find the intersection point
        intersection <- segment_intercept_finder(x1, y1, x2, y2, x3, y3, x4, y4)
        
        # If an intersection exists, store it
        if (!is.null(intersection)) {
          intersections <- append(intersections, list(intersection))
        }
      }
    }
    
    # return as a dataframe
    df_ <- do.call(rbind, intersections)
    df_ <- as.data.frame(df_)
    
    if (length(df_) == 0){
          return(NULL)
        }
        else {
          # Set column names
          colnames(df_) <- c("Q", "C")
          
          intersections <- df_
          return (intersections)
        }
        
      }
  
  scale_the_intercepts <- function(df_intercepts, df_All){
    #' scaled_intercepts
    #' This function scales the calculated intercepts of the loops using the truncated, interpolated dataset
    min_for_Q_scaling <- min(df_All$Qi, na.rm = TRUE)
    max_for_Q_scaling <- max(df_All$Qi, na.rm = TRUE)
    
    min_for_C_scaling <- min(df_All$Ci, na.rm = TRUE)
    max_for_C_scaling <- max(df_All$Ci, na.rm = TRUE)
    
    # scale the intercepts
    # Initialize the intercept_scaled data frame
    intercept_scaled <- data.frame(
      Qs = (df_intercepts$Q - min_for_Q_scaling) / (max_for_Q_scaling - min_for_Q_scaling),
      Cs = (df_intercepts$C - min_for_C_scaling) / (max_for_C_scaling - min_for_C_scaling)
    )
    return(intercept_scaled)
  }
  
  get_sign_of_larger_area <- function(Area1, Area2){  # identify whether the postive or negative loops are collectively larger
    if (abs(Area1) > abs(Area2)){
      return(sign(Area1))
    } else {
      return(sign(Area2))
    }
  }
  
  # identify the intercepts (if any) and then min-max normalise in Q and C
  theIntercepts <- fig_8_intercepts(rising_branch, falling_branch)  # working with original observations
  
  eps = 1e-6 # remove intercepts at the turning point from rising to falling branch
  theIntercepts <- theIntercepts[!(abs(theIntercepts$Q - maxQ) < eps | abs(theIntercepts$Q - max_of_the_mins) < eps), ]
  
  theIntercepts_scaled <- scale_the_intercepts(theIntercepts, df_all)  # scale the intercepts relative to df_all
  
  # OPTIONAL plot to observe location of calculated intercepts
  # ggplot(data = df_obs, aes(x = Q, y = C, color = Etime)) +
  #   geom_path() +
  #   geom_point() + geom_point(data = theIntercepts, aes(x = Q, y = C), color = "red", size = 3)

  
  if (is.null(theIntercepts) || nrow(theIntercepts) == 0){
    fig8 = FALSE}  # no intercepts found
  else
  {fig8 = TRUE} # intercepts found
  
  split_by_intersections <- function(rise_df, fall_df, intersections) {
    intersection_indices_rise <- which(rise_df$Qsi %in% intersections$Qs)  # identify the intersection indices in the rising branch based on the Q values
    intersection_indices_fall <- which(fall_df$Qsi %in% intersections$Qs)  # " falling branch
    
    # rising branch
    segments_rise <- list() 
    
    if (length(intersection_indices_rise) == 1) {
      segments_rise[[1]] <- rise_df[1:intersection_indices_rise[1], ]
    }
    else {
      # Loop over each pair of intersection points to create segments 
      for (i in seq_along(intersection_indices_rise)) {
        if (i == 1) {
          # First segment: from start to first intersection point
          segments_rise[[i]] <- rise_df[1:intersection_indices_rise[i], ]
        } else {
          # Middle segments: from one intersection point to the next
          segments_rise[[i]] <- rise_df[intersection_indices_rise[i-1]:intersection_indices_rise[i], ]
        }
      }
    }
    # Last segment: from last intersection point to end
    segments_rise[[length(intersection_indices_rise) + 1]] <- rise_df[intersection_indices_rise[length(intersection_indices_rise)]:nrow(rise_df), ]
    
    # falling branch
    segments_fall <- list() 
    if (length(intersection_indices_fall) == 1) {
      segments_fall[[1]] <- fall_df[1:intersection_indices_fall[1], ]
    }
    else {
      for (i in seq_along(intersection_indices_fall)) {
        if (i == 1) {
          segments_fall[[i]] <- fall_df[1:intersection_indices_fall[i], ]
        } else {
          segments_fall[[i]] <- fall_df[intersection_indices_fall[i-1]:intersection_indices_fall[i], ]
        }
      }
    }
    segments_fall[[length(intersection_indices_fall) + 1]] <- fall_df[intersection_indices_fall[length(intersection_indices_fall)]:nrow(fall_df), ]
    
    segments <- list(rise = segments_rise, fall = segments_fall)
    
    return(segments)
  }
  
  
  fig8_loop_flag = 0
  if (fig8){ # checks if any intercepts were found
    fig8_loop_flag = 1  # there is at least one loop
    
    for (i in 1:nrow(theIntercepts)){  # add each intercept to the dfi_rise and dfi_fall dataframes
      
      myIntercept <- theIntercepts[i, ]
      myIntercept_scaled <- theIntercepts_scaled[i, ]
      intercept_info <- c(NA, myIntercept$Q, myIntercept$C, myIntercept_scaled$Qs, myIntercept_scaled$Cs)
      
      # add this point to dfi dataframe & reorder
      dfi_rise[(nrow(dfi_rise) + 1),] = intercept_info
      
      dfi_rise = dfi_rise[order(dfi_rise$Qsi),]
      row.names(dfi_rise)<- NULL
      dfi_rise$Itime <- na.approx(dfi_rise$Itime)
      
      rise_rowindex <- which(dfi_rise$Qsi == myIntercept_scaled[1,]$Qs)
      
      dfi_fall[(nrow(dfi_fall) + 1),] = intercept_info 
      dfi_fall = dfi_fall[order(dfi_fall$Qsi),]
      row.names(dfi_fall)<- NULL
      dfi_fall$Itime <- na.approx(dfi_fall$Itime)
      
      fall_rowindex <- which(dfi_fall$Qsi == myIntercept_scaled[1,]$Qs)
      
      dfi[nrow(dfi) + 1,] = dfi_rise[rise_rowindex,]  # add the intercept point for the rise and the fall (different times) to the dfi dataframe  
      dfi[nrow(dfi) + 1,] = dfi_fall[fall_rowindex,]
      
      dfi <- dfi[order(dfi$Itime),] # reorder in increasing time
      row.names(dfi) <- NULL
      
    }
    
    # Split the data frame
    segments <- split_by_intersections(dfi_rise, dfi_fall, theIntercepts_scaled)  # use the interpolated rising and falling branches together with the scaled intercepts to split the df into each individual loop
    
    # Form each loop and calculate the areas
    loops <- list()
    
    for (i in 1: length(segments$rise)){  # create the loops
      
      loop_ <- rbind(segments$rise[[i]], segments$fall[[i]])
      loop_ <- loop_[order(loop_$Itime),]
      
      loops[[i]] <- loop_
    }
    
    negAreas = 0  # need to sum up neg and pos separately
    posAreas = 0
    
    for (i in 1:length(loops)){
      loop <- loops[[i]]
      area <- trapz(loop$Qsi, loop$Csi)
      if (area < 0){
        negAreas = negAreas + area
      } else {
        posAreas = posAreas + area
      }
    }
    
    total_area <- abs(negAreas) + abs(posAreas)  # add the area magnitudes together
    area_sign <- get_sign_of_larger_area(negAreas, posAreas)  # determine if predominantly enriching (-) or diluting (+)
    area <- area_sign * total_area  # include sign with area metric
    
    
  } else {  # single loop 
    area = trapz(dfi$Qsi, dfi$Csi)
  }
  
  metric_df[1, "area"] = area
  
  # radius of equivalent circle
  equivRadMinMax <- sqrt(abs(area) / pi) 
  metric_df[1, "rad_minmax"] = equivRadMinMax
  
  # Residual
  residual = dfi$Csi[nrow(dfi)] - dfi$Csi[1] # positive means enriching
  metric_df[1, "residual"] = residual
  

  # 6. Plot
  circ_equiv_minmax <- genCircle(equivRadMinMax)
  
  mycols <- c("Itime", "Q", "C", "Qs", "Cs")
  obs_df <- df_all %>% select(one_of(mycols)) %>% drop_na
  
  df_all$Qall_o <- ifelse(is.na(df_all$Qi), df_all$Q, df_all$Qi  )
  df_all$Call_o <- ifelse(is.na(df_all$Ci), df_all$C, df_all$Ci  )
  
  df_all$Qall_s <- ifelse(is.na(df_all$Qsi), df_all$Qs, df_all$Qsi  )
  df_all$Call_s <- ifelse(is.na(df_all$Csi), df_all$Cs, df_all$Csi  )
  
  
  interpolated_df <- df_all %>% select(Itime, Qall_o, Call_o, Qall_s, Call_s) 
  
  # set up the metric information for plot
  if(fig8_loop_flag == 1){ # figure eight loop
    plot_title = sprintf("%s: %s",eventname, constituent_name)
    xlabel <- "Normalised flux"
    ylabel <- paste("Normalised", constituent_name)
    line1 = sprintf("$A = %0.2f^\\infty$, $R = %0.3f$",
                    area,
                    residual              
    )
    
    line2 = sprintf("$\\hat{Q} = %0.1f$,$ \\hat{%s} = %0.1f$",
                    peaktime_flux_proportion,
                    constituent_name,
                    peaktime_C_proportion
    )
  } else {  # no loops
    plot_title = sprintf("%s: %s",eventname, constituent_name)
    xlabel <- "Normalised flux"
    ylabel <- paste("Normalised", constituent_name)
    line1 = sprintf("$A = %0.2f$, $R = %0.3f$",
                    area,
                    residual               
    )
    
    line2 = sprintf("$\\hat{Q} = %0.1f$,$ \\hat{%s} = %0.1f$",
                    peaktime_flux_proportion,
                    constituent_name,
                    peaktime_C_proportion
    )
  }
  
  theXlabel = "Normalised discharge Q"
  theYlabel = sprintf("Normalised concentration %s", constituent_name)
  
  legend_label <- "Time"
  
  plt = ggplot(data = interpolated_df, mapping = aes(x = Qall_s, y = Call_s)) +  # scatter plot with path
    geom_path(aes(color = Itime))  + 
    scale_color_viridis(discrete = FALSE) +
    geom_point(data = obs_df, aes(x = Qs, y = Cs, color = Itime)) + # this takes only the original data
    geom_path(data = circ_equiv_minmax,
              mapping = aes(x, y),
              color = 'deeppink3',
              linewidth = 0.3) +  # equivalent area circle
    geom_abline(color = 'gray77', linewidth = 0.1) + 
    coord_fixed() +
    ggtitle(TeX(line1)) +  xlab(theXlabel) + ylab(theYlabel) +
    labs(color = legend_label, subtitle = TeX(line2)) +
    theme(plot.title = element_text(size = 9, hjust = 0, , margin=margin(b=0)), plot.subtitle = element_text(size = 9, hjust = 0, , margin=margin(b=0)), legend.title = element_blank(), legend.margin=margin(l = -2, unit='pt'), axis.text.x=element_blank(),  axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  
  ggsave(plotFileName, plot = plt, width = 10, height = 10, units = c("cm"), dpi = 300, limitsize = TRUE)
  
  return(metric_df)
}


# call the function
harp = HARP(fileName = filename, eventName = eventname, constituentName = constituent_name, plotFileName = plot_filename )

print("HARP has completed")
