#' display_tables
#' @description To show the ranking table by inspecting the results of Gibbs sampler.
#' @param file_path path to the saved BIT Gibbs sampling results.
#' @param output_path path to save the rank table.
#' @param burnin number of samples used for burn-in. If not specify, BIT will use the half of the iterations as burn in.
#' @return a data.frame object contains TR names, theta_i, BIT scores and 95 CIs for each TR.
#' @export
display_tables <- function(file_path, output_path, burnin = NULL) {

  # Validate inputs
  if (!file.exists(file_path)) stop("The specified file path does not exist.")
  if (!dir.exists(output_path)) stop("The specified output path does not exist.")

  # Load data
  cat("Loading data from file...\n")
  dat <- readRDS(file_path)

  TR_names <- dat[["TR_names"]]
  theta_i_mat <- dat$theta_i

  # Determine burnin, defaulting to half of the columns if not provided
  if (is.null(burnin)) {
    burnin <- floor(dim(theta_i_mat)[2] / 2)
  }

  # Remove duplicates and apply burnin
  cat("Processing theta matrix and TR names...\n")
  theta_i_mat <- theta_i_mat[!duplicated(TR_names), (dim(theta_i_mat)[2] - burnin):dim(theta_i_mat)[2]]
  TR_names <- TR_names[!duplicated(TR_names)]

  # Calculate row means and BIT scores
  tr_results_i <- rowMeans(theta_i_mat)
  BIT_score <- logistic(tr_results_i)

  # Calculate confidence intervals
  CI_intervals <- t(apply(theta_i_mat, 1, function(x) quantile(x, probs = c(0.025, 0.975))))

  # Create results dataframe
  cat("Compiling results...\n")
  results_theta_i <- data.frame(
    TR = TR_names,
    Theta_i = tr_results_i,
    lower = CI_intervals[, 1],
    upper = CI_intervals[, 2],
    BIT_score = BIT_score,
    BIT_score_lower = logistic(CI_intervals[, 1]),
    BIT_score_upper = logistic(CI_intervals[, 2])
  )

  # Sort by Theta_i and assign rank
  results_theta_i <- results_theta_i[order(-results_theta_i$Theta_i), ]
  results_theta_i$Rank <- rank(-results_theta_i$Theta_i)

  # Remove row names
  row.names(results_theta_i) <- NULL

  # Write results to CSV
  output_file <- file.path(tools::file_path_as_absolute(output_path),
                           paste0(tools::file_path_sans_ext(basename(file_path)), "_rank_table.csv"))
  write.csv(results_theta_i, output_file, row.names = FALSE)

  cat(paste0("Results saved to ", output_file, "\n"))

  # Return the results
  return(results_theta_i)
}


#' rank_plot
#' @description To draw a barplot for the top n TRs.
#' @param file_path path to the saved BIT Gibbs sampling results.
#' @param output_path path to save the barplot.
#' @param burnin number of samples used for burn-in. If not specify, BIT will use the half of the iterations as burn in.
#' @param n top n TRs will show in the barplot, default: 10.
#' @param colors colors for each bar, default "NPG" for n<=10, has to be manually specified if n>10.
#' @param main main title for the barplot, default: NULL.
#' @param xlab x axis label, default: BIT score.
#' @param ylab y axis label, default: TR symbols.
#' @return a data.frame object contains TR names, theta_i, BIT scores and 95 CIs for each TR.
#' @export
rank_plot <- function(file_path = NULL, output_path = NULL, burnin = NULL, n = 10,
                      colors = "NPG", main = NULL, xlab = "BIT score", ylab = "TR symbols") {

  # Load the dataset from the file path
  dat <- readRDS(file_path)

  TR_names <- dat[["TR_names"]]
  theta_i_mat <- dat$theta_i

  # Set default burn-in period if not provided (half the number of columns)
  if (is.null(burnin)) {
    burnin <- dim(theta_i_mat)[2] %/% 2
  }

  # Remove duplicate TR names and apply burn-in to theta_i_mat
  theta_i_mat <- theta_i_mat[!duplicated(TR_names), (dim(theta_i_mat)[2] - burnin):dim(theta_i_mat)[2]]
  TR_names <- TR_names[!duplicated(TR_names)]

  # Calculate row means (average theta values) after burn-in
  tr_results_i <- rowMeans(theta_i_mat)
  BIT_score <- logistic(tr_results_i)

  # Calculate confidence intervals (2.5% and 97.5% quantiles)
  CI_intervals <- t(apply(theta_i_mat, 1, function(x) quantile(x, probs = c(0.025, 0.975))))

  # Create a data frame of results, including confidence intervals and BIT scores
  results_theta_i <- data.frame(
    TR = TR_names,
    Theta_i = tr_results_i,
    lower = CI_intervals[, 1],
    upper = CI_intervals[, 2],
    BIT_score = BIT_score,
    BIT_score_lower = logistic(CI_intervals[, 1]),
    BIT_score_upper = logistic(CI_intervals[, 2])
  )

  # Order results by Theta_i (descending)
  results_theta_i <- results_theta_i[order(-results_theta_i$Theta_i), ]
  results_theta_i$Rank <- rank(-results_theta_i$Theta_i)
  row.names(results_theta_i) <- NULL

  # Set output file path for the PDF
  output_file_path_complete <- paste0(
    tools::file_path_as_absolute(output_path), "/",
    tools::file_path_sans_ext(basename(file_path)), ".pdf"
  )

  # Select the top n TR names and their corresponding BIT scores
  TR_names_used <- results_theta_i[n:1, "TR"]
  BIT_score_used <- results_theta_i[n:1, "BIT_score"]
  BIT_score_lower_used <- results_theta_i[n:1, "BIT_score_lower"]
  BIT_score_upper_used <- results_theta_i[n:1, "BIT_score_upper"]

  # Set colors based on the provided color palette
  if (colors == "NPG") {
    colors <- ggsci::pal_npg()(n)
  }

  # Generate the bar plot and save it to the PDF
  pdf(output_file_path_complete)
  par(mfrow = c(1, 1), oma = c(1, 1, 1, 1) + 0.1, mar = c(4, 6.5, 2, 1) + 0.1, mgp = c(3, 0.5, 0), font.axis = 2)

  bp <- barplot(BIT_score_used, horiz = TRUE, yaxt = "n",
                xlim = c(0, max(BIT_score_upper_used) * 1.05),
                cex.axis = 1.5, col = colors, tcl = -0.2,
                main = main, cex.main = 1.3, font.axis = 1)

  # Add text labels and points to the bar plot
  text(BIT_score_used - max(BIT_score_used) / 10, bp, round(BIT_score_used, 3), font = 1, cex = 1.2)
  points(BIT_score_used, bp, pch = 16)

  # Draw error bars for the confidence intervals
  segments(BIT_score_used, bp, BIT_score_upper_used, bp, lwd = 2)
  segwidth <- (bp[2] - bp[1]) / 4
  segments(BIT_score_upper_used, bp + segwidth, BIT_score_upper_used, bp - segwidth, lwd = 2)

  # Customize y-axis labels (TR names)
  axis(side = 2, las = 1, at = bp, labels = TR_names_used, tcl = -0.2, cex.axis = 1.2, font.axis = 1)

  # Add axis titles and box
  title(xlab = xlab, line = 3, cex.lab = 1.4, font.lab = 2)
  title(ylab = ylab, line = 5.5, cex.lab = 1.4, font.lab = 2)
  title(main = main, cex.main = 1.4, font.main = 2)
  box()

  # Close the PDF device
  dev.off()

  return()
}


#' compare_scatter_plot
#' @description To draw a barplot for the top n TRs.
#' @param file1_path path to the saved BIT Gibbs sampling results of input 1.
#' @param file2_path path to the saved BIT Gibbs sampling results of input 2.
#' @param output_path path to save the barplot.
#' @param burnin number of samples used for burn-in. If not specify, BIT will use the half of the iterations as burn in.
#' @return NULL
#' @export
compare_scatter_plot<-function(file1_path, file2_path, output_path, burnin=NULL){
  dat1<-readRDS(file1_path)
  dat2<-readRDS(file2_path)

  TR_names_dat1 <- dat1[["TR_names"]]
  dat1_theta_i_mat<-dat1$theta_i

  TR_names_dat2 <- dat2[["TR_names"]]
  dat2_theta_i_mat<-dat2$theta_i

  if(is.null(burnin)){
    burnin=dim(dat1_theta_i_mat)[2]%/%2
  }

  output_file_path_complete<-paste0(tools::file_path_as_absolute(output_path),"/",
                                    tools::file_path_sans_ext(basename(file1_path)),"_",
                                    tools::file_path_sans_ext(basename(file2_path)),"_compare.pdf")

  theta_i_mat_dat1<-dat1_theta_i_mat[!duplicated(TR_names_dat1),(dim(dat1_theta_i_mat)[2]-burnin):dim(dat1_theta_i_mat)[2]]
  TR_names_dat1<-TR_names_dat1[!duplicated(TR_names_dat1)]
  tr_results_i_dat1<-rowMeans(theta_i_mat_dat1)
  BIT_score<-logistic(tr_results_i_dat1)
  CI_intervals<-t(apply(theta_i_mat_dat1, 1, function(x) quantile(x, probs = c(0.025, 0.975))))

  results_theta_i_dat1=data.frame(TR=TR_names_dat1,Theta_i=tr_results_i_dat1,BIT_score=BIT_score,BIT_score_lower=CI_intervals[,1],BIT_score_upper=CI_intervals[,2])
  results_theta_i_dat1=results_theta_i_dat1[order(-results_theta_i_dat1$Theta_i),]
  results_theta_i_dat1$Rank=rank(-results_theta_i_dat1$Theta_i)
  row.names(results_theta_i_dat1) <- NULL

  theta_i_mat_dat2<-dat2_theta_i_mat[!duplicated(TR_names_dat2),(dim(dat2_theta_i_mat)[2]-burnin):dim(dat2_theta_i_mat)[2]]
  TR_names_dat2<-TR_names_dat2[!duplicated(TR_names_dat2)]
  tr_results_i_dat2<-rowMeans(theta_i_mat_dat2)
  BIT_score<-logistic(tr_results_i_dat2)
  CI_intervals<-t(apply(theta_i_mat_dat2, 1, function(x) quantile(x, probs = c(0.025, 0.975))))

  results_theta_i_dat2=data.frame(TR=TR_names_dat2,Theta_i=tr_results_i_dat2,BIT_score=BIT_score,BIT_score_lower=CI_intervals[,1],BIT_score_upper=CI_intervals[,2])
  results_theta_i_dat2=results_theta_i_dat2[order(-results_theta_i_dat2$Theta_i),]
  results_theta_i_dat2$Rank=rank(-results_theta_i_dat2$Theta_i)
  row.names(results_theta_i_dat2) <- NULL

  file1_values<-results_theta_i_dat1[,"BIT_score"]
  file2_values<-results_theta_i_dat2[match(results_theta_i_dat1[,"TR"],results_theta_i_dat2[,"TR"]),"BIT_score"]

  aligned_tables<-data.frame(TR=results_theta_i_dat1[,"TR"],file1_values=file1_values,file2_values=file2_values)

  TR_names1<-results_theta_i_dat1[1:10,"TR"]
  TR_names2<-results_theta_i_dat2[1:10,"TR"]
  TR_union<-union(TR_names1,TR_names2)

  file1_coords<-aligned_tables[which(aligned_tables[,"TR"]%in%TR_union),"file1_values"]
  file2_coords<-aligned_tables[which(aligned_tables[,"TR"]%in%TR_union),"file2_values"]
  TR_plot_names<-aligned_tables[which(aligned_tables[,"TR"]%in%TR_union),"TR"]

  scale<-c()
  scale1_vec<-(file1_values/max(file1_values))
  scale2_vec<-(file2_values/max(file2_values))
  for(i in 1:length(file1_values)){
    scale<-c(scale,max(c(scale1_vec[i],scale2_vec[i])))
  }

  size<-1-median(scale)+scale

  colors<-rgb(file1_values/max(file1_values),0.4,file2_values/max(file2_values),scale)

  pdf(output_file_path_complete)
  par(mfrow=c(1,1),oma = c(1,1,1,1) + 0.1,mar = c(4,4,2,1) + 0.1,mgp=c(3, 0.5, 0),font.axis=2)
  plot(file1_values,file2_values,type="p",pch=19,col=colors,cex=size, xlab="Region set 1 BIT score",ylab="Region set 2 BIT score",font.lab=2,cex.lab=1.4)
  basicPlotteR::addTextLabels(file1_coords, file2_coords, TR_plot_names, cex.label=1.4,lwd=2)
  box()
  dev.off()

  return()
}


#' Add non-overlapping text labels to plot
#' The function addTextLabels is copied and used under GPL 3.0 license,
#' we thank the contribution by original author: Joseph Crispell https://github.com/JosephCrispell/basicPlotteR
#' @param xCoords A vector containing the X coordinates for labels
#' @param yCoords A vector containing the Y coordinates for labels
#' @param labels A vector containing the labels to be plotted
#' @param cex.label A number to scale the size of the plotted labels. Defaults to 1
#' @param col.label The colour of the plotted labels. Defaults to "red". Multiple colours can be provided. If more colours than labels provided colours will be recycled.
#' @param col.line The colour of the line to plot from relocated labels to original location. Defaults to "black". Multiple colours can be provided. If more colours than labels provided colours will be recycled.
#' @param col.background An optional colour for a background polygon plotted behind labels. Defaults to NULL - won't be plotted. Multiple colours can be provided. If more colours than labels provided colours will be recycled.
#' @param lty A number detailing the type of line to plot from relocated labels to original location. 0: blank, 1: solid, 2: dashed, 3: dotted, 4: dotdash, 5: longdash, and 6: twodash. Defaults to 1. Multiple line types can be provided. If more options than labels provided types will be recycled.
#' @param lwd A number to scale the size of line from relocated labels to original location. Defaults to 1. Multiple line widths can be provided. If more options than labels provided widths will be recycled.
#' @param border The colour of the border to be plotted around the polygon. Defaults to NA - won't be plotted. Multiple colours can be provided. If more colours than labels provided colours will be recycled.
#' @param avoidPoints A logical variable indicating whether labels shouldn't be plotted on top of points. Defaults to TRUE
#' @param keepLabelsInside A logical variable indicating whether the labels shouldn't be plotted outside of plotting region. Defaults to TRUE
#' @param cex.pt A number used to scale the points plotted on the graph that labels are to be added to. Defaults to 1
#' @keywords text label plot
#' @examples
#' # Create some random points
#' n <- 50
#' coords <- data.frame(X=runif(n), Y=runif(n), Name="Test Label")
#'
#' # Plot them without labels
#' plot(x=coords$X, y=coords$Y, pch=19, bty="n", xaxt="n", yaxt="n", col="red", xlab="X", ylab="Y")
#'
#' # With potentially overlapping labels
#' plot(x=coords$X, y=coords$Y, pch=19, bty="n", xaxt="n", yaxt="n", col="red", xlab="X", ylab="Y")
#' text(coords$X, coords$Y, labels=coords$Name, xpd=TRUE)
#'
#' # Plot them with non-overlapping labels
#' plot(x=coords$X, y=coords$Y, pch=19, bty="n", xaxt="n", yaxt="n", col="red", xlab="X", ylab="Y")
#' addTextLabels(coords$X, coords$Y, coords$Name, cex.label=1, col.label="black")
#'
#' # Plot them with non-overlapping labels
#' plot(x=coords$X, y=coords$Y, pch=19, bty="n", xaxt="n", yaxt="n", col="red", xlab="X", ylab="Y")
#' addTextLabels(coords$X, coords$Y, coords$Name, cex.label=1, col.background=rgb(0,0,0, 0.75), col.label="white")
addTextLabels <- function(xCoords, yCoords, labels, cex.label=1, col.label="red", col.line="black", col.background=NULL,
                          lty=1, lwd=1, border=NA, avoidPoints=TRUE, keepLabelsInside=TRUE, cex.pt=1){

  #######################################################
  # Check that the input data are in the correct format #
  #######################################################

  # Are each of coordinate vectors the same length?
  if(length(xCoords) != length(yCoords)){
    stop("addTextLabels() The vectors containing the X and Y coodinates must be the same length.")
  }
  if(length(xCoords) != length(labels)){
    stop("addTextLabels() The vector of labels must be the same length as the coordinate vectors.")
  }

  #######################
  # Get the axis limits #
  #######################

  # Get the axis limits
  axisLimits <- graphics::par("usr")

  ############################
  # Check for NA coordinates #
  ############################

  # Check if any NA coordinates present
  indicesOfNAs <- which(is.na(xCoords) | is.na(yCoords))
  if(length(indicesOfNAs) > 0){

    # Send warning
    warning("NA values present in coordinates provided. These are ignored.")

    # Check for each of the parameters that can have multiple parameters
    if(length(col.line) == length(xCoords)){
      col.line = col.line[-indicesOfNAs]
    }
    if(length(col.background) == length(xCoords)){
      col.background = col.background[-indicesOfNAs]
    }
    if(length(lty) == length(xCoords)){
      lty = lty[-indicesOfNAs]
    }
    if(length(lwd) == length(xCoords)){
      lwd = lwd[-indicesOfNAs]
    }
    if(length(border) == length(xCoords)){
      border = border[-indicesOfNAs]
    }

    # Remove the NA coordinates
    xCoords <- xCoords[-indicesOfNAs]
    yCoords <- yCoords[-indicesOfNAs]

    # Remove the respective labels
    labels <- labels[-indicesOfNAs]
  }

  ############################
  # Check if axes are logged #
  ############################

  # Check X axis
  xAxisLogged <- FALSE
  if(graphics::par("xlog")){

    # Note that X axis was logged
    xAxisLogged <- TRUE

    # Log the X coordinates
    xCoords <- log10(xCoords)

    # Reset the X axis logged flag - fools points and polygon commands below
    graphics::par(xlog=FALSE)
  }

  # Check Y axis
  yAxisLogged <- FALSE
  if(graphics::par("ylog")){

    # Note that Y axis was logged
    yAxisLogged <- TRUE

    # Log the Y coordinates
    yCoords <- log10(yCoords)

    # Reset the Y axis logged flag - fools points and polygon commands below
    graphics::par(ylog=FALSE)
  }

  ###############################
  # Store the point information #
  ###############################

  # Store the input coordinates and labels
  pointInfo <- list("X"=xCoords, "Y"=yCoords, "Labels"=labels, "N"=length(xCoords), "cex"=cex.pt)

  # Set the amount to pad onto height and width
  heightPad <- 0.5
  widthPad <- 0.04
  if(is.null(col.background)){
    heightPad <- 0
    widthPad <- 0
  }

  # Calculate the label heights and widths
  pointInfo <- calculateLabelHeightsAndWidths(pointInfo=pointInfo, cex=cex.label,
                                              heightPad=heightPad, widthPad=widthPad)

  ###########################################
  # Produce a list of alternative locations #
  ###########################################

  # Generate the alternative locations
  alternativeLocations <- generateAlternativeLocations(axisLimits)

  # Calculate the distance between the actual and alternative points - rescale X axis remove axis range bias
  distances <- euclideanDistancesWithRescaledXAxis(pointInfo, alternativeLocations, axisLimits)

  ###############################################################
  # Create a list to store the information about plotted points #
  ###############################################################

  # Initialise the list to store the information about plotted labels
  plottedLabelInfo <- list("X"=c(), "Y"=c(), "Height"=c(), "Width"=c(), "N"=0)

  ##############################################################
  # Add labels to plot assigning new locations where necessary #
  ##############################################################

  # Plot the point label
  for(i in seq_len(pointInfo$N)){

    # Set the colours for plotting the label - allows multiple colours and cycling through colours
    labelColour <- setOption(options=col.label, index=i)
    backgroundColour <- setOption(options=col.background, index=i)
    borderColour <- setOption(options=border, index=i)

    # Set the line characteristics
    lineColour <- setOption(options=col.line, index=i)
    lineType <- setOption(options=lty, index=i)
    lineWidth <- setOption(options=lwd, index=i)

    # Get the information for the current point
    x <- pointInfo$X[i]
    y <- pointInfo$Y[i]
    label <- pointInfo$Labels[i]
    height <- pointInfo$Heights[i]
    width <- pointInfo$Widths[i]

    # Get a new location
    newLocationIndex <- chooseNewLocation(pointInfo, i, alternativeLocations, distances, plottedLabelInfo, axisLimits, keepLabelsInside)

    # Is the current point too close to others?
    if(alternativeLocations$N != 0 && newLocationIndex != -1 &&
       (avoidPoints == TRUE || tooClose(x, y, height, width, plottedLabelInfo) || outsidePlot(x, y, height, width, axisLimits))){

      # Get the coordinates for the chosen alternate location
      altX <- alternativeLocations$X[newLocationIndex]
      altY <- alternativeLocations$Y[newLocationIndex]

      # Add line back to previous location
      addLineBackToOriginalLocation(altX=altX, altY=altY, x=x, y=y, label=label,
                                    cex=cex.label, col=lineColour, lty=lineType, lwd=lineWidth, heightPad=heightPad,
                                    widthPad=widthPad)

      # Add label
      addLabel(x=altX, y=altY, label=label,
               cex=cex.label, col=labelColour, bg=backgroundColour, border=borderColour, heightPad=heightPad, widthPad=widthPad)

      # Append the plotted label information
      plottedLabelInfo <- addPlottedLabel(x=altX, y=altY, height=height, width=width,
                                          plottedLabelInfo=plottedLabelInfo)

      # Remove the alternative plotting location used
      alternativeLocations$X <- alternativeLocations$X[-newLocationIndex]
      alternativeLocations$Y <- alternativeLocations$Y[-newLocationIndex]
      alternativeLocations$N <- alternativeLocations$N - 1
      distances <- distances[, -newLocationIndex]

    }else{

      # Add label
      addLabel(x=x, y=y, label=label,
               cex=cex.label, col=labelColour, bg=backgroundColour, border=borderColour,
               heightPad=heightPad, widthPad=widthPad)

      # Append the plotted label information
      plottedLabelInfo <- addPlottedLabel(x=x, y=y, height=height, width=width,
                                          plottedLabelInfo=plottedLabelInfo)
    }
  }

  #####################################################################################
  # Return axes logged flags to original state - for if person makes any future plots #
  #####################################################################################

  graphics::par(xlog=xAxisLogged)
  graphics::par(ylog=yAxisLogged)

}

#' function used to compute the logistic transformation
#' @param x vectorized object
#'
#' @return None
#'
logistic <- function(x) {
  1 / (1 + exp(-x))
}
