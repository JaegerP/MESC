######################################################
#         0. initialization                          #
######################################################

## Edit these three parameters to the location of code / image file
codepath <- "C:/Users/User/Philipp_Fuhrmann/roughness"
filepath <- "C:/Users/User/Philipp_Fuhrmann/scan_data/11_06"
filename <- "4_Mgo_Fe 2012_11_06007 [Z_fwd] image_inset2.bmp"
#gold sample 2012_10_30004 [Z_fwd].asc


## edit distance and height according to the scale of the image
distance <- 2500/1500 	# distance between two scan points in nm (1000/350 means 350 pixel equal 1 micron)
height <- 12.344/256 		# height conversion factor from color scale to nm (256 is total color depth)
						# local_height(x,y) = (color_value(x,y) / 256) * height_in_nm
						
## edit these parameters to optimize the fitting functions
x_search_correlation <- 400  	#	approximation for the correlation length in pixel
x_start_fit <- 3          		# 1 or higher, leave out the first couple of pixels for the roughness fit (sometimes improves results a lot)
x_drop_last_pixels <- 400   	# number of pixels that are left out at the end for calculation of interface layer thickness


search_correlation <- 400  	#	approximation for the correlation length in pixel
start_fit <- 3          		# 1 or higher, leave out the first couple of pixels for the roughness fit (sometimes improves results a lot)
drop_last_pixels <- 400   	# number of pixels that are left out at the end for calculation of interface layer thickness

## initialize variables
direction <- "x"
Nx <- 0 				# number of values along the fast scan axis
Ny <- 0 				# number of values along the slow scan axis
a <- 0 	 	    		# roughness
w <- 0	 	  			# interface thickness
xi <- 0   				# correlation length
avg_height <- 0 		# average height

## read the function that imports images
setwd(codepath)
source("read-bmp.R")
setwd(filepath)


#file_parameters <- function() {
#	rel_filepath <- readline(prompt = "relative file path: (leave empty for no change)")
#	setwd(paste(filename,"/", rel_filepath,sep=""))
#	filename <- readline(prompt = "file name:(leave empty for no change)")
#	distance <- readline(prompt = "Distance between two neighbouring points in nm:(leave empty for no change)")
#	height <- readline(prompt = "Height between two color levels in nm:(leave empty for no change)")
#	new_paramters <- list(rel_filepath, filename, distance, height)
#	return (new_parameters)
#}


#parameters <- function() {
#	search_correlation <- readline(prompt = "Approximate upper limit for correlation length in pixels:")
#	start_fit <- readline(prompt = "Omit first # points for roughness fit: ")
#	drop_last_pixles <- readline(prompt = "Omit last # points for interface layer thickness fit")
#}

######################################################
#   1. determin the height correlation function      #
######################################################
#
# see phd thesis of tim mewes, anhang a, for further description of the method
#
## convert color scale to height

loadfile <- function(option) {

	############ initialize
	if (option == 0) { ## if the file is a bmp
		image_data <- read.bmp(filename) # pixel data of the STM-Scan:
										# 3 dimensions [x, y,(r, g, b)] r,g,b are values from 0 to 255 indicating red, green and blue
										# that are always equal for an individual pixel (gray scale)
		n <- dim(image_data) 			 # array with 3 elements: (x-length, y-length, 3), sometimes 2 elements (x-length, y-length)
		if (length(n) == 3) { 
			sample_matrix <- t(image_data[,,1]) # 2D array with one number from 0 to 255 for each pixel
		} else {
			sample_matrix <- t(as.matrix(image_data))
		}
		Nx <- n[2] 				# x-length
		Ny <- n[1] 				# y-length
		Hexp_x <- rep(0, Nx)    # height correlation function in x direction
		Hexp_y <- rep(0, Ny)    # height correlation function in y direction
		print(n[1])
		print(n[2])
		########## convert height information

		sample_matrix <- sample_matrix * height
	} else {
		image_data <- read.table(filename, skip = 33)
		Nx <- dim(image_data[1])
		Ny <- dim(image_data[1])
	}
	## get average height
	avg_height <- mean(sample_matrix)
	peak_to_valley <- max(sample_matrix) - min(sample_matrix)
	rms_roughness <- sqrt(1/(Nx^2) * sum(sample_matrix^2))

	## calculate height correlation function in x direction
	for (i in 1:round((Nx-1)/2)) {	
			for (k in 1:(Nx-i)) {
					Hexp_x[i] <- Hexp_x[i] + sum((sample_matrix[i+k,] - sample_matrix[k,])^2)
			}
		# normalize
		Hexp_x[i] <- Hexp_x[i] / (Ny * (Nx - i))
		print(paste("Point", i , "of", Nx))
	}
	
	## calculate height correlation function in y direction
	for (i in 1:round((Ny-1)/2)) {		
			for (k in 1:(Ny-i)) {
					Hexp_y[i] <- Hexp_y[i] + sum((sample_matrix[,i+k] - sample_matrix[,k])^2)
			}
		# normalize
		Hexp_y[i] <- Hexp_y[i] / (Nx * (Ny - i))
		print(paste("Point", i , "of", Ny))
	}
	returnlist <- list(Hexp_x, Hexp_y, avg_height, peak_to_valley, rms_roughness, Nx, Ny, sample_matrix)
	return(returnlist)
}

######################################################
#   2. calculate characteristic parameters           #
######################################################
#
# calculate the roughness a, interface layer thickness w, and correlation length xi
# a is calculated by determining the slope of a linear fit to the double logarithmic 
# representation of the height correlation function for small x
#
# w is the saturation value for the height correlation function is determined by fitting 
# the height correlation function to a constant for large x
#
# xi is then calculated using the y-intercept of the first graph.
# for small x the height correlation function is proportional to 2*w^2 * (x/xi)^2a
# the y-intercept in the double logarithmic plot is thus 2*w^2/xi^2a and xi can be determined.

determine_err <- function(Hexp, search_correlation, w, xi, a) {
	#print(paste(a, w, xi, search_correlation, distance))
	#print(Hexp)
	x <- 1:search_correlation * distance
	#err <- sum((log(Hexp[1:search_correlation]) - log(w*(1-exp(-(x/xi)^(2*a)))))^2*t(1/1:search_correlation))
	err <- sum((log(Hexp[1:search_correlation]) - log(w*(1-exp(-(x/xi)^(2*a)))))^2)
	#if (err < 0.02) {
	#	print(paste(w, xi, a, err))
	#}
	#print(err)
	return(err)
}

mod <- function( a, b) {
	d <- floor(a/b)
	a <- a- b * d
}

opti <- function(Hexp, a, w, xi, precision, search_correlation) {
	## precision: 1 or 0.1 or 0.01
	print(paste("a",a, "w",w, "xi",xi))
	a_array <- a - 0.1*precision + 0.01*precision*0:20
	w_array <- w - 1*precision + 0.1*precision*0:20
	xi_array <- xi - 1*precision + 0.1*precision*0:20
	
	err_array <- rep(0,21^3)
	dim(err_array) <- c(21, 21, 21)
	for (i in 1:21) {
		for(j in 1:21) {
			for(k in 1:21) {
				err_array[i, j, k] <- determine_err(Hexp, search_correlation, w_array[i], xi_array[j], a_array[k])
			}
		}
	}
	total <- which.min(err_array)
	print(which.min(err_array))
	first <- mod(total, 21)
	total <- (total - first)/21
	second <- mod(total, 21) + 1
	third <- (total - second + 1) / 21 + 1
	print(min(err_array))
	print(paste(first, second, third))
	a_local <- a - 0.1 * precision + (third-1) * 0.01 * precision
	xi_local <- xi - 1 * precision + (second-1) * 0.1 * precision
	w_local <- w - 1 * precision + (first-1) * 0.1 * precision
	print(paste(a_local, w_local, xi_local))
	return(c(a_local, w_local, xi_local))
}

analysis <- function(direction, sav=FALSE){
	if (direction == "x") {
		N <- Nx
		Hexp <- Hexp_x
		if (sav) {
			png(paste(filename, "_x_analysis.png", sep=""), width=1024, height=768)
		}
	}
	if (direction == "y") {
		N <- Ny
		Hexp <- Hexp_y
		if (sav) {
			png(paste(filename, "_y_analysis.png", sep=""), width=1024, height=768)
		}
	}
	if (direction != "y" && direction != "x") {
		return
	}
	
	plot(x = 1:N, Hexp, log="xy", main = paste(direction, "-Direction-Height-Correlation-Function", sep=""), pch = 46, xlab = paste(direction, "[nm]"), ylab = "<z(r)-z(0)>^2[nm^2]", cex = 5, lab = c(20,5,5), xaxp=c(1,1000,3))
	search_correlation <- as.integer(readline("approximation for the correlation length in pixel: search_correlation = "))
	uncertainty <- as.integer(readline("uncertainty of this approximation (about 100): uncertainty = "))
	start_fit <- as.integer(readline("1 or higher, leave out the first couple of pixels for the roughness fit (sometimes improves results a lot): start_fit = "))
	#drop_last_pixels <- as.integer(readline("number of pixels that are left out at the end for calculation of interface layer thickness: drop_last_pixels = "))

	## get roughness a
	
	# plot(x = 1:Ny, Hexp_y, log="xy") ## x axis in pixels
	old = Inf
	new = 1000
	count <- -1
	
	## find the best fit (optimization of mean square residuals)
	residual_array <- rep(0, search_correlation)
	while(search_correlation - count > 3 + start_fit) {
		count <- count + 1
		x <- log(start_fit:(search_correlation - count) * distance)
		y <- log(Hexp[start_fit:(search_correlation - count)])
		fit <-  lm(y ~ x)
		residual_array[search_correlation - count] <- mean(fit[[2]]^2)
		#print(new)
	}
	count <- search_correlation - which.min(residual_array[(3 + start_fit):search_correlation])
	# BACKUP while(old > new) {
	#	count <- count + 1
	#	x <- log(start_fit:(search_correlation - count) * distance)
	#	y <- log(Hexp[start_fit:(search_correlation - count)])
	#	fit <-  lm(y ~ x)
	#	old <- new
	#	new <- mean(fit[[2]]^2)
	#	#print(new)
	#}
	count <- count - 1
	x <- log(start_fit:(search_correlation - count) * distance)
	y <- log(Hexp[start_fit:(search_correlation - count)])
	fit <-  lm(y ~ x)
	#plot(x, y) #plot with logarithmic data and linear axis
	#curve(fit[[1]][1] + fit[[1]][2]*x, add=TRUE) # fit curve with logarithmic data and linear axis
	aoffset <- fit[[1]][1] 
	
	a <- fit[[1]][2]/2 # store roughness value
	xi <- which.max(Hexp[(search_correlation - uncertainty):(search_correlation + uncertainty)]) + search_correlation - uncertainty # wrong logic. should be: where linear fit hits interlayer thickness
	w <- Hexp[xi]

	
	xi <- exp((log(w)-aoffset)/(2*a))
	
	curve(w*(1-exp(-(x/xi)^(2*a))), add=TRUE, col = "pink", lwd = 2, sub = paste(direction, "- Direction: roughness:", a, "interface layer thickness:", sqrt(w/2), "correlation length:", xi)) # fit curve
	#print(paste("curve with a: ", a , "w:", w, "xi:", xi))
	print(paste(a, w, xi))
	## optimize fit parameters
	temp <- c(a,w,xi)
	if (a > 0.1 && w > 1) {temp <- opti(Hexp, temp[1], temp[2], temp[3],1, search_correlation)}
	temp <- opti(Hexp, temp[1], temp[2], temp[3],0.1, search_correlation)
	temp <- opti(Hexp, temp[1], temp[2], temp[3],0.01, search_correlation)
	a <- temp[1]
	w <- temp[2]
	xi <- temp[3]
	
	plot(x = 1:round(N/2) * distance, Hexp[1:round(N/2)], log="xy", main = paste(direction, "-Direction-Height-Correlation-Function", sep=""), sub = paste("roughness:", a, "interface layer thickness:", sqrt(w/2), "correlation length:", xi), pch = 46, xlab = paste(direction, "[nm]"), ylab = "<z(r)-z(0)>^2[nm^2]", cex = 4)
	curve(exp(aoffset + a*2*log(x)), add=TRUE, col = "red", lwd = 2) # fit curve
	curve(w + x * 0, add=TRUE, col = "blue", lwd = 2) # fit curve
	curve(w*(1-exp(-(x/xi)^(2*a))), add=TRUE, col = "green", lwd = 2) # fit curve
	print(paste("curve with a: ", a , "w:", w, "xi:", xi))

	## store x and y coordinates for the calculation of xi later
#	x1 <- x
#	y1 <- y
#
#	## get interface layer thickness w
#	old = Inf
#	new = 1000
#	count <- search_correlation
#	## find the best fit (optimization of mean square residuals)
#	while(old > new) {
#		count <- count + 1
#		x <- (count):(N - drop_last_pixels)
#		y <- Hexp[(count):(N - drop_last_pixels)]
#		fit <-  lm(y ~ 1) # fit with linear regression
#		old <- new
#		new <- mean(fit[[2]]^2) # get mean residuals
#		#print(new)
#	}
#	count <- count - 1
#	x <- (count):(N - drop_last_pixels)
#	y <- Hexp[(count):(N - drop_last_pixels)]
#	fit <-  lm(y ~  1) # fit with linear regression
#	
#	w <- sqrt(fit[[1]][1]/2) # store interface layer thickness
#
#
#	## get correlation length
#	fit <- lm(y1~x1) # fit with linear regression
#	xi <-(2*w^2/exp(fit[[1]][1])) ^(1/(2*a)) # store correlation length
#	curve(2*w^2*(1-exp(-(x/xi)^(2*a))), add=TRUE, col = "green", lwd = 2, sub = paste(direction, "- Direction: roughness:", a, "interface layer thickness:", w, "correlation length:", xi)) # fit curve
#	## if this last curve is off, it is usually to a bad fit for the roughness. has to be manually adjusted by deselecting a few points.
#
#	## Output
	print(paste("Results: Fit using the following fit modifiers: search_correlation:", search_correlation, "start_fit:", start_fit, "drop_last_pixels:", drop_last_pixels, "uncertainty:", uncertainty))
	print(paste(direction, "- Direction: roughness:", a, "interface layer thickness:", sqrt(w/2), "correlation length:", xi))
	print(paste("average height:", avg_height, "peak to valley:", peak_to_valley, "root-mean-square roughness:", rms_roughness))
	if (sav) {
		dev.off()
		## save textfile of settings
		savelist <- list(search_correlation, start_fit, drop_last_pixels, distance, height, avg_height, peak_to_valley, rms_roughness)
		cat("search_correlation: ", search_correlation, "start_fit: ", start_fit, "drop_last_pixels: ", drop_last_pixels, "distance: ", distance, "height: ", height, "avg_height: ", avg_height, "peak_to_valley: ", peak_to_valley, "rms_roughness: ", rms_roughness, file = paste(filename, "_info.txt", sep=""))
	}
}

######################################################
#   3. plot surface                                  #
######################################################
surface <- function() {
	persp(seq(from = 1, to = Nx, by = 10), seq(from = 1, to = Ny, by = 10), sample_matrix[seq(from = 1, to = Nx, by = 10),seq(from = 1, to = Ny, by = 10)], col="lightgreen",
		theta=40, phi=30, r=100, d=0.1, expand=0.5, ltheta=90, lphi=180, shade=0.2, ticktype="detailed", nticks=5, xlab = "x[nm]", ylab = "y[nm]", zlab = "z[nm]")
}

######################################################
#   3. Execute                                       #
######################################################

returnlist <- loadfile(0)
Hexp_x <- returnlist[[1]]
Hexp_y <- returnlist[[2]]
avg_height <- returnlist[[3]]
peak_to_valley <- returnlist[[4]]
rms_roughness <- returnlist[[5]]
Nx <- returnlist[[6]]
Ny <- returnlist[[7]]
sample_matrix <- returnlist[[8]]

######################################################
#   3. Save  to file                                 #
######################################################
saveall <- function() {
	## save x plot
	png(paste(filename, "_x_analysis.png", sep=""))
	analysis("x")
	dev.off()

	## save y plot
	png(paste(filename, "_y_analysis.png", sep=""))
	analysis("y")
	dev.off()

	## save 3D surface plot
	png(paste(filename, "_surface.png", sep=""))
	surface()
	dev.off()

	## save textfile of settings
	savelist <- list(search_correlation, start_fit, drop_last_pixels, distance, height, avg_height, peak_to_valley, rms_roughness)
	cat("search_correlation: ", search_correlation, "start_fit: ", start_fit, "drop_last_pixels: ", drop_last_pixels, "distance: ", distance, "height: ", height, "avg_height: ", avg_height, "peak_to_valley: ", peak_to_valley, "rms_roughness: ", rms_roughness, file = paste(filename, "_info.txt", sep=""))
}