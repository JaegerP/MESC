######################################################
#         0. initialization                          #
######################################################

## Edit these three parameters to the location of code / image file
codepath <- "C:/Users/Joe.joelap/Desktop/Google Drive/Uni/roughness"
filepath <- "C:/Users/Joe.joelap/Desktop/Google Drive/Uni/roughness/no5"
filename <- "5_Mgo_Fe_2012_11_07002 [Z_fwd] image.bmp"
#gold sample 2012_10_30004 [Z_fwd].asc


## edit distance and height according to the scale of the image
distance <- 200/800 	*5.1 # distance between two scan points in nm (1000/350 means 350 pixel equal 1 micron). Factor of 5,1 is applied to images taken before the xy calibration on jan. 23 2013
height <-0.61634556  /256 	#*0.1077 # height conversion factor from color scale to nm (256 is total color depth). be careful to use the correct z-calibration. Factor of 0.01077 is applied to images taken before the z calibration in nov 2012
						# local_height(x,y) = (color_value(x,y) / 256) * height_in_nm
						
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
		sample_matrix[,1:Ny] <- sample_matrix[,Ny:1]
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

## Returns a measure for the error made by the fit function relative to the actual measurement points
determine_err <- function(Hexp, search_correlation, w, xi, a, start_fit) {
	x <- start_fit:search_correlation * distance
	#err <- sum((log(Hexp[1:search_correlation]) - log(w*(1-exp(-(x/xi)^(2*a)))))^2*t(1/1:search_correlation))
	err <- sum((log(Hexp[start_fit:search_correlation]) - log(w*(1-exp(-(x/xi)^(2*a)))))^2)#*t(1/sqrt(1:search_correlation))
	#err <- sum((Hexp[1:search_correlation] - w*(1-exp(-(x/xi)^(2*a))))^2*t(1/1:search_correlation))
	#err <- sum((Hexp[1:search_correlation] - w*(1-exp(-(x/xi)^(2*a))))^2)
	return(err)
}

## modulo function
mod <- function( a, b) {
	d <- floor(a/b)
	a <- a- b * d
}

## Varies the fit parameters a, w, xi to minimize the error of the fit function
opti <- function(Hexp, a, w, xi, precision, search_correlation, start_fit) {
	## precision: 1 or 0.1 or 0.01
	print(paste("a",a, "w",w, "xi",xi))
	a_array <- a - 0.1*precision + 0.01*precision*0:20
	w_array <- w - 0.05*precision + 0.005*precision*0:20
	xi_array <- xi - 5*precision + 0.5*precision*0:20
	
	err_array <- rep(0,21^3)
	dim(err_array) <- c(21, 21, 21)
	for (i in 1:21) {
		for(j in 1:21) {
			for(k in 1:21) {
				err_array[i, j, k] <- determine_err(Hexp, search_correlation, w_array[i], xi_array[j], a_array[k], start_fit)
			}
		}
	}
	total <- which.min(err_array)
	#print(which.min(err_array))
	first <- mod(total, 21)
	if(first == 0) {first <- 21}
	total <- (total - first)/21
	second <- mod(total, 21) + 1
	third <- (total - second + 1) / 21 + 1
	print(min(err_array))
	#print(paste(first, second, third))
	a_local <- a - 0.1 * precision + (third-1) * 0.01 * precision
	xi_local <- xi - 5 * precision + (second-1) * 0.5 * precision
	w_local <- w - 0.05 * precision + (first-1) * 0.005 * precision
    #print(paste(a_local, w_local, xi_local))
	#print(determine_err(Hexp, search_correlation, w_local, xi_local, a_local, start_fit))
	return(c(a_local, w_local, xi_local))
}

## Roughness analysis via Height-Height Correlation Function
## this function queries the user for approximate fit data and the varies the parameters for the best possible fit
analysis <- function(direction, sav=FALSE){
	## Initialize
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
	
	## plot initial data and ask user for approximations
	plot(x = 1:N, Hexp, log="xy", main = paste(direction, "-Direction-Height-Correlation-Function", sep=""), pch = 46, xlab = paste(direction, "[nm]"), ylab = "<z(r)-z(0)>^2[nm^2]", cex = 5, lab = c(20,5,5), xaxp=c(1,1000,3))
	search_correlation <- as.integer(readline("approximation for the correlation length in pixel: search_correlation = "))
	start_fit <- as.integer(readline("1 or higher, start of the range of the roughness exponent fit: start_fit = "))
	end_fit <- as.integer(readline("higher than start_fit, end of the range for the roughness exponent fit: end_fit = "))
	approx_w <- as.double(readline("approximate saturation value: approx_w = "))
	
	## get roughness a
	
	# plot(x = 1:Ny, Hexp_y, log="xy") ## x axis in pixels
	old = Inf
	new = 1000
	count <- -1
	
	## find the best fit (optimization of mean square residuals) for roughness exponent
	residual_array <- rep(0, end_fit)
	while(end_fit - count > 3 + start_fit) {
		count <- count + 1
		x <- log(start_fit:(end_fit - count) * distance)
		y <- log(Hexp[start_fit:(end_fit - count)])
		fit <-  lm(y ~ x)
		residual_array[end_fit - count] <- mean(fit[[2]]^2)
		#print(new)
	}
	count <- end_fit - which.min(residual_array[(3 + start_fit):end_fit])

	count <- count - 1
	x <- log(start_fit:(end_fit - count) * distance)
	y <- log(Hexp[start_fit:(end_fit - count)])
	fit <-  lm(y ~ x)
	#plot(x, y) #plot with logarithmic data and linear axis
	#curve(fit[[1]][1] + fit[[1]][2]*x, add=TRUE) # fit curve with logarithmic data and linear axis
	aoffset <- fit[[1]][1] 
	
	## Calculate xi, set w
	approx_a <- fit[[1]][2]/2 # store roughness value
	a <- approx_a
	w <- approx_w	
	approx_xi <- exp((log(w)-aoffset)/(2*a))
	xi <- approx_xi
	
	print(paste("Approximation Values: ", a , "w:", w, "xi:", xi,". Start fit:"))
	
	## optimize fit parameters
	temp <- c(a,w,xi)
	if (a > 0.1 && w > 0.05) {temp <- opti(Hexp, temp[1], temp[2], temp[3],1, search_correlation, start_fit)}
	temp <- opti(Hexp, temp[1], temp[2], temp[3],0.1, search_correlation, start_fit)
	temp <- opti(Hexp, temp[1], temp[2], temp[3],0.01, search_correlation, start_fit)
	a <- temp[1]
	w <- temp[2]
	xi <- temp[3]
	
	## plot
	plot(x = 1:round(N/2) * distance,Hexp[1:round(N/2)], log="xy", main = paste(direction, "-Direction-Height-Correlation-Function", sep=""), sub = paste("roughness:", a, "interface layer thickness:", sqrt(w/2), "correlation length:", xi), pch = 46, xlab = paste(direction, "[nm]"), ylab = "<z(r)-z(0)>^2[nm^2]", cex = 4)
	curve(exp(aoffset + a*2*log(x)), add=TRUE, col = "red", lwd = 2) # fit curve
	print(aoffset)
	curve(w + x * 0, add=TRUE, col = "blue", lwd = 2) # fit curve
	curve(w*(1-exp(-(x/xi)^(2*a))), add=TRUE, col = "green", lwd = 2) # fit curve
	curve(approx_w*(1-exp(-(x/approx_xi)^(2*approx_a))), add=TRUE, col = "pink", lwd = 2) # fit curve
	print(paste("Optimized Values: ", a , "w:", w, "xi:", xi))


	## Output
	print(paste("Results: Fit using the following fit modifiers: search_correlation:", search_correlation, "start_fit:", start_fit, "end_fit:", end_fit))
	print(paste(direction, "- Direction: roughness:", a, "interface layer thickness:", sqrt(w/2), "correlation length:", xi))
	print(paste("average height:", avg_height, "peak to valley:", peak_to_valley, "root-mean-square roughness:", rms_roughness))
	
	## save textfile of settings
	if (sav) {
		dev.off()
		savelist <- list(search_correlation, start_fit, distance, height, avg_height, peak_to_valley, rms_roughness)
		cat(" Sample geometry: \n distance between adjacent pixels[nm]: ", distance, "\n height for color scale [nm]: ", height*256, 
		"\n \n Statistical features: avg_height[nm]: ", avg_height, "\n peak_to_valley[nm]: ", peak_to_valley, "\n rms_roughness[nm]: ", rms_roughness, 
		"\n \n User defined values for the fit routine: \n search_correlation: ", search_correlation, "\n start_fit: ", start_fit, "\n end_fit:  ", end_fit, "\n approximate 2*w^2 (w = interface layer thickness)[nm^2]: ", approx_w, 
		"\n \n Results of the correlation analysis: \n roughness exponent a: ", a, "\n interface layer thickness w[nm]: ", sqrt(w/2), "\n correlation length[nm]: ", xi, file = paste(filename, "_", direction, "_info.txt", sep=""))
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
#   4. cross section                                 #
######################################################

cross_section <- function() {
	query = "n"
	while (query !="y") {
		image(1:Nx, 1:Ny, sample_matrix, col=heat.colors(256))
		dummy <- sample_matrix
		dummy[,] <- NA
		## user input: start and end points for cross section
		start_x <- as.integer(readline("Start_x = "))
		start_y <- as.integer(readline("Start_y = "))
		end_x <- as.integer(readline("End_x = "))
		end_y <- as.integer(readline("End_y = "))
	
		## check if input is valid
		if (!start_x %in% 1:Nx || !start_y %in% 1:Ny ||!end_x %in% 1:Nx ||!end_y %in% 1:Ny) {
			print("Out of bounds.")
			return
		}
	
		## calculate lenght, number of points
		cross_length <- sqrt((end_x - start_x)^2 + (end_y - start_y)^2)*distance
		number_of_points <- max(abs(end_x - start_x), abs(end_y - start_y))
	
		## check if length is zero
		if (number_of_points == 0) {
			print("Number of points is zero.")
			return
		}
	
		## get principle direction (x or y)
		## Assuming the principle direction is x, the algorithm steps either in x-direction or in both x and y direction.
		if (abs(end_x - start_x) > abs(end_y - start_y)) {
			direction = "x"
		} else {
			direction = "y"
		}
	
		## determine direction of scan (x, y) and target ratio (x,y)
		target_xy_ratio = abs(end_x - start_x) / abs(end_y - start_y) ## Inf for x direction, 0 for y direction. Else: move so that ratio is best maintained
		x_dir = ((end_x > start_x)-0.5) *2 # 1 for positive, -1 for negative
		y_dir = ((end_y > start_y)-0.5) *2 # 1 for positive, -1 for negative
		print(target_xy_ratio)
	
		## initialize
		position_array <- rep(0,number_of_points )
		cross_sec <- data.frame(x=position_array, y = position_array, z =position_array)
		current_x <- start_x
		current_y <- start_y
	
		## MAIN LOOP
		if (direction == "x") {
			for (i in 1:number_of_points) {
				if (target_xy_ratio == Inf || abs(((current_x + x_dir-start_x) / (current_y-start_y) - target_xy_ratio)) < abs((current_x +x_dir-start_x)/ (current_y + y_dir-start_y) - target_xy_ratio)) {
					current_x <- current_x + x_dir
				} else {
					current_y <- current_y + y_dir
					current_x <- current_x + x_dir
				}
				current_xy_ratio <- (current_x - start_x) / (current_y - start_y)
				print(current_xy_ratio)
				cross_sec[i,] <- c(current_x, current_y, sample_matrix[current_x,current_y])
				dummy[current_x, current_y] <- 1
			}
		}
		if (direction == "y") {
			for (i in 1:number_of_points) {
				if (target_xy_ratio == Inf || abs(((current_x + x_dir - start_x) / (current_y + y_dir - start_y) - target_xy_ratio)) < abs(((current_x - start_x) / (current_y + y_dir - start_y) - target_xy_ratio))) {
					current_x <- current_x + x_dir
					current_y <- current_y + y_dir
				} else {
					current_y <- current_y + y_dir
				}
				current_xy_ratio <- (current_x - start_x) / (current_y - start_y)
				print(current_xy_ratio)
				cross_sec[i,] <- c(current_x, current_y, sample_matrix[current_x,current_y])
				dummy[current_x, current_y] <- 1
			}
		}
		
		image(1:Nx, 1:Ny, dummy, col=rainbow(256), add= TRUE)
		query <- as.character(readline("OK?(y/n):"))
	}
	plot(1:number_of_points * distance, cross_sec[,3], main=paste("Cross Section from (", start_x,",",start_y,") to (", end_x, ",", end_y,")", sep=""), xlab="Distance[nm]", ylab="Height[nm]", type = "l")	
}

######################################################
#   5. Execute                                       #
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
#   6. Save  to file                                 #
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
	savelist <- list(search_correlation, start_fit, distance, height, avg_height, peak_to_valley, rms_roughness)
	cat("search_correlation: ", search_correlation, "start_fit: ", start_fit, "drop_last_pixels: ", drop_last_pixels, "distance: ", distance, "height: ", height, "avg_height: ", avg_height, "peak_to_valley: ", peak_to_valley, "rms_roughness: ", rms_roughness, file = paste(filename, "_info.txt", sep=""))
}