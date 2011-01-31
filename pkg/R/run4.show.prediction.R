run4.show.prediction <-
function(filename)
{
	# Returns the list with the following parameters:
	# tpr, fpr, rocarea

	data <- read.table(file=filename, as.is=TRUE);

	nrows <- nrow(data)
	ncols <- ncol(data)
	
	# ***********************************************
	# Gather statistics prior to plotting
	# ***********************************************
	data <- t(data)
	nrows <- ncols # old columns are new rows

	d <- matrix(0,nrows,4)
	d[,1] <- data[,1]
	data <- data[, -1]	

	for(i in (1:nrows)) 
	{	
		d[i,2] <- mean(data[i,])
		d[i,3] <- quantile(data[i,], 0.1, type=5)
		d[i,4] <- quantile(data[i,], 0.9, type=5)	
	}

	lenx <- nrow(d)
	x <- array(0, lenx)	#matrix(0,1,lenx) or array(0, lenx)
	for(i in (1:lenx))
		x[i] <- log(d[i,2]/(1-d[i,2])) + rnorm(1)

	# pdf(file=paste(filename,".plot.pdf",sep="")) # prepare for saving
	
	# **********************************************
	# Plot the predicted values
	# **********************************************

	par(mfrow=c(1,2)) # subplot 1.
	yrange <- c(0,1)  #range(d[,2]) The probability axis
	xrange <- range(x)

	ind <- which(d[,1]==0)
	plot(x[ind],d[ind,2], col="blue", pch=20, ylim=yrange, xlim=xrange, xaxs="i", yaxs="i", xlab="Linear Predictor", ylab="Probability")  # blue "filled small circles"; large enough axis for both blue and red points; axis fit tightly.
	
	for(i in (1:length(ind)))
		matlines(c(x[ind[i]], x[ind[i]]), c(d[ind[i],3], d[ind[i],4]), type="l", col="blue", lty="longdash", lwd=1) # lty: "dashed", "twodash"
	
	par(new=T)	# prevent it from erasing first plot.

	# ------------------------
	ind <- which(d[,1]==1)
	plot(x[ind], d[ind,2], col="red", pch=20, ylim=yrange, xlim=xrange, xaxs="i", yaxs="i", xlab="", ylab="")

	for(i in (1:length(ind)))
		matlines(c(x[ind[i]], x[ind[i]]), c(d[ind[i],3], d[ind[i],4]), type="l", col="red", lty="longdash", lwd=1)
	

	# Horizontal line:
	matlines(xrange, c(0.5,0.5), type="l", color="black", lty="solid", lwd=1)
	par(new=F) # reset back

	# **********************************************
	# Compute for ROC curve
	# **********************************************

	prates <- seq(from=1, to=0, by=-0.01)
	plen <- length(prates)
	tpr <- array(0,c(plen,1))
	fpr <- array(0,c(plen,1))

	for(k in (1:plen))
	{
		for(i in (1:nrow(d)))
		{
			if(d[i,1]==0 && d[i,2]>prates[k])
				fpr[k] <- fpr[k]+1
			if(d[i,1]==1 && d[i,2]>prates[k])
				tpr[k] <- tpr[k]+1
		}
		fpr[k] <- fpr[k]/length(which(d[,1]==0))
		tpr[k] <- tpr[k]/length(which(d[,1]==1))
	}	

	rocarea <- 0;
	for(k in (1:(plen-1)))
		rocarea <- rocarea + (fpr[k+1]-fpr[k]) * (1-tpr[k]+1-tpr[k+1])/2

	# ***********************************************
	# Plot ROC Curve
	# ***********************************************

	plot(0,0, col="white", main="ROC curve", xlab="False positive rate", ylab="True positive rate", xlim=c(-0.05, 1.05), ylim=c(-0.05, 1.05), xaxs="i", yaxs="i") 
	matlines(c(0,1), c(0,1), type="l", lwd=1, lty="solid", col="black")
	matlines(fpr, tpr, type="l", lty="dotdash", col="red", lwd=1.5)

	par(mfrow=c(1,1))	# reset back to one plot per page
	
	#dev.off()	

	return(list(tpr=tpr, fpr=fpr, rocarea=rocarea))
}

