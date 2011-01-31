run5.brier <-
function(filename)
{
	# Computes Brier score (takes time to compute).
	# Returns the list with the following parameters:
	# briermean, brierstd

	data <- read.table(file=filename, as.is=TRUE);
	nrows <- nrow(data)

	brier <- apply(data[2:nrows,], 1, 'diff.square.sum', data[1,])
	brier <- c(0,brier)		# Since 1st element must be 0

	briermean <- mean(brier)
	brierstd <- sd(brier)
	
	return(list(briermean=briermean, brierstd=brierstd))
}

