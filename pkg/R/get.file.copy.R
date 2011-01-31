get.file.copy <- function(dir.in, dir.out, fname="", prefix="", key="", ending="", verbal=TRUE) {
# From given directory dir.in, copies files into dir.out: 
# - either list of file names in fname,
# - or all files from dir.in that start from given prefix and end with ending 
#   and contain keyword key. 
# verbal: only matters if fname="", thus function needs to find files. If files can not be
#   found and verbal=TRUE then a message would be printed. 
# 
# Example:
#
# 

# TODO: delete this later:
#source("get.file.name.R")

# If no list was provided, get files from the directory.
if(all(fname == "")) {
	#print("getting file names")
	fname <- get.file.name(dir=dir.in, prefix=prefix, key=key, ending=ending, verbal=verbal)
}

# If any files are found to copy, copy them one by one if they don't already exist:
if(fname != "" && length(fname) > 0) {

	i <- 1
	while(i <= length(fname)) {

		if(!file.exists(paste(dir.out, fname[i], sep="/")))
			ans <- file.copy(paste(dir.in, fname[i], sep="/"), paste(dir.out, fname[i], sep="/"))
		i <- i + 1
	}
}

}
