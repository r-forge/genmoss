.First.lib <- function(lib, pkg)
{
  library.dynam("libf2c", pkg, lib, getOption("verbose"), .Platform$dynlib.ext, local=FALSE)
  library.dynam("libblas", pkg, lib, getOption("verbose"), .Platform$dynlib.ext, local=FALSE)
  library.dynam("libgslcblas", pkg, lib, getOption("verbose"), .Platform$dynlib.ext, local=FALSE)
  library.dynam("libgsl", pkg, lib, getOption("verbose"), .Platform$dynlib.ext, local=FALSE)
  library.dynam("liblapack", pkg, lib, getOption("verbose"), .Platform$dynlib.ext, local=FALSE)
  library.dynam("libtmglib", pkg, lib, getOption("verbose"), .Platform$dynlib.ext, local=FALSE)
  library.dynam("ttsplits", pkg, lib, getOption("verbose"), .Platform$dynlib.ext, local=TRUE)
  library.dynam("shotgun", pkg, lib, getOption("verbose"), .Platform$dynlib.ext, local=TRUE)
  library.dynam("cvv", pkg, lib, getOption("verbose"), .Platform$dynlib.ext, local=TRUE)
  library.dynam("traintest", pkg, lib, getOption("verbose"), .Platform$dynlib.ext, local=TRUE)

}

.onLoad <- function(lib, pkg)
{
  # This function will be called if NAMESPACE is present.
  # Otherwise .First.lib will be executed.
  # It is run whenever you load up the library in R:
  # library(testPack, lib.loc="/home/briollaislab/olga/try_pack/myinstalled/")

  .First.lib(lib, pkg)
}

.Last.lib <- function(libpath) 
{
  library.dynam.unload("traintest", libpath)
  library.dynam.unload("cvv", libpath)
  library.dynam.unload("shotgun", libpath)
  library.dynam.unload("ttsplits", libpath)
  library.dynam.unload("libtmglib", libpath)
  library.dynam.unload("liblapack", libpath)
  library.dynam.unload("libgsl", libpath)
  library.dynam.unload("libgslcblas", libpath)
  library.dynam.unload("libblas", libpath)
  library.dynam.unload("libf2c", libpath)

}

.onUnload <- function(libpath)
{
  .Last.lib(libpath)
}

