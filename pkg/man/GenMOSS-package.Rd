\name{GenMOSS-package}
\Rdversion{1.1}
\alias{GenMOSS-package}
\alias{GenMOSS}
\docType{package}
\title{
Application of MOSS algorithm to dense SNP array data
}
\description{
Performs genome-wide analysis of dense SNP array data using the mode oriented stochastic search (MOSS) algorithm in a case-control design. Finds combination of best predictive SNPs associated with the response. The identified regression models are then tested by peforming cross-validation and prediction in a test set. Includes function for visualization of the obtained results. Includes preprocessing of the data from Plink format to the format required by the MOSS algorithm.
}
\details{
\tabular{ll}{
Package: \tab GenMOSS\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2009-08-12\cr
License: \tab GPL (>=2)\cr
}

System Requirements:
\preformatted{
 * Linux
 * 64 bit machine
 * MaCH software (http://www.sph.umich.edu/csg/abecasis/MACH/download/)
}

This package contains the source code of the following open-source softwares:
\preformatted{
 * GSL C++ library (http://www.gnu.org/software/gsl/)
 * CLAPACK library (http://www.netlib.org/clapack/)
}



The package consists of three groups of files: preprocessing functions, main \acronym{MOSS} functions, and helper functions. The name of the the first two groups of functions begins with "pre" and "run", respectively. The preprocessing ("pre") functions are necessary for converting data from Plink format to required binary \acronym{MOSS} format. The main \acronym{MOSS} functions ("run") are needed to perform the model selection, cross-validation, data prediction, as well as plotting the results. The helper functions are available for user's convenience to check things out for their datasets. We describe basic steps for "pre" and "run" functions below.

\preformatted{
Preprocessing Functions
-----------------------
}


The preprocessing step converts data from Plink format (\code{\link{ex2plink}} describes the Plink format) to the format required by the \acronym{MOSS} algorithm. Frequently geno data has missing values, for their imputation we use \acronym{MaCH} software (\url{http://www.sph.umich.edu/csg/abecasis/MACH/download/}). This imputation may require to run \acronym{MaCH} algorithm on one chromosome at a time, thus all preprocessing steps deal with multiple files: one for each chromosome. There is a total of 10 preprocessing steps that should be run in their proper order (the names of these functions begin with "pre" followed by the sequence number, followed by short description of what it does). Thus the number of intermediate files generated will be very large, for which good organization of files into directories is necessary. It is recommended to use the directory structure of the format created by \code{\link{pre0.dir.create}}. 


Almost every preprocessing function has two versions: normal mode and batch mode. In normal mode, users are requested to provide input and output directory names, full names of the required files, and some other additional parameters specific to the task. Whereas the batch mode is designed to run the function for ALL the files in the input directory that satisfy a naming criterion. This batch mode saves the user from having to call the same function 22-25 times for each chromosome. The naming criterion is as follows:

\preformatted{
 * prefix - The beginning string of the file name up until the chromosome number. 
     Here the assumption is that when a dataset is split into 22-25 files, 
     one chromosome in each, then the beginning of the file name is usually the
     same, followed by the chromosome number.
     E.x. Files with names:
        ~ "geno.data_chr1.my.ped"
        ~ "geno.data_chr2.my.ped"
        ~ "geno.data_chr3.my.ped"
          ...
        ~ "geno.data_chr22.my.ped"
     They all share the same beginning string:
        "geno.data_chr" - this is the 'prefix' for the above example.
        Note that it must be immediately followed by chromosome number.
	Also the chromosome number is expected to be a 1- or 2-digit number.
        Rename all X, Y, M, etc, to some 2-digit number.

 * key - Any string that appears in the file name. In case that the input directory
     contains files that begin with the same prefix, but should not be processed
     by the function, this parameter gives additional flexibility to filter 
     such files out.
     E.x. Suppose input directory contains the following files:
        ~ "geno.data_1.CASE.ped"
        ~ "geno.data_2.CASE.ped"
          ...
        ~ "geno.data_22.CASE.ped"
        ~ "geno.data_1.CONTROL.ped"
        ~ "geno.data_2.CONTROL.ped"
          ...
        ~ "geno.data_22.CONTROL.ped"
        ~ "geno.data_1.short_try.ped"
     First note that they all have the same prefix = "geno.data_".
     Now if you wish to specify that only CASE files should be processed,
     set key="CASE" - this will ignore all CONTROL files. Also it will ignore all 
     those testing files like "geno.data_1.short_try.ped", which might have
     been manually created by users for testing purposes.
     Note: this key is usually optional: if the input directory contains ONLY the
     files that need to be processed, then key can be set to an empty string "".

 * ending - A string that appears at the end of the file name. Normally this does
     not have to be the filename extension, unless specifically stated. The ending
     should not include chromosome number. If preprocessing functions are run in
     their proper order, then the suggested default values for endings in the 
     preprocessing functions should apply.
     Ex. 
        ~ "geno.data_1.CASE.ped" - ".ped" or "d" or "CASE.ped" or "E.ped", etc.
        ~ "geno.data_2.CASE" - "CASE" or ".CASE" or "" or "E" or "SE", etc.
        ~ "geno.data_15.CONTROL.dat" - ".dat" or "t" or "CONTROL.dat", etc.

 * Note: it is preferable to name files such that they have a filename extension,
    Ex.
        ~ good: "geno.data_1.CASE.ped";  bad: "geno.data_1.CASE"
        ~ good: "CGEM.chr11CONTROL.dat"; bad: "CGEM.chr11CONTROL"
   Sometimes preprocessing functions name their output functions by slightly 
   modifying the name of the input file. When this is done, filename extension
   is usually removed. For example, suppose function wants to add word 
   "_cleaned.txt" to the end of your filename "CGEM.chr_12CONTROL.ped"
   Resultant filename would be: "CGEM.chr_12CONTROL_cleaned.txt", 
   since ".ped" will be identified as filename extension and will be lost.

   Consider what happens if you are not using filename extensions:
   then filename "CGEM.chr_12CONTROL" will be renamed as "CGEM_cleaned.txt",
   since the entire ".chr_12CONTROL" will be identified as file name extension,
   but it contains valuable chromosome information that will be lost.

   Thus always use file name extensions: ".ped", ".dat", ".txt", ".map", etc.
}  

It is recommended to run the preprocessing functions in the following order:

\preformatted{

 * pre0.dir.create - creates a set of empty directories d0 to d10. 
 * get.file.copy - copy original format files to dir d0.
 * ex2plink - modify this function, or write something similar to 
                      convert your format into Plink, this may involve
                      splitting dataset into multiple files: one per 
                      chromosome; place the result into dir d1.
 * pre1.plink2mach.batch - converts Plink format to MaCH's input format,
                      which splits each chromosome into CASE and CONTROL
                      files; store result into dir d2.
 * pre2.remove.genos.batch - remove all SNPs that have too many missing 
                      values, store result into dir d3.
 * pre3.call.mach.batch - imputes missing values using MaCH1, store
                      results in dir d5 (current version does not use d4).
 * pre4.combine.case.control.batch - combines CASE and CONTROL files,
                      place result into dir d6.
 * pre5.genos2numeric.batch - convert data from "A/G", "C/T", "G/G", etc
                      format to 3 levels: 1, 2, 3; store into dir d7.
 * pre6.discretize.batch - convert 3 levels: 1, 2, 3 into binary: 0, 1;
                      store result into dir d8.
 * pre7.merge.genos - merges all files across all chromosomes into one,
                      result should go into dir d9.
 * pre9.split.train.test.batch - split the full dataset into train and 
                      test files; save the result into d10. 
}

\preformatted{
MOSS Functions
---------------
}


After the preprocessing steps are complete, continue running the main steps in the order described below. Note that main functions always require the full file name, including the directory, and their output will always go into that same directory. Thus if you start from directory d10, then all the output will end up going to dir d10 as well. 

\preformatted{

 * run1.moss.regression - perform MOSS search for log-linear models.
 * run2.prediction.cvv - does prediction by cross-validation using the 
                      regression models identified in "run1" step.
 * run3.prediction.train.test - performs prediction on the test file, using 
                      the regression models identified in "run1" step.
 * run4.save.prediction - saves a plot (as .pdf file) of the predicted values 
                      and the corresponding ROC curve for the resulting 
                      predictions from "run2" and/or "run3". To see the plots,
                      open the .pdf files (they should be in the same input 
                      directory, d10).
 * run4.show.prediction - shows the plot, without saving it - useful only if
                      you have graphical interface to see R's plots.
 * run5.brier - computes the Brier score, its mean and standard deviation.

}

To see the functionality of preprocessing and MOSS algorithm, try running:


demo("gendemo")

}
\author{
Author: Olga Vesselova, Laurent Briollais, Adrian Dobra, Helene Massam.

Maintainer: <laurent@lunenfeld.ca> 
}
\references{
\cite{Dobra, A., Briollais, L., Jarjanazi, H., Ozcelic, H. and Massam, H. (2008). Applications of the mode oriented stochastic search (MOSS) agorithm for discrete multi-way data to genomewide studies. Bayesian Modelling in Bioinformatics (D. Dey, S. Ghosh and B. Mallick, eds.), Taylor & Francis. To appear.}
}
\keyword{package}
\examples{
write(rbinom(200,1,0.5), file="randbinary.txt", append=FALSE, sep=" ", ncolumns=50)
run1.moss.regression("randbinary.txt")
try(system("rm randbinary.txt*"))
}
