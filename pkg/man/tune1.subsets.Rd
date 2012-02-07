\name{tune1.subsets}
\Rdversion{1.1}
\alias{tune1.subsets}
\title{
	Imputes dense map of SNPs on chromosome regions with MaCH
}
\description{
	For chromosomes and their small regions specified, run MaCH1 with hapmap to get more detailed sampling of SNPs in the region, and prepares this subset of data to be processed by MOSS algorithm.
}
\usage{
tune1.subsets(dir.dat, dir.ped, dir.ann, dir.pos.snp, dir.pos.ann, dir.pos.hap, dir.out, prefix.dat, prefix.ped, prefix.ann, prefix.pos.snp, prefix.pos.ann, prefix.pos.hap, key.dat = "", key.ann = "", key.pos.ann = "", key.pos.hap = "", ending.dat = ".dat", ending.ped = ".ped", ending.ann = ".map", ending.pos.snp = ".snps", ending.pos.ann = "annotation.txt", ending.pos.hap = ".hap.gz", pos.list.triple, ped.nonsnp = 5, ann.header=FALSE, ann.snpcol=2, ann.poscol=4, ann.chrcol=0, pos.ann.header = TRUE, pos.ann.snpcol = 5, pos.ann.poscol = 2, pos.hap.nonsnp = 2, out.name.subdir = "seg1", out.prefix = "subdata", rsq.thresh = 0.5, num.iters = 2, hapmapformat = FALSE, mach.loc = "/software/mach1")
}
\arguments{
  \item{dir.dat}{
	The name of directory where file listing SNPs of the dataset can be found.
}
  \item{dir.ped}{
	The name of directory where file with data of the dataset can be found.
}
  \item{dir.ann}{
	The name of directory where SNP position information for the dataset can be found. Note: this file must contain position information about all SNPs that are listed in .dat; all other SNPs will be ignored.
}
  \item{dir.pos.snp}{
	The name of directory where hapmap SNP list can be found.
}
  \item{dir.pos.ann}{
	The name of directory where hapmap annotation file containing position information can be found.
}
  \item{dir.pos.hap}{
	The name of directory where the hapmap zipped data can be found.
}
  \item{dir.out}{
	The name of directory to which output folder should be placed.
}
  \item{prefix.dat}{
	The beginning of the file name for dataset's list of SNPs.
}
  \item{prefix.ped}{
	The beginning of the file name for dataset's data.
}
  \item{prefix.ann}{
	The beginning of the file name for dataset's SNP position information.
}
  \item{prefix.pos.snp}{
	The beginning of the file name for hapmap's list of SNPs.
}
  \item{prefix.pos.ann}{
	The beginning of the file name for hapmap's SNP position information.
}
  \item{prefix.pos.hap}{
	The beginning of the file name for hapmap's data.
}
  \item{key.dat}{
	Any keyword in the name of dataset's list of SNPs.
}
  \item{key.ann}{
	Any keyword in the name of dataset's SNP position information.
}
  \item{key.pos.ann}{
	Any keyword in the name of hapmap's SNP position information.
}
  \item{key.pos.hap}{
	Any keyword in the name of hapmap's data.
}
  \item{ending.dat}{
	The ending of dataset's list of SNPs filename.
}
  \item{ending.ped}{
	The ending of dataset's data filename.
}
  \item{ending.ann}{
	The ending of dataset's SNP position information filename.
}
  \item{ending.pos.snp}{
	The ending of hapmap's list of SNPs filename.
}
  \item{ending.pos.ann}{
	The ending of hapmap's SNP position information filename.
}
  \item{ending.pos.hap}{
	The ending of hapmap's data filename.
}
  \item{pos.list.triple}{
	A list of chromosomes and position boundaries to be expanded upon. The list should contain information in the order: (chromosome number, start position, end position, chromosome number, start position, end position, etc.). This allows users to specify multiple chromosomes with multiple regions within each chromosome. For example, specifying region of positions 6000-19000 and 111000-222000 in chrom 15, together with positions 55000-77000 in chrom 21, can be listed as: c(15, 6000, 19000, 15, 111000, 222000, 21, 55000, 77000). Note that MaCH will be run on one chromosome at a time, and for all its specified regions. 
}
  \item{ped.nonsnp}{
	The number of non-snp leading columns in dataset's data file. For example input to MaCH format has 5 columns, Plink has 6 columns.
}
  \item{ann.header}{
        Whether or not hapmap's SNP position information file has a header. Ex. .annotation.txt = TRUE, .legend.txt = TRUE. Since format of the hapmap file is not hard-coded, specify the format of your prefered hapmap library; the defaults are set to the 1000 Genome data (from MaCH website).
}
  \item{pos.ann.snpcol}{
        The column number in hapmap's SNP position information file that contains SNP names/ids. For example in .annotation.txt it's column 5; in .legend.txt it's column 1.
}
  \item{pos.ann.poscol}{
        The column number in hapmap's SNP position information file that contains position information. For example in .annotation.txt it's column 2; in .legend.txt it's also column 2.
}
  \item{pos.ann.header}{
	Whether or not dataset's SNP position information file has a header. Ex. .map = FALSE, but other formats might have a header. Since format of this file is not hard-coded, specify the format that your dataset comes with.
}
  \item{ann.snpcol}{
	The column number in dataset's SNP position information file that contains SNP names/ids. For example in .map it is column 2.
}
  \item{ann.poscol}{
	The column number in dataset's SNP position information file that contains position information. For example in .map it is column 4.
}
  \item{ann.chrcol}{
        The column number in dataset's SNP position information file that contains chromosome number information. In .map format there is no such column, since there is a unique file per chromosome, thus default for this parameter is 0. In case if all position information is included in one single file for all/many chromosomes, specify which column corresponds to chromosome number.
}
  \item{pos.hap.nonsnp}{
	The number of non-SNP leading columns in hapmap's data file. In .hap.gz it is 2.
}
  \item{out.name.subdir}{
	The name of subdirectory structure to be created for output for this sequence of chromosomes and positions. Note: this folder name MUST be different for each different set of chromosome and position boundaries triplets. 
}
  \item{out.prefix}{
	The beginning of output file names.
}
  \item{rsq.thresh}{
	Threshold for RSQ of MaCH's imputation. Recommended default is 0.5.
}
  \item{num.iters}{
        The number of iterations MaCH should make in its first step to estimate its model parameters.
}
  \item{hapmapformat}{
        The type of haplotype data format: 1000G haplotype dataset has .snps file with one column, so \code{\var{hapmapformat}} defaults to FALSE. Another dataset format listing SNPs (.legend.txt) has 4 columns - change \code{\var{hapmapformat}} to TRUE.
}

  \item{mach.loc}{
	The location directory where "mach" executable can be found.
}
}
\details{
	The input files for this function are inteded to be from different folders of the subdirectory structure used in preprocessing steps (see \code{\link{pre0.dir.create}}). The dataset's SNP (.dat) and data (.ped) information are intended to come from d3 (d03_removed); whereas the dataset's position information (.map) can be obtained from d1 (d01_plink) subdirectory. The hapmap files are huge and can be used by many datasets, thus there is no need to keep a copy of them in our subdirectory structure for each dataset. Note: if the hapmap file that specifies SNP information ALSO lists their position information, simply provide that file (and it's column format) to this function twice (as \code{\var{prefix.pos.snp}} and \code{\var{prefix.pos.ann}}).
	This function is meant to begin from early pre-processing steps, re-run MaCH with hapmap on desired regions, then combine CASE with CONTROL, and call all the pre-processing functions in sequence up until \code{\link{pre7.merge.genos}}. At the end, the output will be a single file ready to be called by MOSS \code{\link{run1.moss.regression}}. 
	A new convenient subdirectory structure will be created, similar to \code{\link{pre0.dir.create}} within new directory \code{\var{out.name.subdir}}. 
	This function requires two sets of data: user's dataset and reference haplotypes. There are many hapmap libraries for download from the web, so this function tries to be as general as possible to allow users to give column information about the format. MaCH also needs to understand the given hapmap format. The defaults are set for 1000G Phase I(a) from MaCH's website: \url{http://www.sph.umich.edu/csg/abecasis/MaCH/download/1000G-PhaseI-Interim.html}. Note: the data file (.hap.gz) is expected to be zipped. However please unzip the .annotation.txt file before calling this function. 
	The first thing this function would do is extract the given position intervals from user's datafiles and from haplotype files. This would make both files smaller so that running MaCH is feasible. 
	MaCH will be run on CASE and CONTROL data files separately. After MaCH is run with hapmap, most of the predicted SNPs would have very low RSQ score, thus out of thousands of SNPs that are within the region in hapmap file, only hundreds will be actually reliable. This function prunes out all the SNPs with RSQ score lower than \code{\var{rsq.thresh}}. Then CASE and CONTROL will be combined based on common remaining SNPs. 
	Then the function will run the three preprocessing functions (\code{\link{pre5.genos2numeric.batch}},  \code{\link{pre6.discretize.batch}}, \code{\link{pre7.merge.genos}}) to output the final ready-to-use file. 
	 
}
\value{
	The FULL name of the combined result geno file (including the directory).
}
\references{
	MaCH website: \url{http://www.sph.umich.edu/csg/abecasis/MACH/download/}
}
\author{
Olia Vesselova
}

\seealso{
\code{\link{pre2.remove.genos}}, \code{\link{pre2.remove.genos.batch}},
\code{\link{pre3.call.mach}}, \code{\link{pre4.combine.case.control}},
\code{\link{pre4.combine.case.control.batch}},
\code{\link{pre5.genos2numeric}}, \code{\link{pre5.genos2numeric.batch}},
\code{\link{pre6.discretize}}, \code{\link{pre6.discretize.batch}}, 
\code{\link{pre7.merge.genos}}, \code{\link{run1.moss.regression}}
}
\examples{
print("See the demo 'gendemo'.")
}
