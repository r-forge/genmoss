\name{pre8.add.conf.var.unix}
\Rdversion{1.1}
\alias{pre8.add.conf.var.unix}
\title{
        Append confounding variables using Linux
}
\description{
        Uses Linux functions to append confounding variables listed in \code{\var{file.conf}} to the end of the \code{\var{file.name}}, right before the disease status (last) column. The output will contain only the patients for which confounding variables exist (other patients will be omitted), so new family file will be written. This function is similar to \code{\link{pre8.add.conf.var}}, only it avoids having to load up into memory the \code{\var{file.name}} (since this file can be very large).
}
\usage{
pre8.add.conf.var.unix(file.name, dir.file, file.fam, dir.fam = dir.file, 
file.conf, dir.conf = dir.file, file.out, fam.out = file.fam, dir.out)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file.name}{
        The name of the binary data file. The format of this file should have last column as the disease status, tab separated, no header.
}
  \item{dir.file}{
        The name of directory where \code{\var{file.name}} can be found.
}
  \item{file.fam}{
        The name of the family file. Format: one column - one patient ID per line.
}
  \item{dir.fam}{
        The name of directory where \code{\var{file.fam}} can be found.
}
  \item{file.conf}{
        The name of the file that contains confounding variable information. The file should be in the following format:
\preformatted{
    patientID1 1       2 ...
    patientID2 3       1 ...
    patientID3 2       2 ...
    ...

    - Column 1: patient ID, exactly the same names should appear in 
	file.fam;
               * order does not matter;
               * some patients may be missing;
               * no new patients should appear in file.conf (if they 
		don't exist in file.fam)
    - Column 2: the confounding variable must have no more than 3 
	different values.
    - Other columns are optional, may be included if there are 
	more confounding variables (3 categories each)
    - No header
    - Tab separated
    - No missings or NAs
}
}
  \item{dir.conf}{
        The name of directory where \code{\var{file.conf}} can be found.
}
  \item{file.out}{
        The name of the output file, which will contain all information of \code{\var{file.name}}, plus confounding variables, only for the patients mentioned in \code{\var{file.conf}}.
}
  \item{fam.out}{
        The name of the family output file.
}
  \item{dir.out}{
        The name of directory to which \code{\var{file.out}} and \code{\var{fam.out}} should be saved.
}
}
\author{
Olia Vesselova
}


\seealso{
\code{\link{pre8.add.conf.var}}, \code{\link{pre7.merge.genos}}, 
\code{\link{pre9.split.train.test}}, \code{\link{pre9.split.train.test.batch}}
}
\examples{
print("See the demo 'gendemo'.")
}
