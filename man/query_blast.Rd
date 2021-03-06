% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/blast_search.R
\name{query_blast}
\alias{query_blast}
\title{query_blast}
\usage{
query_blast(db, query, blast, format = 6, evalue = 1e-06, remote = F, ...)
}
\arguments{
\item{db}{path to database(s) to query - dbs of nucl, prot or cdd (fast local searches for conserved domains)}

\item{query}{path to fasta formated query sequence(s)}

\item{blast}{type of blast search - one of blastn, blastp, blastx, tblastn, tblastx, rpsblast, rpstblastn.}

\item{format}{output format of blast query - default to six, leave at this
value if you want this function to work properly}

\item{evalue}{evalue cut-off of the blast search}

\item{remote}{default to FALSE - searches the local database specified in 'db' argument.
If 'remote = TRUE' a remote db (provided by NCBI) is searched. In the latter case the
'db' argument needs to match one of the database specifiers coined by NCBI!
See also 'rentrez::entrez_dbs()'}

\item{...}{to allow for additional parameters}
}
\description{
Local or remote BLAST+ query. Requires standalone blast executive (https://www.ncbi.nlm.nih.gov/books/NBK279690/) and proper setup.
}
\examples{
query_blast()
}
\keyword{blast}
