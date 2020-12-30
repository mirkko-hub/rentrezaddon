# run ncbi query ----------------------------------------------------------

#' query_blast
#'
#' Local or remote BLAST+ query. Requires standalone blast executive (https://www.ncbi.nlm.nih.gov/books/NBK279690/) and proper setup.
#'
#' @param db path to database(s) to query - dbs of nucl, prot or cdd (fast local searches for conserved domains)
#' @param query path to fasta formated query sequence(s)
#' @param blast type of blast search - one of blastn, blastp, blastx, tblastn, tblastx, rpsblast, rpstblastn.
#' @param format output format of blast query - default to six, leave at this
#' value if you want this function to work properly
#' @param evalue evalue cut-off of the blast search
#' @param remote default to FALSE - searches the local database specified in 'db' argument.
#' If 'remote = TRUE' a remote db (provided by NCBI) is searched. In the latter case the
#' 'db' argument needs to match one of the database specifiers coined by NCBI!
#' See also 'rentrez::entrez_dbs()'
#'
#' @param ... to allow for additional parameters
#'
#' @keywords blast
#' @export
#' @importFrom magrittr %>% %$%
#' @examples
#' query_blast()

query_blast <- function(db, query, blast, format = 6, evalue = 1e-6, remote = F, ... ) {
  blast = rlang::enexpr(blast)
  blastn = "/opt/ncbi-blast-2.10.1+/bin/blastn"
  blastp = "/opt/ncbi-blast-2.10.1+/bin/blastp"
  blastx = "/opt/ncbi-blast-2.10.1+/bin/blastx"
  tblastn = "/opt/ncbi-blast-2.10.1+/bin/tblastn"
  tblastx = "/opt/ncbi-blast-2.10.1+/bin/tblastx"
  rpsblast = "/opt/ncbi-blast-2.10.1+/bin/rpsblast"
  rpstblastn = "/opt/ncbi-blast-2.10.1+/bin/rpstblastn"

  colnames <- c("qseqid",
                "sseqid",
                "pident",
                "length",
                "mismatch",
                "gapopen",
                "qstart",
                "qend",
                "sstart",
                "send",
                "evalue",
                "bitscore")

  args_local <- c("-db", db,
                  "-query", query,
                  "-outfmt", format,
                  "-evalue", evalue)
  args_remote <- c("-remote",
                   "-db", db,
                   "-query", query,
                   "-outfmt", format,
                   "-evalue", evalue)


  if (remote == F) {
    out <- system2(command = eval(blast),
                   args = args_local,
                   wait = TRUE,
                   stdout = TRUE)
  }
  if (remote == T) {
    out <- system2(command = eval(blast),
                   args = args_remote,
                   wait = TRUE,
                   stdout = TRUE)
  }
  out  %>%
    tibble::enframe(name = NULL) %>% # change to enframe(name = NULL)
    tidyr::separate(col = value,
                    into = colnames,
                    sep = "\t",
                    convert = TRUE) %>%
    dplyr::mutate(database = as.character(db),
                  qseqid = as.character(qseqid),
                  sseqid = as.character(sseqid),
                  pident = as.double(pident),
                  length = as.integer(length),
                  mismatch = as.integer(mismatch),
                  gapopen = as.integer(gapopen),
                  qstart = as.integer(qstart),
                  qend = as.integer(qend),
                  sstart = as.integer(sstart),
                  send = as.integer(send),
                  evalue = as.double(evalue),
                  bitscore = as.integer(bitscore)
                  )
}


# query_blast <- function(db, query, blast, format = 6, evalue = 1e-6, ... ) {
#   blast = rlang::enexpr(blast)
#   blastn = "/home/mirkko/bin/ncbi-blast/ncbi-blast-2.8.1+/bin/blastn"
#   blastp = "/home/mirkko/bin/ncbi-blast/ncbi-blast-2.8.1+/bin/blastp"
#   blastx = "/home/mirkko/bin/ncbi-blast/ncbi-blast-2.8.1+/bin/blastx"
#   tblastn = "/home/mirkko/bin/ncbi-blast/ncbi-blast-2.8.1+/bin/tblastn"
#   tblastx = "/home/mirkko/bin/ncbi-blast/ncbi-blast-2.8.1+/bin/tblastx"
#
#   colnames <- c("qseqid",
#                 "sseqid",
#                 "pident",
#                 "length",
#                 "mismatch",
#                 "gapopen",
#                 "qstart",
#                 "qend",
#                 "sstart",
#                 "send",
#                 "evalue",
#                 "bitscore")
#
#   system2(command = eval(blast),
#           args = c("-db", db,
#                    "-query", query,
#                    "-outfmt", format,
#                    "-evalue", evalue),
#           wait = TRUE,
#           stdout = TRUE) %>%
#     tibble::enframe(name = NULL) %>% # change to enframe(name = NULL)
#     tidyr::separate(col = value,
#              into = colnames,
#              sep = "\t",
#              convert = TRUE) %>%
#     dplyr::mutate(database = db)
#
# }

# should work just need proper files!

# query_blast(db = "/home/mirkko/Documents/NCBI/ncbi_nuccore/db_nuccore/1011336839",
#             blast = tblastn, query = "/home/mirkko/Documents/NCBI/queries/pcc7120_rbcl.faa")



# blastcmd ----------------------------------------------------------------

#' query_blastdbcmd
#'
#' This function extracts sequences of a given range from a local database - important - with PARSED IDs!
#'
#' @param db path to database
#' @param entry any identifier which is recognized by the local database
#' @param dbtype type of database - one of 'guess', 'nucl', 'prot'.
#' @param range optional; sequence range, given as '1-10'
#' @param strand optional; one of 'plus', 'minus'.
#' @param filepath optional; with append = TRUE
#' @param ... to allow for additional parameters
#'
#' @keywords blast
#' @export
#' @importFrom magrittr %>%
#' @importFrom Biostrings AAStringSet writeXStringSet
#' @importFrom stringr str_remove str_flatten
#' @examples
#' query_blastdbcmd()
#'
#' query_blastdbcmd(db = "/home/mirkko/Documents/NCBI/ncbi_databases/db_genomes_parsed_id/GCF_000012525.1_ASM1252v1_protein",
#' entry = "WP_011242444.1",
#' dbtype = "prot",
#' filepath = "/home/mirkko/Documents/NCBI/testtesttest.faa")
#'
#' c("WP_071818122.1", "WP_011242444.1") %>% map(., query_blastdbcmd,
#' db = "/home/mirkko/Documents/NCBI/ncbi_databases/db_genomes_parsed_id/GCF_000012525.1_ASM1252v1_protein",
#' dbtype = "prot", filepath = "/home/mirkko/Documents/NCBI/testtesttest.faa")

query_blastdbcmd <- function(db, entry, dbtype, range, strand, filepath, ... ) {
  blastdbcmd = "/opt/ncbi-blast-2.10.1+/bin/blastdbcmd"

  args_local <- c("-db", db,
                  "-entry", entry,
                  "-dbtype", dbtype,
                  "-outfmt", c('%f'))

  if(missing(range) == FALSE) {
    arg_range <- c("-range", range)
    args_local <- c(args_local, arg_range)}

  if(missing(strand) == FALSE) {
    arg_strand <- c("-strand", strand)
    args_local <- c(args_local, arg_strand)}

  out <- system2(command = blastdbcmd,
                 args = args_local,
                 wait = TRUE,
                 stdout = TRUE)

  aastringset <- out[2:length(out)] %>% str_flatten() %>% AAStringSet(.)
  names(aastringset) <- out[1] %>% str_remove(. ,">")

  if(missing(filepath) == FALSE){
    aastringset %>% writeXStringSet(., filepath = filepath, append = TRUE, format = "fasta")
  } else {aastringset}
}


# query_blastdbcmd_cdd ----------------------------------------------------

#' query_blastdbcmd_cdd
#'
#' This function extracts cdd descriptions from cdd motifs in a local cdd database.
#'
#' @param db path to cdd_database (eg. "/home/mirkko/Documents/NCBI/ncbi_databases/cdd/Cdd_NCBI")
#' @param entry any identifier which is recognized by the local database, ideally obtaind by querying this local cdd database (query_blast, rpsblast)
#' @param ... to allow for additional parameters
#'
#' @keywords blast
#' @export
#' @importFrom tibble tibble
#' @examples
#' query_blastdbcmd_cdd(db = "/home/mirkko/Documents/NCBI/ncbi_databases/cdd/Cdd_NCBI", entry = "CDD:197200")


query_blastdbcmd_cdd <- function(db, entry, ... ) {
  blastdbcmd = "/opt/ncbi-blast-2.10.1+/bin/blastdbcmd"

  args_local <- c("-db", db,
                  "-entry", entry)

  out <- system2(command = blastdbcmd,
                 args = args_local,
                 wait = TRUE,
                 stdout = TRUE)

  tibble(cdd = entry,
         cdd_title = out[1] %>% str_extract(., paste0("(?<=", entry, " ).*?(?=\\.)")),
         description = out[1])
}
