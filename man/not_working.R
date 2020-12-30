


# blast_result <- pcc7120_groel_nuccore[1,]
# xstring <- readDNAStringSet(filepath = "/home/mirkko/Documents/NCBI/ncbi_nuccore/fasta/17227497.fna", format = "fasta")
# path <- "/tmp/seqs.fna"
# nuc_up = 20000
# nuc_down = 300


ex_1 <- tibble(qseqid = c("a", "b", "c"),
               sstart = c(4, 8, 25),
               send = c(8, 4, 25))

ex_2 <- tibble(qseqid = c("a"),
               sstart = c(12),
               send = c(11))



bliblablupp <- function(blast_result) {

  if(dim(blast_result)[1] > 1) {

    blast_result <- blast_result %>%
      mutate(., row = c(1:dim(blast_result)[1]))

    1:3 %>% map_dfr(blast_result, define_strand)

  } else {
    hit_position <- blast_result %>% define_strand(.)
  }
  hit_position
}


bliblablupp(ex_1)


define_strand <- function(blast_result) {



  if(blast_result$sstart == blast_result$send) return(
    "Invalid blast hit: sstart equals send"
  )

  if(blast_result$sstart > blast_result$send) {
    hit_position <- blast_result %>%
      mutate(., sstart2 = send, send2 = sstart,
             sstart = sstart2, send = send2,
             strand = '-') %>%
      select(-c(sstart2, send2))

  } else {
    hit_position <- blast_result %>%
      mutate(strand = '+')
  }
  hit_position
}


define_strand(ex_2)


set_subseq_limits <- function(hit_position) {

  # do we need to define add and sub differently for +/- strand situation ???
  extract <- hit_position %>% mutate(sstart_extr = case_when(sstart >= sub_5 ~ sstart - sub_5,
                                                             sstart < sub_5 ~ 1),
                                     send_extr = case_when(width(contig) >= send + add_3 ~ send + add_3,
                                                           width(contig) < send ~ width(contig) %>% as.double()) )

}


extract_subseq <- function(xstring, start, end, rev_comp = F) {

  # do input checking

  subsequence <- xstring %>% subseq(., start = start, end = end)
  if (rev_comp == T) {
    subsequence <- subsequence %>% reverseComplement(.)
    subsequence
  } else {
    subsequence
  }
}






fooo <- function(blast_result, xstring, nuc_buffer, path) {

  if (identical(blast_result$sseqid, names(xstring)) == TRUE) {
  } else {
    contig_names <- names(xstring) %>% map_dfr(., extract_id_from_xstring)
    names(xstring) <- contig_names
  }

  # todo: useful error mesagges are required: like in the case that sseqid has not been found in names(xstring)

  hit_position <- blast_result %>% select(sseqid, sstart, send)
  contig <- xstring[xstring@ranges@NAMES == hit_position$sseqid]
  nuc_down <- nuc_buffer[1]
  nuc_up <- nuc_buffer[2]

  # todo: include option for space +/- in function params
  # test positioning of nuc_buffer, evtl nuc_down as negative value ?
  # test what to do if buffer not given?!

  if (hit_position$sstart > hit_position$send) {
    hit_position <- hit_position %>% mutate(., sstart2 = send, send2 = sstart, sstart = sstart2, send = send2)
    extract <- hit_position %>% mutate(sstart_extr = case_when(sstart >= nuc_down ~ sstart - nuc_down,
                                                               sstart < nuc_down ~ 1),
                                       send_extr = case_when(width(contig) >= send + nuc_up ~ send + nuc_up,
                                                             width(contig) < send ~ width(contig) %>% as.double()) )
    fasta <- contig %>% subseq(., start = extract$sstart_extr, end = extract$send_extr) %>% reverseComplement(.)
  } else {
    extract <- hit_position %>% mutate(sstart_extr = case_when(sstart >= nuc_down ~ sstart - nuc_down,
                                                               sstart < nuc_down ~ 1),
                                       send_extr = case_when(width(contig) >= send + nuc_up ~ send + nuc_up,
                                                             width(contig) < send ~ width(contig) %>% as.double()) )
    fasta <- contig %>% subseq(., start = extract$sstart_extr, end = extract$send_extr)
  }
  writeXStringSet(fasta, filepath = path, append = T)
}

fooo(blast_result = pcc7120_groel_nuccore[2,],
     xstring = readDNAStringSet(filepath = "/home/mirkko/Documents/NCBI/ncbi_nuccore/fasta/17227497.fna", format = "fasta"),
     path = "/tmp/seqs.fna",
     nuc_buffer = c(345, 200))

fasta <- readAAStringSet(filepath = "/tmp/seqs.fna", format = "fasta")
