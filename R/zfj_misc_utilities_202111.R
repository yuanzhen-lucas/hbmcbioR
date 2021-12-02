#' Check if a sequence is from Truseq adaptor
#'
#' @param test_seq the sequence or motif consensus to compare with Truseq ligand
#' @param lib_seq the sequence of full length Truseq ligand
#'
#' @return the result of alignment of test_seq, lib_seq, and the reverse complement of test_seq
#' @export
#'
#' @examples
is_from_adaptor<-function(test_seq="ACGTGTGCTCTTC",lib_seq="AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnAGATCGGAAGAGCACACGTCTGAACTCCAGTCACnnnnnnnnATCTCGTATGCCGTCTTCTGCTTG")
{
  seqs <- c(test_seq,fjComm::revComp(test_seq),lib_seq) %>% DNAStringSet()
  names(seqs)=c("test","test_rc","lib")
  aln_result=DECIPHER::AlignSeqs(myXStringSet = seqs)
  DECIPHER::BrowseSeqs(aln_result)
}


