% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/public_functions.R
\name{alignSequences}
\alias{alignSequences}
\title{Align sequences in a seqlist}
\usage{
alignSequences(
  monolist,
  subset,
  alignment_directory_path,
  sequences_of_interest_directory_path,
  input_sequence_type = c("nucl", "amin"),
  mode = c("nucl_align", "amin_align", "codon_align", "fragment_align"),
  base_fragment = NULL
)
}
\arguments{
\item{monolist}{Monolist of the sequences to be aligned. First column should be "accession"}

\item{subset}{TRUE/FALSE column in monolist that specifies which sequences should be included in the alignment}

\item{alignment_directory_path}{Path to where the alignment should be written}

\item{sequences_of_interest_directory_path}{Path to a directory where blast hits should be written out as fasta files}

\item{input_sequence_type}{One of "nucl" or "amin"}

\item{mode}{One of "nucl_align", "amin_align", or "codon_align"}

\item{base_fragment}{TROUBLESHOOTING: 
"ERROR: inconsistency between the following pep and nuc seqs" - usually means there are duplicate accessions numbers in the input monolist}
}
\description{
Align sequences in a seqlist
}
