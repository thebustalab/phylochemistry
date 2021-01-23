% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/public_functions.R
\name{blastTranscriptomes}
\alias{blastTranscriptomes}
\title{BLAST search local transcriptomes}
\usage{
blastTranscriptomes(transcriptomes, query_in_path,
  sequences_of_interest_directory_path, blast_module_directory_path,
  blast_type = c("blastn", "dc-megablast"), remove = NULL,
  monolist_out_path)
}
\arguments{
\item{transcriptomes}{A list of paths to the transcriptomes that should be searched, named by taxonomic identifier (e.g. species names)}

\item{query_in_path}{Path to a fasta file containing the query}

\item{sequences_of_interest_directory_path}{Path to a directory where blast hits should be written out as fasta files}

\item{blast_type}{One of "blastn" or "dc-megablast"}

\item{remove}{Names of columns to remove from the output monolist}

\item{monolist_out_path}{Path to where the output monolist should be written}
}
\description{
Search locally stored transcriptomes for query sequences.
}