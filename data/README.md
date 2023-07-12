# Input data
This directory contains input data for the analysis.

## Basic information about sequences and alignments

These files are used for the basic processing of the deep sequencing data to call variants by barcode and count barcodes:

   - `PacBio_amplicons_*_.gb`: the amplicons being sequenced by PacBio.

   - `feature_parse_specs_*_.yaml`: how to parse the amplicon when handling the PacBio data.

   - [PacBio_runs.csv](PacBio_runs.csv): list of the PacBio runs used to call the variants.

   - [barcode_runs.csv](barcode_runs.csv): list of the Illumina runs used to count the barcodes for different samples.

   - [RBD_sites.csv](RBD_sites.csv): table giving indexing concordance between integer-numbering and spike residue numbering for each of the sarbecovirus backgrounds