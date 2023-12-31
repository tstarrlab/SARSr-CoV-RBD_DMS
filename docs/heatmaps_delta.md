---
layout: heatmaps_delta
permalink: /RBD-heatmaps_delta/
---

---

*NOTE data are preliminary*

### Overview

You can use this tool to explore the experimentally determined impacts of amino acid mutations on ACE2 binding affinity (delta-log10Ka) in sarbecovirus receptor-binding domains (RBD). 

#### Instruction

To use this tool, select the sarbecovirus variants and the metrics that you wish to display in the heatmap (change in ACE2 binding affinity (delta-log10Ka), or change in mean fluorescence intensity (delta-MFI) for RBD expression, by selecting that metric in the corresponding drop down menu. Hover over individual mutations to see exact numerical details. Click and drag on the site zoom bar to make it easier to scroll along the linear sequence.

*Note that sites are listed according to aligned SARS-CoV-2 (Wuhan-Hu-1) spike indexing for all sarbecoviruses, instead of their own indices. Columns also do not align perfectly (though tiptools/hover cursor does align to corresponding site)*

#### Technical Details

The impact on ACE2 receptor-binding of every single amino-acid mutation in six sarbecovirus RBDs, as determined by high-throughput FACS-seq assays. Wildtype amino acids are indicated by an 'x', and gray squares indicate missing mutations from each library. The number of internally replicated barcodes with which a mutation was measured is visible as `Barcode Count` in the tooltips, where higher numbers indicate higher-confidence measurements.


### Data

Raw data  can be found [here](https://github.com/tstarrlab/SARSr-CoV-RBD_DMS/blob/main/results/final_variant_scores/final_variant_scores_lib40_41.csv). The code used to make these plots can be found [here](https://github.com/tstarrlab/SARSr-CoV-RBD_DMS/blob/main/RBD-Heatmaps-Interactive-Visualization_delta.ipynb).
