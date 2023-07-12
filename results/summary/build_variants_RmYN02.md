<h1>Table of Contents<span class="tocSkip"></span></h1>
<div class="toc"><ul class="toc-item"><li><span><a href="#Process-CCSs" data-toc-modified-id="Process-CCSs-1">Process CCSs</a></span><ul class="toc-item"><li><span><a href="#Setup" data-toc-modified-id="Setup-1.1">Setup</a></span></li><li><span><a href="#PacBio-amplicons" data-toc-modified-id="PacBio-amplicons-1.2">PacBio amplicons</a></span></li><li><span><a href="#CCS-stats-for-PacBio-runs" data-toc-modified-id="CCS-stats-for-PacBio-runs-1.3">CCS stats for PacBio runs</a></span></li><li><span><a href="#Align-CCSs-to-amplicons" data-toc-modified-id="Align-CCSs-to-amplicons-1.4">Align CCSs to amplicons</a></span></li><li><span><a href="#Write-valid-CCSs" data-toc-modified-id="Write-valid-CCSs-1.5">Write valid CCSs</a></span></li></ul></li></ul></div>

# Build Variants
This Python Jupyter notebook processes the nt-level barcode lookup table for the RshSTT182 background to generate a protein-level barcode-variant lookup table and analyze mutation coverage.

## Setup

Import Python modules

Plotting is done with [plotnine](https://plotnine.readthedocs.io/en/stable/), which uses ggplot2-like syntax.

The analysis uses the Bloom lab's [alignparse](https://jbloomlab.github.io/alignparse) and [dms_variants](https://jbloomlab.github.io/dms_variants) packages.


```python
import collections
import math
import os
import re
import time
import warnings

import alignparse
import alignparse.ccs
from alignparse.constants import CBPALETTE
import alignparse.minimap2
import alignparse.targets
import alignparse.consensus

import dms_variants
import dms_variants.codonvarianttable
import dms_variants.plotnine_themes
import dms_variants.utils

from IPython.display import display, HTML

import numpy

import pandas as pd

from plotnine import *

import yaml

```

Define the background for variant mapping.


```python
background = "RmYN02"
```

Set [plotnine](https://plotnine.readthedocs.io/en/stable/) theme to the one defined in [dms_variants](https://jbloomlab.github.io/dms_variants):


```python
theme_set(dms_variants.plotnine_themes.theme_graygrid())
```

Versions of key software:


```python
print(f"Using alignparse version {alignparse.__version__}")
print(f"Using dms_variants version {dms_variants.__version__}")
```

    Using alignparse version 0.2.4
    Using dms_variants version 0.8.9


Ignore warnings that clutter output:


```python
warnings.simplefilter('ignore')
```

Read the configuration file:


```python
with open('config.yaml') as f:
    config = yaml.safe_load(f)
```

Make output directory for figures:


```python
os.makedirs(config['figs_dir'], exist_ok=True)
os.makedirs(config['variants_dir'], exist_ok=True)
```

## Consensus sequences
Read the CSV file giving consensus sequences (in terms of nt mutations) and barcodes. Subset on the background analyzed in this notebook:


```python
consensus = pd.read_csv(config['nt_variant_table_file'], na_filter=None)

consensus = consensus[consensus['target'] == background]
consensus = consensus[consensus['library'] != "lib46"]

nlibs = consensus['library'].nunique()  # number of unique libraries

ntargets = consensus['target'].nunique()  # number of unique targets

print(f"Read {len(consensus)} consensus sequences from {nlibs} libraries and {ntargets} targets.")

#output bg-specific nt table for reading in as codon variant table
(consensus
 [['target', 'library', 'barcode', 'substitutions', 'variant_call_support','number_of_indels']]
 .to_csv(config['nt_variant_table_file_' + background], index=False)
 )

```

    Read 56669 consensus sequences from 2 libraries and 1 targets.


## Create barcode-variant table
We now create a [CodonVariantTable](https://jbloomlab.github.io/dms_variants/dms_variants.codonvarianttable.html#dms_variants.codonvarianttable.CodonVariantTable) that stores and processes all the information about the variant consensus sequences.
Below we initialize such a table, and then analyze information about its composition.

### Initialize codon variant table
In order to initialize the codon variant table, we need two pieces of information:
  1. The wildtype gene sequence.
  2. The list of nucleotide mutations for each variant as determined in the consensus calling above.

Read "wildtype" gene sequence to which we made the alignments (in order to do this, initialize an `alignparse.Targets` and get the gene sequence from it):


```python
targets = alignparse.targets.Targets(seqsfile=config['amplicons_' + background],
                                     feature_parse_specs=config['feature_parse_specs_' + background])
geneseq = targets.get_target(background).get_feature('gene').seq

print(f"Read gene of {len(geneseq)} nts for {background} from {config['amplicons_' + background]}")
```

    Read gene of 546 nts for RmYN02 from data/PacBio_amplicon_RmYN02.gb


Now initialize the codon variant table using this wildtype sequence and our list of nucleotide mutations for each variant:


```python
variants = dms_variants.codonvarianttable.CodonVariantTable(
                barcode_variant_file=config['nt_variant_table_file_' + background],
                geneseq=geneseq,
                primary_target=background,
                )
```

### Basic stats on variants
We now will analyze the variants.
In this call and in the plots below, we set `samples=None` as we aren't looking at variant counts in specific samples, but are simply looking at properties of the variants in the table.

Here are the number of variants for each target:


```python
display(HTML(
    variants
    .n_variants_df(samples=None)
    .pivot_table(index=['target'],
                 columns='library',
                 values='count')
    .to_html()
    ))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>lib40</th>
      <th>lib41</th>
      <th>all libraries</th>
    </tr>
    <tr>
      <th>target</th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>RmYN02</th>
      <td>24980</td>
      <td>31689</td>
      <td>56669</td>
    </tr>
  </tbody>
</table>


Plot the number of variants supported by each number of CCSs:


```python
max_support = 10  # group variants with >= this much support

p = variants.plotVariantSupportHistogram(max_support=max_support,
                                         widthscale=1.1,
                                         heightscale=0.9)
p = p + theme(panel_grid_major_x=element_blank())  # no vertical grid lines
_ = p.draw()
```


    
![png](build_variants_RmYN02_files/build_variants_RmYN02_25_0.png)
    


### Mutations per variant
Plot the number of barcoded variants with each number of amino-acid and codon mutations.
This is for the primary target only, and doesn't include the spiked-in secondary targets:


```python
max_muts = 7  # group all variants with >= this many mutations

for mut_type in ['aa', 'codon']:
    p = variants.plotNumMutsHistogram(mut_type, samples=None, max_muts=max_muts,
                                      widthscale=1.1,
                                      heightscale=0.9)
    p = p + theme(panel_grid_major_x=element_blank())  # no vertical grid lines
    _ = p.draw()
    plotfile = os.path.join(config['figs_dir'], f"n_{mut_type}_muts_per_variant_"+background+".pdf")
    print(f"Saving plot to {plotfile}")
    p.save(plotfile)
```

    Saving plot to results/figures/n_aa_muts_per_variant_RmYN02.pdf
    Saving plot to results/figures/n_codon_muts_per_variant_RmYN02.pdf



    
![png](build_variants_RmYN02_files/build_variants_RmYN02_27_1.png)
    



    
![png](build_variants_RmYN02_files/build_variants_RmYN02_27_2.png)
    


Plot the frequencies of different codon mutation types among **all** variants (any number of mutations), again only for primary target:


```python
p = variants.plotNumCodonMutsByType(variant_type='all', samples=None,
                                    ylabel='mutations per variant',
                                    heightscale=0.8)
p = p + theme(panel_grid_major_x=element_blank())  # no vertical grid lines
_ = p.draw()
plotfile = os.path.join(config['figs_dir'], f"avg_muts_per_variant_"+background+".pdf")
print(f"Saving plot to {plotfile}")
p.save(plotfile)
```

    Saving plot to results/figures/avg_muts_per_variant_RmYN02.pdf



    
![png](build_variants_RmYN02_files/build_variants_RmYN02_29_1.png)
    


Variants supported by multiple PacBio CCSs should have fewer spurious mutations since sequencing errors are very unlikely to occur on two CCSs.
Below we plot the number of codon mutations per variant among variants with at least two CCSs supporting their call.
The difference in mutation rates here and in the plot above (that does not apply the `min_support=2` filter) gives some estimate of the frequency of mutations in our variants our spurious.
In fact, we see the numbers are very similar, indicating that few of the mutations are spurious:


```python
p = variants.plotNumCodonMutsByType(variant_type='all', samples=None,
                                    ylabel='mutations per variant', 
                                    min_support=2, heightscale=0.8)
p = p + theme(panel_grid_major_x=element_blank())  # no vertical grid lines
_ = p.draw()
```


    
![png](build_variants_RmYN02_files/build_variants_RmYN02_31_0.png)
    


### Completeness of mutation sampling
We examine how completely amino-acid mutations are sampled by the variants for the primary target, looking at single-mutant variants only and all variants.
The plot below shows that virtually every mutation is found in a variant in each library, even if we just look among the single mutants.
Things look especially good if we aggregate across libraries:


```python
for variant_type in ['all', 'single']:
    p = variants.plotCumulMutCoverage(variant_type, mut_type='aa', samples=None)
    _ = p.draw()
    plotfile = os.path.join(config['figs_dir'],
                            f"variant_cumul_{variant_type}_mut_coverage_"+background+".pdf")
    print(f"Saving plot to {plotfile}")
    p.save(plotfile)
```

    Saving plot to results/figures/variant_cumul_all_mut_coverage_RmYN02.pdf
    Saving plot to results/figures/variant_cumul_single_mut_coverage_RmYN02.pdf



    
![png](build_variants_RmYN02_files/build_variants_RmYN02_33_1.png)
    



    
![png](build_variants_RmYN02_files/build_variants_RmYN02_33_2.png)
    


To get more quantitative information like that plotted above, we determine how many mutations are found 0, 1, or >1 times both among single and all mutants for the primary target:


```python
count_dfs = []
for variant_type in ['all', 'single']:
    i_counts = (variants.mutCounts(variant_type, mut_type='aa', samples=None)
                .assign(variant_type=variant_type)
                )
    count_dfs += [i_counts.assign(include_stops=True),
                  i_counts
                  .query('not mutation.str.contains("\*")', engine='python')
                  .assign(include_stops=False)
                  ]
    
display(HTML(
    pd.concat(count_dfs)
    .assign(count=lambda x: (numpy.clip(x['count'], None, 2)
                             .map({0: '0', 1: '1', 2:'>1'}))
            )
    .groupby(['variant_type', 'include_stops', 'library', 'count'])
    .aggregate(number_of_mutations=pd.NamedAgg(column='mutation', aggfunc='count'))
    .to_html()
    ))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th>number_of_mutations</th>
    </tr>
    <tr>
      <th>variant_type</th>
      <th>include_stops</th>
      <th>library</th>
      <th>count</th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th rowspan="18" valign="top">all</th>
      <th rowspan="9" valign="top">False</th>
      <th rowspan="3" valign="top">lib40</th>
      <th>0</th>
      <td>64</td>
    </tr>
    <tr>
      <th>1</th>
      <td>173</td>
    </tr>
    <tr>
      <th>&gt;1</th>
      <td>3221</td>
    </tr>
    <tr>
      <th rowspan="3" valign="top">lib41</th>
      <th>0</th>
      <td>39</td>
    </tr>
    <tr>
      <th>1</th>
      <td>87</td>
    </tr>
    <tr>
      <th>&gt;1</th>
      <td>3332</td>
    </tr>
    <tr>
      <th rowspan="3" valign="top">all libraries</th>
      <th>0</th>
      <td>14</td>
    </tr>
    <tr>
      <th>1</th>
      <td>30</td>
    </tr>
    <tr>
      <th>&gt;1</th>
      <td>3414</td>
    </tr>
    <tr>
      <th rowspan="9" valign="top">True</th>
      <th rowspan="3" valign="top">lib40</th>
      <th>0</th>
      <td>223</td>
    </tr>
    <tr>
      <th>1</th>
      <td>187</td>
    </tr>
    <tr>
      <th>&gt;1</th>
      <td>3230</td>
    </tr>
    <tr>
      <th rowspan="3" valign="top">lib41</th>
      <th>0</th>
      <td>192</td>
    </tr>
    <tr>
      <th>1</th>
      <td>101</td>
    </tr>
    <tr>
      <th>&gt;1</th>
      <td>3347</td>
    </tr>
    <tr>
      <th rowspan="3" valign="top">all libraries</th>
      <th>0</th>
      <td>161</td>
    </tr>
    <tr>
      <th>1</th>
      <td>45</td>
    </tr>
    <tr>
      <th>&gt;1</th>
      <td>3434</td>
    </tr>
    <tr>
      <th rowspan="18" valign="top">single</th>
      <th rowspan="9" valign="top">False</th>
      <th rowspan="3" valign="top">lib40</th>
      <th>0</th>
      <td>103</td>
    </tr>
    <tr>
      <th>1</th>
      <td>215</td>
    </tr>
    <tr>
      <th>&gt;1</th>
      <td>3140</td>
    </tr>
    <tr>
      <th rowspan="3" valign="top">lib41</th>
      <th>0</th>
      <td>54</td>
    </tr>
    <tr>
      <th>1</th>
      <td>107</td>
    </tr>
    <tr>
      <th>&gt;1</th>
      <td>3297</td>
    </tr>
    <tr>
      <th rowspan="3" valign="top">all libraries</th>
      <th>0</th>
      <td>21</td>
    </tr>
    <tr>
      <th>1</th>
      <td>41</td>
    </tr>
    <tr>
      <th>&gt;1</th>
      <td>3396</td>
    </tr>
    <tr>
      <th rowspan="9" valign="top">True</th>
      <th rowspan="3" valign="top">lib40</th>
      <th>0</th>
      <td>275</td>
    </tr>
    <tr>
      <th>1</th>
      <td>221</td>
    </tr>
    <tr>
      <th>&gt;1</th>
      <td>3144</td>
    </tr>
    <tr>
      <th rowspan="3" valign="top">lib41</th>
      <th>0</th>
      <td>230</td>
    </tr>
    <tr>
      <th>1</th>
      <td>108</td>
    </tr>
    <tr>
      <th>&gt;1</th>
      <td>3302</td>
    </tr>
    <tr>
      <th rowspan="3" valign="top">all libraries</th>
      <th>0</th>
      <td>192</td>
    </tr>
    <tr>
      <th>1</th>
      <td>47</td>
    </tr>
    <tr>
      <th>&gt;1</th>
      <td>3401</td>
    </tr>
  </tbody>
</table>


### Mutation frequencies along gene
We plot the frequencies of mutations along the gene among the variants for the primary target.
Ideally, this would be uniform.
We make the plot for both all variants and single-mutant / wildtype variants:


```python
for variant_type in ['all', 'single']:
    p = variants.plotMutFreqs(variant_type, mut_type='codon', samples=None)
    p.draw()
```


    
![png](build_variants_RmYN02_files/build_variants_RmYN02_37_0.png)
    



    
![png](build_variants_RmYN02_files/build_variants_RmYN02_37_1.png)
    


We can also use heat maps to examine the extent to which specific amino-acid or codon mutations are over-represented.
These heat maps are large, so we make them just for all variants and the merge of all libraries:


```python
for mut_type in ['aa', 'codon']:
    p = variants.plotMutHeatmap('all', mut_type, samples=None, #libraries='all_only',
                                widthscale=2)
    p.draw()
```


    
![png](build_variants_RmYN02_files/build_variants_RmYN02_39_0.png)
    



    
![png](build_variants_RmYN02_files/build_variants_RmYN02_39_1.png)
    


### Write codon-variant table
We write the codon variant table to a CSV file.
This table looks like this:


```python
display(HTML(
    variants.barcode_variant_df
    .head()
    .to_html(index=False)
    ))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>target</th>
      <th>library</th>
      <th>barcode</th>
      <th>variant_call_support</th>
      <th>codon_substitutions</th>
      <th>aa_substitutions</th>
      <th>n_codon_substitutions</th>
      <th>n_aa_substitutions</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>RmYN02</td>
      <td>lib40</td>
      <td>AAAAAAAAACCTATCT</td>
      <td>4</td>
      <td>TCT181AAT</td>
      <td>S181N</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <td>RmYN02</td>
      <td>lib40</td>
      <td>AAAAAAAAATTAAAAT</td>
      <td>1</td>
      <td>ACA3TGT TCC134ATG TTG154TTT</td>
      <td>T3C S134M L154F</td>
      <td>3</td>
      <td>3</td>
    </tr>
    <tr>
      <td>RmYN02</td>
      <td>lib40</td>
      <td>AAAAAAAACACGATTA</td>
      <td>1</td>
      <td>GAT112CAT</td>
      <td>D112H</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <td>RmYN02</td>
      <td>lib40</td>
      <td>AAAAAAAACCTCGTAT</td>
      <td>3</td>
      <td></td>
      <td></td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <td>RmYN02</td>
      <td>lib40</td>
      <td>AAAAAAAATCTAACGT</td>
      <td>1</td>
      <td>GTT77ACT</td>
      <td>V77T</td>
      <td>1</td>
      <td>1</td>
    </tr>
  </tbody>
</table>


Note how this table differs from the nucleotide variant table we generated above and used to initialize the [CodonVariantTable](https://jbloomlab.github.io/dms_variants/dms_variants.codonvarianttable.html#dms_variants.codonvarianttable.CodonVariantTable) in that it gives **codon** substitutions and associated amino-acid substitutions.

Write it to CSV file:


```python
print(f"Writing codon-variant table to {config['codon_variant_table_file_'+background]}")

variants.barcode_variant_df.to_csv(config['codon_variant_table_file_'+background], index=False)
```

    Writing codon-variant table to results/variants/codon_variant_table_RmYN02.csv



```python

```
