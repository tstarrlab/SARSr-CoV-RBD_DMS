<h1>Table of Contents<span class="tocSkip"></span></h1>
<div class="toc"><ul class="toc-item"><li><span><a href="#Process-CCSs" data-toc-modified-id="Process-CCSs-1">Process CCSs</a></span><ul class="toc-item"><li><span><a href="#Setup" data-toc-modified-id="Setup-1.1">Setup</a></span></li><li><span><a href="#PacBio-amplicons" data-toc-modified-id="PacBio-amplicons-1.2">PacBio amplicons</a></span></li><li><span><a href="#CCS-stats-for-PacBio-runs" data-toc-modified-id="CCS-stats-for-PacBio-runs-1.3">CCS stats for PacBio runs</a></span></li><li><span><a href="#Align-CCSs-to-amplicons" data-toc-modified-id="Align-CCSs-to-amplicons-1.4">Align CCSs to amplicons</a></span></li><li><span><a href="#Write-valid-CCSs" data-toc-modified-id="Write-valid-CCSs-1.5">Write valid CCSs</a></span></li></ul></li></ul></div>

# Process CCSs
This Python Jupyter notebook processes the PacBio circular consensus sequences (CCSs) to extract barcodes and call mutations in the gene.

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
import dms_variants.plotnine_themes
import dms_variants.utils

from IPython.display import display, HTML

import numpy

import pandas as pd

from plotnine import *

import yaml

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
os.makedirs(config['process_ccs_dir'], exist_ok=True)
os.makedirs(config['variants_dir'], exist_ok=True)
```

## PacBio amplicons
Get the amplicons sequenced by PacBio as the alignment target along with the specs on how to parse the features:


```python
print(f"Reading amplicons from {config['amplicons']}")
print(f"Reading feature parse specs from {config['feature_parse_specs']}")

targets = alignparse.targets.Targets(
                seqsfile=config['amplicons'],
                feature_parse_specs=config['feature_parse_specs'])
```

    Reading amplicons from data/PacBio_amplicons.gb
    Reading feature parse specs from data/feature_parse_specs.yaml


Draw the target amplicons:


```python
fig = targets.plot(ax_width=7,
                   plots_indexing='biopython',  # numbering starts at 0
                   ax_height=2,  # height of each plot
                   hspace=1.2,  # vertical space between plots
                   )

plotfile = os.path.join(config['figs_dir'], 'amplicons.pdf')
print(f"Saving plot to {plotfile}")
fig.savefig(plotfile, bbox_inches='tight')
```

    Saving plot to results/figures/amplicons.pdf



    
![png](process_ccs_files/process_ccs_17_1.png)
    


Write out the specs used to parse the features (these are the same specs provided as `feature_parse_specs` when initializing `targets`, but with defaults filled in):


```python
print(targets.feature_parse_specs('yaml'))
```

    SARS-CoV-2_WH1: &id001
      query_clip5: 4
      query_clip3: 4
      termini5:
        filter:
          clip5: 4
          mutation_nt_count: 1
          mutation_op_count: null
          clip3: 0
        return: []
      gene:
        filter:
          mutation_nt_count: 18
          mutation_op_count: null
          clip5: 0
          clip3: 0
        return:
        - mutations
        - accuracy
      spacer:
        filter:
          mutation_nt_count: 1
          mutation_op_count: null
          clip5: 0
          clip3: 0
        return: []
      barcode:
        filter:
          mutation_nt_count: 0
          mutation_op_count: null
          clip5: 0
          clip3: 0
        return:
        - sequence
        - accuracy
      termini3:
        filter:
          clip3: 4
          mutation_nt_count: 1
          mutation_op_count: null
          clip5: 0
        return: []
    AncSarbecovirus_MAP: *id001
    AncSarbecovirus_alt: *id001
    AncSarbecovirus_alt1_ins117ins118: *id001
    AncSarbecovirus_tree1: *id001
    AncEurAf_alt: *id001
    AncEurAf_tree1: *id001
    AncAsia_MAP: *id001
    AncAsia_alt: *id001
    AncAsia_tree1: *id001
    AncAsia_tree2: *id001
    AncClade2_MAP: *id001
    AncClade2_alt1_subs-only: *id001
    AncClade2_alt2_del1-only: *id001
    AncClade2_alt3_del2-only: *id001
    AncClade2_alt4_dels-only: *id001
    AncClade2_alt: *id001
    AncClade2_tree2: *id001
    AncSARS2a_MAP: *id001
    AncSARS2a_alt: *id001
    AncSARS2c_MAP: *id001
    AncSARS1a_MAP: *id001
    AncSARS1a_alt: *id001
    AncSARS1a_tree1: *id001
    AncSARS1a_tree2: *id001
    AncSARS1c_MAP: *id001
    AncSARS-CoV-1_MAP: *id001
    AncSARS-CoV-1_alt: *id001
    GD-Pangolin: *id001
    RaTG13: *id001
    GX-Pangolin: *id001
    SARS-CoV-1_Urbani_HP03L: *id001
    SARS-CoV-1_HGZ8L1-A_HP03E: *id001
    SARS-CoV-1_GD01_HP03L: *id001
    SARS-CoV-1_GZ-C_HP03L: *id001
    SARS-CoV-1_Sino1-11_HP03L: *id001
    SARS-CoV-1_Sin852_HP03L: *id001
    SARS-CoV-1_SZ1_PC03: *id001
    SARS-CoV-1_GD03T0013_HP04: *id001
    SARS-CoV-1_GZ0402_HP04: *id001
    SARS-CoV-1_PC4-127_PC04: *id001
    SARS-CoV-1_PC4-137_PC04: *id001
    SARS-CoV-1_PC4-13_PC04: *id001
    WIV1: *id001
    Rs7327: *id001
    LYRa11: *id001
    Rs4231: *id001
    RsSHC014: *id001
    Rs4084: *id001
    BM48-31: *id001
    BtKY72: *id001
    ZXC21: *id001
    ZC45: *id001
    JL2012: *id001
    Rf1: *id001
    HeB2013: *id001
    273-2005: *id001
    Rf4092: *id001
    YN2013: *id001
    RmYN02: *id001
    As6526: *id001
    Rs4237: *id001
    Rs4081: *id001
    Rp3: *id001
    279-2005: *id001
    Shaanxi2011: *id001
    Yunnan2011: *id001
    Rs4247: *id001
    HKU3-1: *id001
    GX2013: *id001
    Longquan-140: *id001
    HKU3-8: *id001
    HuB2013: *id001
    SARS-CoV-2_2649: *id001
    SARS-CoV-1_2693: *id001
    Rc-o319: *id001
    RshSTT182: *id001
    PRD-0038: *id001
    PDF-2370: *id001
    RsYN04: *id001
    BB9904: *id001
    Khosta-1: *id001
    Khosta-2: *id001
    RhGB01: *id001
    BANAL-20-103: *id001
    BANAL-20-236: *id001
    YN2020B: *id001
    SARS-CoV-2_Alpha: *id001
    SARS-CoV-2_Beta: *id001
    SARS-CoV-2_Delta: *id001
    


## CCS stats for PacBio runs
Read data frame with information on PacBio runs:


```python
pacbio_runs = (
    pd.read_csv(config['pacbio_runs'], dtype=str)
    .drop(columns=['input_file','fastq_prebuilt'])
    .assign(name=lambda x: x['library'] + '_' + x['run'],
            fastq=lambda x: config['ccs_dir'] + '/' + x['name'] + '_ccs.fastq.gz'
            )
    )

display(HTML(pacbio_runs.to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>run</th>
      <th>name</th>
      <th>fastq</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>lib40</td>
      <td>210917</td>
      <td>lib40_210917</td>
      <td>results/ccs/lib40_210917_ccs.fastq.gz</td>
    </tr>
    <tr>
      <td>lib41</td>
      <td>210917</td>
      <td>lib41_210917</td>
      <td>results/ccs/lib41_210917_ccs.fastq.gz</td>
    </tr>
    <tr>
      <td>lib40</td>
      <td>211118_A</td>
      <td>lib40_211118_A</td>
      <td>results/ccs/lib40_211118_A_ccs.fastq.gz</td>
    </tr>
    <tr>
      <td>lib40</td>
      <td>211118_B</td>
      <td>lib40_211118_B</td>
      <td>results/ccs/lib40_211118_B_ccs.fastq.gz</td>
    </tr>
    <tr>
      <td>lib40</td>
      <td>211118_C</td>
      <td>lib40_211118_C</td>
      <td>results/ccs/lib40_211118_C_ccs.fastq.gz</td>
    </tr>
    <tr>
      <td>lib40</td>
      <td>211118_D</td>
      <td>lib40_211118_D</td>
      <td>results/ccs/lib40_211118_D_ccs.fastq.gz</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>200916_1</td>
      <td>lib46_200916_1</td>
      <td>results/ccs/lib46_200916_1_ccs.fastq.gz</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>200916_2</td>
      <td>lib46_200916_2</td>
      <td>results/ccs/lib46_200916_2_ccs.fastq.gz</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>201110</td>
      <td>lib46_201110</td>
      <td>results/ccs/lib46_201110_ccs.fastq.gz</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>211108_1</td>
      <td>lib46_211108_1</td>
      <td>results/ccs/lib46_211108_1_ccs.fastq.gz</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>211108_2</td>
      <td>lib46_211108_2</td>
      <td>results/ccs/lib46_211108_2_ccs.fastq.gz</td>
    </tr>
  </tbody>
</table>


Create an object that summarizes the `ccs` runs:


```python
ccs_summaries = alignparse.ccs.Summaries(pacbio_runs,
                                         report_col=None,
                                         ncpus=config['max_cpus'],
                                         )
```

If available, plot statistics on the number of ZMWs for each run:


```python
if ccs_summaries.has_zmw_stats():
    p = ccs_summaries.plot_zmw_stats()
    p = p + theme(panel_grid_major_x=element_blank())  # no vertical grid lines
    _ = p.draw()
else:
    print('No ZMW stats available.')
```

    No ZMW stats available.


Plot statistics on generated CCSs: their length, number of subread passes, and accuracy (as reported by the `ccs` program):


```python
for variable in ['length', 'passes', 'accuracy']:
    if ccs_summaries.has_stat(variable):
        p = ccs_summaries.plot_ccs_stats(variable, maxcol=7, bins=25)
        p = p + theme(panel_grid_major_x=element_blank())  # no vertical grid lines
        _ = p.draw()
    else:
        print(f"No {variable} statistics available.")
```

    No passes statistics available.



    
![png](process_ccs_files/process_ccs_27_1.png)
    



    
![png](process_ccs_files/process_ccs_27_2.png)
    


## Align CCSs to amplicons
We now align the CCSs to the amplicon and parse features from the resulting alignments using the specs above.

First, we initialize an `alignparse.minimap2.Mapper` to align the reads to SAM files:


```python
mapper = alignparse.minimap2.Mapper(alignparse.minimap2.OPTIONS_CODON_DMS)

print(f"Using `minimap2` {mapper.version} with these options:\n" +
      ' '.join(mapper.options))
```

    Using `minimap2` 2.18-r1015 with these options:
    -A2 -B4 -O12 -E2 --end-bonus=13 --secondary=no --cs


Next, we use `Targets.align_and_parse` to create the alignments and parse them:


```python
readstats, aligned, filtered = targets.align_and_parse(
        df=pacbio_runs,
        mapper=mapper,
        outdir=config['process_ccs_dir'],
        name_col='run',
        group_cols=['name', 'library'],
        queryfile_col='fastq',
        overwrite=True,
        ncpus=config['max_cpus'],
        )
```

First, examine the read stats from the alignment / parsing, both extracting alignment target name and getting stats aggregated by target:


```python
readstats = (
    readstats
    .assign(category_all_targets=lambda x: x['category'].str.split().str[0],
            target=lambda x: x['category'].str.split(None, 1).str[1],
            valid=lambda x: x['category_all_targets'] == 'aligned')
    )
```

Now plot the read stats by run (combining all targets and libraries within a run):


```python
ncol = 7
p = (
    ggplot(readstats
           .groupby(['name', 'category_all_targets', 'valid'])
           .aggregate({'count': 'sum'})
           .reset_index(),
           aes('category_all_targets', 'count', fill='valid')) +
    geom_bar(stat='identity') +
    facet_wrap('~ name', ncol=ncol) +
    theme(axis_text_x=element_text(angle=90),
          figure_size=(1.85 * min(ncol, len(pacbio_runs)),
                       2 * math.ceil(len(pacbio_runs) / ncol)),
          panel_grid_major_x=element_blank(),
          legend_position='none',
          ) +
    scale_fill_manual(values=CBPALETTE)
    )
_ = p.draw()
```


    
![png](process_ccs_files/process_ccs_35_0.png)
    


And the read stats by library (combining all targets and runs within a library):


```python
p = (
    ggplot(readstats
           .groupby(['library', 'category_all_targets', 'valid'])
           .aggregate({'count': 'sum'})
           .reset_index(), 
           aes('category_all_targets', 'count', fill='valid')) +
    geom_bar(stat='identity') +
    facet_wrap('~ library', nrow=1) +
    theme(axis_text_x=element_text(angle=90),
          figure_size=(1.5 * pacbio_runs['library'].nunique(), 2),
          panel_grid_major_x=element_blank(),
          legend_position='none',
          ) +
    scale_fill_manual(values=CBPALETTE)
    )
_ = p.draw()
```


    
![png](process_ccs_files/process_ccs_37_0.png)
    


And the number of reads by target (combining all libraries and runs for a target):


```python
p = (
    ggplot(readstats
           .groupby(['target'])
           .aggregate({'count': 'sum'})
           .reset_index(), 
           aes('count','target')) +
    geom_point(stat='identity', size=3) +
    theme(axis_text_x=element_text(angle=90),
          figure_size=(2, 0.3 * readstats['target'].nunique()),
          panel_grid_major_y=element_blank(),
          ) +
    scale_x_log10(name='number of reads')
    )
_ = p.draw()
```


    
![png](process_ccs_files/process_ccs_39_0.png)
    


And read stats by target (combining all libraries and runs for a target):


```python
p = (
    ggplot(readstats
           .groupby(['target', 'valid'])
           .aggregate({'count': 'sum'})
           .reset_index()
           .assign(total=lambda x: x.groupby('target')['count'].transform('sum'),
                   frac=lambda x: x['count'] / x['total'],
                   ), 
           aes('target', 'frac', fill='valid')) +
    geom_bar(stat='identity') +
    theme(axis_text_x=element_text(angle=90),
          figure_size=(0.5 * readstats['target'].nunique(), 2),
          panel_grid_major_x=element_blank(),
          ) +
    scale_fill_manual(values=CBPALETTE)
    )
_ = p.draw()
```


    
![png](process_ccs_files/process_ccs_41_0.png)
    


Now let's see **why** we filtered the reads.
First, we do some transformations on the `filtered` dict returned by `Targets.align_and_parse`.
Then we count up the number of CCSs filtered for each reason, and group together "unusual" reasons that represent less than some fraction of all filtering.
For now, we group together all targets to the stats represent all targets combined:


```python
other_cutoff = 0.02  # group as "other" reasons with <= this frac

filtered_df = (
    pd.concat(df.assign(target=target) for target, df in filtered.items())
    .groupby(['library', 'name', 'run', 'filter_reason'])
    .size()
    .rename('count')
    .reset_index()
    .assign(tot_reason_frac=lambda x: (x.groupby('filter_reason')['count']
                                       .transform('sum')) / x['count'].sum(),
            filter_reason=lambda x: numpy.where(x['tot_reason_frac'] > other_cutoff,
                                                x['filter_reason'],
                                                'other')
            )
    )
```

Now plot the filtering reason for all runs:


```python
ncol = 7
nreasons = filtered_df['filter_reason'].nunique()

p = (
    ggplot(filtered_df, aes('filter_reason', 'count')) +
    geom_bar(stat='identity') +
    facet_wrap('~ name', ncol=ncol) +
    theme(axis_text_x=element_text(angle=90),
          figure_size=(0.25 * nreasons * min(ncol, len(pacbio_runs)),
                       2 * math.ceil(len(pacbio_runs) / ncol)),
          panel_grid_major_x=element_blank(),
          )
    )
_ = p.draw()
```


    
![png](process_ccs_files/process_ccs_45_0.png)
    


Now make a similar plot to above, but combine all the runs for each library:


```python
p = (
    ggplot(filtered_df
           .groupby(['library', 'filter_reason'])
           .aggregate({'count': 'sum'})
           .reset_index(),
           aes('filter_reason', 'count')) +
    geom_bar(stat='identity') +
    facet_wrap('~ library', nrow=1) +
    theme(axis_text_x=element_text(angle=90),
          figure_size=(0.3 * nreasons * pacbio_runs['library'].nunique(), 2),
          panel_grid_major_x=element_blank(),
          )
    )
_ = p.draw()
```


    
![png](process_ccs_files/process_ccs_47_0.png)
    


Finally, we take the successfully parsed alignments and read them into a data frame, keeping track of the target that each CCS aligns to.
We also drop the pieces of information we won't use going forward, and rename a few columns:


```python
aligned_df = (
    pd.concat(df.assign(target=target) for target, df in aligned.items())
    .drop(columns=['query_clip5', 'query_clip3', 'run','name'])
    .rename(columns={'barcode_sequence': 'barcode'})
    )

print(f"First few lines of information on the parsed alignments:")
display(HTML(aligned_df.head().to_html(index=False)))
```

    First few lines of information on the parsed alignments:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>query_name</th>
      <th>gene_mutations</th>
      <th>gene_accuracy</th>
      <th>barcode</th>
      <th>barcode_accuracy</th>
      <th>target</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>lib46</td>
      <td>m54228_200925_184138/4194978/ccs</td>
      <td></td>
      <td>1.0</td>
      <td>GCTCCCTGGATGAATA</td>
      <td>1.0</td>
      <td>AncSarbecovirus_MAP</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>m54228_200925_184138/4195041/ccs</td>
      <td>T371C A372T</td>
      <td>1.0</td>
      <td>GCGGTTACCGCACAGG</td>
      <td>1.0</td>
      <td>AncSarbecovirus_MAP</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>m54228_200925_184138/4260338/ccs</td>
      <td>A499C C500A A501T ins506C</td>
      <td>1.0</td>
      <td>ACCATGAAGCCGACGG</td>
      <td>1.0</td>
      <td>AncSarbecovirus_MAP</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>m54228_200925_184138/4260429/ccs</td>
      <td>T370G T371G A372T</td>
      <td>1.0</td>
      <td>GCGTGACTAATATTAA</td>
      <td>1.0</td>
      <td>AncSarbecovirus_MAP</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>m54228_200925_184138/4260622/ccs</td>
      <td>T487A C488G T489A</td>
      <td>1.0</td>
      <td>ATACTAGATGAGATAA</td>
      <td>1.0</td>
      <td>AncSarbecovirus_MAP</td>
    </tr>
  </tbody>
</table>


## Write valid CCSs

Write the processed CCSs to a file:


```python
aligned_df.to_csv(config['processed_ccs_file'], index=False)

print("Barcodes and mutations for valid processed CCSs "
      f"have been written to {config['processed_ccs_file']}.")
```

    Barcodes and mutations for valid processed CCSs have been written to results/process_ccs/processed_ccs.csv.


# Build barcode variant table
Next, we build consensus sequences for barcoded variants from the mutations called in the processed PacBio CCSs.

## Build variant consensus sequences
We now use the processed CCSs to build a consensus sequence (list of mutations) for each barcoded variant in each library.

### Get processed CCSs
Read the CSV file with the processed CCSs into a data frame:


```python
processed_ccs = pd.read_csv(config['processed_ccs_file'], na_filter=None)

nlibs = processed_ccs['library'].nunique()  # number of unique libraries

ntargets = processed_ccs['target'].nunique()  # number of unique targets

print(f"Read {len(processed_ccs)} CCSs from {nlibs} libraries and {ntargets} targets.")
```

    Read 1237173 CCSs from 3 libraries and 90 targets.


Overall statistics on number of total CCSs and number of unique barcodes:


```python
display(HTML(
    processed_ccs
    .groupby(['target', 'library'])
    .aggregate(total_CCSs=('barcode', 'size'),
               unique_barcodes=('barcode', 'nunique'))
    .assign(avg_CCSs_per_barcode=lambda x: x['total_CCSs'] / x['unique_barcodes'])
    .round(2)
    .to_html()
    ))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th></th>
      <th>total_CCSs</th>
      <th>unique_barcodes</th>
      <th>avg_CCSs_per_barcode</th>
    </tr>
    <tr>
      <th>target</th>
      <th>library</th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>273-2005</th>
      <th>lib46</th>
      <td>1206</td>
      <td>436</td>
      <td>2.77</td>
    </tr>
    <tr>
      <th>279-2005</th>
      <th>lib46</th>
      <td>1361</td>
      <td>454</td>
      <td>3.00</td>
    </tr>
    <tr>
      <th>AncAsia_MAP</th>
      <th>lib46</th>
      <td>15531</td>
      <td>6165</td>
      <td>2.52</td>
    </tr>
    <tr>
      <th>AncAsia_alt</th>
      <th>lib46</th>
      <td>965</td>
      <td>321</td>
      <td>3.01</td>
    </tr>
    <tr>
      <th>AncAsia_tree1</th>
      <th>lib46</th>
      <td>621</td>
      <td>198</td>
      <td>3.14</td>
    </tr>
    <tr>
      <th>AncAsia_tree2</th>
      <th>lib46</th>
      <td>1192</td>
      <td>376</td>
      <td>3.17</td>
    </tr>
    <tr>
      <th>AncClade2_MAP</th>
      <th>lib46</th>
      <td>12903</td>
      <td>5246</td>
      <td>2.46</td>
    </tr>
    <tr>
      <th>AncClade2_alt</th>
      <th>lib46</th>
      <td>1464</td>
      <td>500</td>
      <td>2.93</td>
    </tr>
    <tr>
      <th>AncClade2_alt1_subs-only</th>
      <th>lib46</th>
      <td>1305</td>
      <td>448</td>
      <td>2.91</td>
    </tr>
    <tr>
      <th>AncClade2_alt2_del1-only</th>
      <th>lib46</th>
      <td>1438</td>
      <td>453</td>
      <td>3.17</td>
    </tr>
    <tr>
      <th>AncClade2_alt3_del2-only</th>
      <th>lib46</th>
      <td>1207</td>
      <td>389</td>
      <td>3.10</td>
    </tr>
    <tr>
      <th>AncClade2_alt4_dels-only</th>
      <th>lib46</th>
      <td>1204</td>
      <td>409</td>
      <td>2.94</td>
    </tr>
    <tr>
      <th>AncClade2_tree2</th>
      <th>lib46</th>
      <td>1561</td>
      <td>510</td>
      <td>3.06</td>
    </tr>
    <tr>
      <th>AncEurAf_alt</th>
      <th>lib46</th>
      <td>1331</td>
      <td>466</td>
      <td>2.86</td>
    </tr>
    <tr>
      <th>AncEurAf_tree1</th>
      <th>lib46</th>
      <td>1257</td>
      <td>394</td>
      <td>3.19</td>
    </tr>
    <tr>
      <th>AncSARS-CoV-1_MAP</th>
      <th>lib46</th>
      <td>603</td>
      <td>201</td>
      <td>3.00</td>
    </tr>
    <tr>
      <th>AncSARS-CoV-1_alt</th>
      <th>lib46</th>
      <td>1348</td>
      <td>411</td>
      <td>3.28</td>
    </tr>
    <tr>
      <th>AncSARS1a_MAP</th>
      <th>lib46</th>
      <td>10948</td>
      <td>4500</td>
      <td>2.43</td>
    </tr>
    <tr>
      <th>AncSARS1a_alt</th>
      <th>lib46</th>
      <td>1585</td>
      <td>494</td>
      <td>3.21</td>
    </tr>
    <tr>
      <th>AncSARS1a_tree1</th>
      <th>lib46</th>
      <td>872</td>
      <td>247</td>
      <td>3.53</td>
    </tr>
    <tr>
      <th>AncSARS1a_tree2</th>
      <th>lib46</th>
      <td>1323</td>
      <td>430</td>
      <td>3.08</td>
    </tr>
    <tr>
      <th>AncSARS1c_MAP</th>
      <th>lib46</th>
      <td>575</td>
      <td>182</td>
      <td>3.16</td>
    </tr>
    <tr>
      <th>AncSARS2a_MAP</th>
      <th>lib46</th>
      <td>14890</td>
      <td>5662</td>
      <td>2.63</td>
    </tr>
    <tr>
      <th>AncSARS2a_alt</th>
      <th>lib46</th>
      <td>714</td>
      <td>243</td>
      <td>2.94</td>
    </tr>
    <tr>
      <th>AncSARS2c_MAP</th>
      <th>lib46</th>
      <td>8904</td>
      <td>3715</td>
      <td>2.40</td>
    </tr>
    <tr>
      <th>AncSarbecovirus_MAP</th>
      <th>lib46</th>
      <td>12728</td>
      <td>5100</td>
      <td>2.50</td>
    </tr>
    <tr>
      <th>AncSarbecovirus_alt</th>
      <th>lib46</th>
      <td>581</td>
      <td>201</td>
      <td>2.89</td>
    </tr>
    <tr>
      <th>AncSarbecovirus_alt1_ins117ins118</th>
      <th>lib46</th>
      <td>597</td>
      <td>195</td>
      <td>3.06</td>
    </tr>
    <tr>
      <th>AncSarbecovirus_tree1</th>
      <th>lib46</th>
      <td>1409</td>
      <td>418</td>
      <td>3.37</td>
    </tr>
    <tr>
      <th>As6526</th>
      <th>lib46</th>
      <td>1044</td>
      <td>341</td>
      <td>3.06</td>
    </tr>
    <tr>
      <th>BANAL-20-103</th>
      <th>lib46</th>
      <td>1625</td>
      <td>351</td>
      <td>4.63</td>
    </tr>
    <tr>
      <th>BANAL-20-236</th>
      <th>lib46</th>
      <td>1544</td>
      <td>333</td>
      <td>4.64</td>
    </tr>
    <tr>
      <th>BB9904</th>
      <th>lib46</th>
      <td>1666</td>
      <td>410</td>
      <td>4.06</td>
    </tr>
    <tr>
      <th>BM48-31</th>
      <th>lib46</th>
      <td>12750</td>
      <td>4828</td>
      <td>2.64</td>
    </tr>
    <tr>
      <th>BtKY72</th>
      <th>lib46</th>
      <td>38226</td>
      <td>6215</td>
      <td>6.15</td>
    </tr>
    <tr>
      <th>GD-Pangolin</th>
      <th>lib46</th>
      <td>13048</td>
      <td>4583</td>
      <td>2.85</td>
    </tr>
    <tr>
      <th>GX-Pangolin</th>
      <th>lib46</th>
      <td>861</td>
      <td>271</td>
      <td>3.18</td>
    </tr>
    <tr>
      <th>GX2013</th>
      <th>lib46</th>
      <td>2240</td>
      <td>692</td>
      <td>3.24</td>
    </tr>
    <tr>
      <th>HKU3-1</th>
      <th>lib46</th>
      <td>822</td>
      <td>240</td>
      <td>3.42</td>
    </tr>
    <tr>
      <th>HKU3-8</th>
      <th>lib46</th>
      <td>1747</td>
      <td>610</td>
      <td>2.86</td>
    </tr>
    <tr>
      <th>HeB2013</th>
      <th>lib46</th>
      <td>1310</td>
      <td>467</td>
      <td>2.81</td>
    </tr>
    <tr>
      <th>HuB2013</th>
      <th>lib46</th>
      <td>1455</td>
      <td>505</td>
      <td>2.88</td>
    </tr>
    <tr>
      <th>JL2012</th>
      <th>lib46</th>
      <td>1773</td>
      <td>590</td>
      <td>3.01</td>
    </tr>
    <tr>
      <th>Khosta-1</th>
      <th>lib46</th>
      <td>1216</td>
      <td>287</td>
      <td>4.24</td>
    </tr>
    <tr>
      <th>Khosta-2</th>
      <th>lib46</th>
      <td>1800</td>
      <td>426</td>
      <td>4.23</td>
    </tr>
    <tr>
      <th>LYRa11</th>
      <th>lib46</th>
      <td>1221</td>
      <td>374</td>
      <td>3.26</td>
    </tr>
    <tr>
      <th>Longquan-140</th>
      <th>lib46</th>
      <td>1452</td>
      <td>442</td>
      <td>3.29</td>
    </tr>
    <tr>
      <th>PDF-2370</th>
      <th>lib46</th>
      <td>1364</td>
      <td>306</td>
      <td>4.46</td>
    </tr>
    <tr>
      <th rowspan="3" valign="top">PRD-0038</th>
      <th>lib40</th>
      <td>78047</td>
      <td>37109</td>
      <td>2.10</td>
    </tr>
    <tr>
      <th>lib41</th>
      <td>103824</td>
      <td>41846</td>
      <td>2.48</td>
    </tr>
    <tr>
      <th>lib46</th>
      <td>1878</td>
      <td>408</td>
      <td>4.60</td>
    </tr>
    <tr>
      <th>RaTG13</th>
      <th>lib46</th>
      <td>12354</td>
      <td>4305</td>
      <td>2.87</td>
    </tr>
    <tr>
      <th>Rc-o319</th>
      <th>lib46</th>
      <td>925</td>
      <td>253</td>
      <td>3.66</td>
    </tr>
    <tr>
      <th>Rf1</th>
      <th>lib46</th>
      <td>952</td>
      <td>325</td>
      <td>2.93</td>
    </tr>
    <tr>
      <th>Rf4092</th>
      <th>lib46</th>
      <td>1722</td>
      <td>580</td>
      <td>2.97</td>
    </tr>
    <tr>
      <th>RhGB01</th>
      <th>lib46</th>
      <td>1785</td>
      <td>411</td>
      <td>4.34</td>
    </tr>
    <tr>
      <th rowspan="3" valign="top">RmYN02</th>
      <th>lib40</th>
      <td>55458</td>
      <td>27348</td>
      <td>2.03</td>
    </tr>
    <tr>
      <th>lib41</th>
      <td>87001</td>
      <td>34000</td>
      <td>2.56</td>
    </tr>
    <tr>
      <th>lib46</th>
      <td>1382</td>
      <td>485</td>
      <td>2.85</td>
    </tr>
    <tr>
      <th>Rp3</th>
      <th>lib46</th>
      <td>1998</td>
      <td>596</td>
      <td>3.35</td>
    </tr>
    <tr>
      <th>Rs4081</th>
      <th>lib46</th>
      <td>1347</td>
      <td>445</td>
      <td>3.03</td>
    </tr>
    <tr>
      <th>Rs4084</th>
      <th>lib46</th>
      <td>1038</td>
      <td>323</td>
      <td>3.21</td>
    </tr>
    <tr>
      <th>Rs4231</th>
      <th>lib46</th>
      <td>1510</td>
      <td>515</td>
      <td>2.93</td>
    </tr>
    <tr>
      <th>Rs4237</th>
      <th>lib46</th>
      <td>1361</td>
      <td>475</td>
      <td>2.87</td>
    </tr>
    <tr>
      <th>Rs4247</th>
      <th>lib46</th>
      <td>2073</td>
      <td>696</td>
      <td>2.98</td>
    </tr>
    <tr>
      <th>Rs7327</th>
      <th>lib46</th>
      <td>9695</td>
      <td>3916</td>
      <td>2.48</td>
    </tr>
    <tr>
      <th>RsSHC014</th>
      <th>lib46</th>
      <td>894</td>
      <td>280</td>
      <td>3.19</td>
    </tr>
    <tr>
      <th rowspan="3" valign="top">RsYN04</th>
      <th>lib40</th>
      <td>80876</td>
      <td>38130</td>
      <td>2.12</td>
    </tr>
    <tr>
      <th>lib41</th>
      <td>104323</td>
      <td>42167</td>
      <td>2.47</td>
    </tr>
    <tr>
      <th>lib46</th>
      <td>1254</td>
      <td>294</td>
      <td>4.27</td>
    </tr>
    <tr>
      <th rowspan="3" valign="top">RshSTT182</th>
      <th>lib40</th>
      <td>97111</td>
      <td>39658</td>
      <td>2.45</td>
    </tr>
    <tr>
      <th>lib41</th>
      <td>96830</td>
      <td>40284</td>
      <td>2.40</td>
    </tr>
    <tr>
      <th>lib46</th>
      <td>1045</td>
      <td>272</td>
      <td>3.84</td>
    </tr>
    <tr>
      <th rowspan="3" valign="top">SARS-CoV-1_2693</th>
      <th>lib40</th>
      <td>110932</td>
      <td>40006</td>
      <td>2.77</td>
    </tr>
    <tr>
      <th>lib41</th>
      <td>108352</td>
      <td>36411</td>
      <td>2.98</td>
    </tr>
    <tr>
      <th>lib46</th>
      <td>23031</td>
      <td>5778</td>
      <td>3.99</td>
    </tr>
    <tr>
      <th>SARS-CoV-1_GD01_HP03L</th>
      <th>lib46</th>
      <td>784</td>
      <td>237</td>
      <td>3.31</td>
    </tr>
    <tr>
      <th>SARS-CoV-1_GD03T0013_HP04</th>
      <th>lib46</th>
      <td>1357</td>
      <td>483</td>
      <td>2.81</td>
    </tr>
    <tr>
      <th>SARS-CoV-1_GZ-C_HP03L</th>
      <th>lib46</th>
      <td>607</td>
      <td>195</td>
      <td>3.11</td>
    </tr>
    <tr>
      <th>SARS-CoV-1_GZ0402_HP04</th>
      <th>lib46</th>
      <td>1044</td>
      <td>334</td>
      <td>3.13</td>
    </tr>
    <tr>
      <th>SARS-CoV-1_HGZ8L1-A_HP03E</th>
      <th>lib46</th>
      <td>1017</td>
      <td>332</td>
      <td>3.06</td>
    </tr>
    <tr>
      <th>SARS-CoV-1_PC4-127_PC04</th>
      <th>lib46</th>
      <td>1603</td>
      <td>492</td>
      <td>3.26</td>
    </tr>
    <tr>
      <th>SARS-CoV-1_PC4-137_PC04</th>
      <th>lib46</th>
      <td>8929</td>
      <td>3313</td>
      <td>2.70</td>
    </tr>
    <tr>
      <th>SARS-CoV-1_PC4-13_PC04</th>
      <th>lib46</th>
      <td>1369</td>
      <td>441</td>
      <td>3.10</td>
    </tr>
    <tr>
      <th>SARS-CoV-1_SZ1_PC03</th>
      <th>lib46</th>
      <td>1244</td>
      <td>394</td>
      <td>3.16</td>
    </tr>
    <tr>
      <th>SARS-CoV-1_Sin852_HP03L</th>
      <th>lib46</th>
      <td>1238</td>
      <td>394</td>
      <td>3.14</td>
    </tr>
    <tr>
      <th>SARS-CoV-1_Sino1-11_HP03L</th>
      <th>lib46</th>
      <td>1385</td>
      <td>504</td>
      <td>2.75</td>
    </tr>
    <tr>
      <th>SARS-CoV-1_Urbani_HP03L</th>
      <th>lib46</th>
      <td>11520</td>
      <td>4339</td>
      <td>2.65</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">SARS-CoV-2_2649</th>
      <th>lib40</th>
      <td>1</td>
      <td>1</td>
      <td>1.00</td>
    </tr>
    <tr>
      <th>lib46</th>
      <td>2997</td>
      <td>899</td>
      <td>3.33</td>
    </tr>
    <tr>
      <th>SARS-CoV-2_Alpha</th>
      <th>lib46</th>
      <td>2539</td>
      <td>485</td>
      <td>5.24</td>
    </tr>
    <tr>
      <th>SARS-CoV-2_Beta</th>
      <th>lib46</th>
      <td>1865</td>
      <td>375</td>
      <td>4.97</td>
    </tr>
    <tr>
      <th>SARS-CoV-2_Delta</th>
      <th>lib46</th>
      <td>1972</td>
      <td>376</td>
      <td>5.24</td>
    </tr>
    <tr>
      <th>SARS-CoV-2_WH1</th>
      <th>lib46</th>
      <td>6764</td>
      <td>2519</td>
      <td>2.69</td>
    </tr>
    <tr>
      <th>Shaanxi2011</th>
      <th>lib46</th>
      <td>1479</td>
      <td>508</td>
      <td>2.91</td>
    </tr>
    <tr>
      <th>WIV1</th>
      <th>lib46</th>
      <td>1265</td>
      <td>391</td>
      <td>3.24</td>
    </tr>
    <tr>
      <th>YN2013</th>
      <th>lib46</th>
      <td>1162</td>
      <td>419</td>
      <td>2.77</td>
    </tr>
    <tr>
      <th>YN2020B</th>
      <th>lib46</th>
      <td>1732</td>
      <td>401</td>
      <td>4.32</td>
    </tr>
    <tr>
      <th>Yunnan2011</th>
      <th>lib46</th>
      <td>1658</td>
      <td>567</td>
      <td>2.92</td>
    </tr>
    <tr>
      <th>ZC45</th>
      <th>lib46</th>
      <td>1467</td>
      <td>468</td>
      <td>3.13</td>
    </tr>
    <tr>
      <th>ZXC21</th>
      <th>lib46</th>
      <td>2390</td>
      <td>774</td>
      <td>3.09</td>
    </tr>
  </tbody>
</table>


### Filter processed CCSs
We have the PacBio `ccs` program's estimated "accuracy" for both the barcode and the gene sequence for each processed CCS.
We will filter the CCSs to only keep ones of sufficiently high accuracy.

First, we want to plot the accuracies.
It is actually visually easier to look at the error rate, which is one minus the accuracy.
Because we want to plot on a log scale (which can't show error rates of zero), we define a *error_rate_floor*, and set all error rates less than this to that value:


```python
error_rate_floor = 1e-7  # error rates < this set to this
if error_rate_floor >= config['max_error_rate']:
    raise ValueError('error_rate_floor must be < max_error_rate')

processed_ccs = (
    processed_ccs
    .assign(barcode_error=lambda x: numpy.clip(1 - x['barcode_accuracy'],
                                               error_rate_floor, None),
            gene_error=lambda x: numpy.clip(1 - x['gene_accuracy'],
                                            error_rate_floor, None)
            )
    )
```

Now plot the error rates, drawing a dashed vertical line at the threshold separating the CCSs we retain for consensus building versus those that we discard:


```python
_ = (
 ggplot(processed_ccs
        .melt(value_vars=['barcode_error', 'gene_error'],
              var_name='feature_type', value_name='error rate'),
        aes('error rate')) +
 geom_histogram(bins=25) +
 geom_vline(xintercept=config['max_error_rate'],
            linetype='dashed',
            color=CBPALETTE[1]) +
 facet_wrap('~ feature_type') +
 theme(figure_size=(4.5, 2)) +
 ylab('number of CCSs') +
 scale_x_log10()
 ).draw()
```


    
![png](process_ccs_files/process_ccs_60_0.png)
    


Flag the CCSs to retain, and indicate how many we are retaining and purging due to the accuracy filter:


```python
processed_ccs = (
    processed_ccs
    .assign(retained=lambda x: ((x['gene_error'] < config['max_error_rate']) &
                                (x['barcode_error'] < config['max_error_rate'])))
    )
```

Here are number of retained CCSs:


```python
_ = (
 ggplot(processed_ccs.assign(xlabel=lambda x: x['target'] + ', ' + x['library'])
                     .groupby(['xlabel', 'retained'])
                     .size()
                     .rename('count')
                     .reset_index(),
        aes('count', 'xlabel', color='retained', label='count')) +
 geom_point(size=3) +
 geom_text(va='bottom', size=7, ha='center',format_string='{:.3g}', nudge_x=0.2) +
 theme(figure_size=(3, 0.5 * nlibs * ntargets),
       panel_grid_major_y=element_blank(),
       axis_text_x=element_text(angle=90),
       ) +
 scale_x_log10(name='number of CCSs') +
 ylab('') +
 scale_color_manual(values=CBPALETTE[1:])
 ).draw()
```


    
![png](process_ccs_files/process_ccs_64_0.png)
    


### Sequences per barcode
How many times is each barcode sequenced?
This is useful to know for thinking about building the barcode consensus.

First, plot the distribution of the number of times each **barcode** is observed among the retained CCSs:


```python
max_count = 8 # in plot, group all barcodes with >= this many counts

p = (
 ggplot(
    processed_ccs
     .query('retained')
     .groupby(['library', 'barcode'])
     .size()
     .rename('nseqs')
     .reset_index()
     .assign(nseqs=lambda x: numpy.clip(x['nseqs'], None, max_count)),
    aes('nseqs')) +
 geom_bar() +
 facet_wrap('~ library', nrow=1) +
 theme(figure_size=(1.75 * nlibs, 2),
       panel_grid_major_x=element_blank(),
       ) +
 ylab('number of barcodes') +
 xlab('CCSs for barcode')
 )

_ = p.draw()
```


    
![png](process_ccs_files/process_ccs_66_0.png)
    


### Empirical accuracy of CCSs
We want to directly estimate the accuracy of the gene-barcode link rather than relying on the PacBio `ccs` accuracy, which doesn't include inaccuracies due to things like strand exchange or the same barcode on different sequences.

One way to do this is to examine instances when we have multiple sequences for the same barcode. 
We can calculate the empirical accuracy of the sequences by looking at all instances of multiple sequences of the same barcode and determining how often they are identical.
This calculation is performed by `alignparse.consensus.empirical_accuracy` using the equations described in the docs for that function.

We will do this four for sets of sequences:

 1. All of the CCSs retained above.
 2. CCSs retained by applying a PacBio `ccs` accuracy filter 10-fold more stringent than the one above.
    The rationale is that if this improves the concordance (real accuracy) of the CCSs substantially then maybe we should make the accuracy filter more stringent.
 3. Like (1) but excluding all CCSs with indels.
    the rationale is that we only really care about substitutions, and will exclude sequences with indels anyway.
 4. Like (2) but excluding all CCSs with indels.
 
First, we annotate the sequences with the number of indels and whether they have an indel to enable categorization into the aforementioned sets:


```python
processed_ccs = alignparse.consensus.add_mut_info_cols(processed_ccs,
                                                       mutation_col='gene_mutations',
                                                       n_indel_col='n_indels')

processed_ccs = processed_ccs.assign(has_indel=lambda x: x['n_indels'] > 0)
```

Plot how many sequences have indels:


```python
_ = (
 ggplot(processed_ccs,
        aes('retained', fill='has_indel')) +
 geom_bar(position='dodge') +
 geom_text(aes(label='..count..'), stat='count', va='bottom', size=7,
           position=position_dodge(width=0.9), format_string='{:.2g}') +
 theme(figure_size=(2.5 * nlibs, 3),
       panel_grid_major_x=element_blank(),
       ) +
 ylab('number of CCSs') +
 scale_fill_manual(values=CBPALETTE[1:]) +
 facet_wrap('~ library', nrow=1)
 ).draw()
```


    
![png](process_ccs_files/process_ccs_70_0.png)
    


Now get the empirical accuracy for each of the CCS groups mentioned above:


```python
high_acc = config['max_error_rate'] / 10
empirical_acc = []

for desc, query_str in [
        ('retained', 'retained'),
        ('retained, no indel', 'retained and not has_indel'),
        ('10X accuracy',
         f"(gene_error < {high_acc}) and (barcode_error < {high_acc})"),
        ('10X accuracy, no indel',
         f"(gene_error < {high_acc}) and (barcode_error < {high_acc}) and not has_indel")
        ]:
    # get just CCSs in that category
    df = processed_ccs.query(query_str)
    
    # compute empirical accuracy
    empirical_acc.append(
        alignparse.consensus.empirical_accuracy(df,
                                                mutation_col='gene_mutations')
        .assign(description=desc)
        .merge(df
               .groupby('library')
               .size()
               .rename('number_CCSs')
               .reset_index()
               )
        )

# make description categorical to preserve order, and annotate as "actual"
# the category ("retained, no indel") that we will use for building variants.
empirical_acc = (
    pd.concat(empirical_acc, ignore_index=True, sort=False)
    .assign(description=lambda x: pd.Categorical(x['description'],
                                                 x['description'].unique(),
                                                 ordered=True),
            actual=lambda x: numpy.where(x['description'] == 'retained, no indel',
                                         True, False),
            )
    )
```

Display table of the empirical accuracies:


```python
display(HTML(empirical_acc.to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>accuracy</th>
      <th>description</th>
      <th>number_CCSs</th>
      <th>actual</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>lib40</td>
      <td>0.967651</td>
      <td>retained</td>
      <td>380810</td>
      <td>False</td>
    </tr>
    <tr>
      <td>lib41</td>
      <td>0.971843</td>
      <td>retained</td>
      <td>454965</td>
      <td>False</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>0.966077</td>
      <td>retained</td>
      <td>291422</td>
      <td>False</td>
    </tr>
    <tr>
      <td>lib40</td>
      <td>0.987410</td>
      <td>retained, no indel</td>
      <td>370289</td>
      <td>True</td>
    </tr>
    <tr>
      <td>lib41</td>
      <td>0.989833</td>
      <td>retained, no indel</td>
      <td>443069</td>
      <td>True</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>0.990398</td>
      <td>retained, no indel</td>
      <td>275046</td>
      <td>True</td>
    </tr>
    <tr>
      <td>lib40</td>
      <td>0.974973</td>
      <td>10X accuracy</td>
      <td>358171</td>
      <td>False</td>
    </tr>
    <tr>
      <td>lib41</td>
      <td>0.978484</td>
      <td>10X accuracy</td>
      <td>428290</td>
      <td>False</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>0.975930</td>
      <td>10X accuracy</td>
      <td>278066</td>
      <td>False</td>
    </tr>
    <tr>
      <td>lib40</td>
      <td>0.987744</td>
      <td>10X accuracy, no indel</td>
      <td>350743</td>
      <td>False</td>
    </tr>
    <tr>
      <td>lib41</td>
      <td>0.989879</td>
      <td>10X accuracy, no indel</td>
      <td>419905</td>
      <td>False</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>0.991242</td>
      <td>10X accuracy, no indel</td>
      <td>264329</td>
      <td>False</td>
    </tr>
  </tbody>
</table>


Plot the empirical accuracies, using a different color to show the category that we will actually use:


```python
p = (
    ggplot(empirical_acc,
           aes('description', 'accuracy', color='actual', label='accuracy')
           ) +
    geom_point(size=3) +
    geom_text(va='bottom', size=9, format_string='{:.3g}', nudge_y=0.003) +
    facet_wrap('~ library') +
    theme(figure_size=(1.75 * nlibs, 2.25),
          axis_text_x=element_text(angle=90),
          panel_grid_major_x=element_blank(),
          ) +
    xlab('') +
    scale_y_continuous(name='empirical accuracy', limits=(0.95, 1.005)) +
    scale_color_manual(values=CBPALETTE, guide=False)
    )

plotfile = os.path.join(config['figs_dir'], 'empirical_CCS_accuracy.pdf')
print(f"Saving plot to {plotfile}")
_ = p.draw()
```

    Saving plot to results/figures/empirical_CCS_accuracy.pdf



    
![png](process_ccs_files/process_ccs_76_1.png)
    


The above analysis shows that if we exclude sequences with indels (which we plan to do among our consensus sequences), then the accuracy of each CCS is around 99%. 
We do **not** get notably higher empirical accuracy by imposing a more stringent filter from the PacBio `ccs` program, indicating that the major sources of error are due to processes that are not modeled in this program's accuracy filter (perhaps strand exchange or barcode sharing).

Note that this empirical accuracy is for a **single** CCS.
When we build the consensus sequences for each barcode below, we will take the consensus of CCSs within a barcode.
So for barcodes with multiple CCSs, the actual accuracy of the consensus sequences will be higher than the empirical accuracy above due to capturing information from multiple CCSs.

### Consensus sequences for barcodes
We call the consensus sequence for each barcode using the simple method implemented in [alignparse.consensus.simple_mutconsensus](https://jbloomlab.github.io/alignparse/alignparse.consensus.html?highlight=simple_mutconsensus#alignparse.consensus.simple_mutconsensus).
The documentation for that function explains the method in detail, but basically it works like this:
 1. When there is just one CCS per barcode, the consensus is just that sequence.
 2. When there are multiple CCSs per barcode, they are used to build a consensus--however, the entire barcode is discarded if there are many differences between CCSs with the barcode, or high-frequency non-consensus mutations. The reason that barcodes are discarded in such cases as many differences between CCSs or high-frequency non-consensus mutations suggest errors such as barcode collisions or strand exchange.
 
First, call the consensus for each barcode including **all** retained sequences, even those with indels:


```python
consensus, dropped = alignparse.consensus.simple_mutconsensus(
                        processed_ccs.query('retained'),
                        group_cols=('library', 'barcode', 'target'),
                        mutation_col='gene_mutations',
                        )
```

Here are the first few lines of the data frame of consensus sequences for each barcode.
In addition to giving the library, barcode, target, and mutations, it also has a column indicating how many CCSs support the variant call:


```python
display(HTML(consensus.head().to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>barcode</th>
      <th>target</th>
      <th>gene_mutations</th>
      <th>variant_call_support</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>lib40</td>
      <td>AAAAAAAAAAAAATGC</td>
      <td>PRD-0038</td>
      <td>T212G T213G</td>
      <td>1</td>
    </tr>
    <tr>
      <td>lib40</td>
      <td>AAAAAAAAAAAGTCAT</td>
      <td>SARS-CoV-1_2693</td>
      <td>C391G T392C G393T</td>
      <td>1</td>
    </tr>
    <tr>
      <td>lib40</td>
      <td>AAAAAAAAAACATGAG</td>
      <td>SARS-CoV-1_2693</td>
      <td></td>
      <td>2</td>
    </tr>
    <tr>
      <td>lib40</td>
      <td>AAAAAAAAAATCTAGA</td>
      <td>SARS-CoV-1_2693</td>
      <td>A340T del337to339</td>
      <td>2</td>
    </tr>
    <tr>
      <td>lib40</td>
      <td>AAAAAAAAAATTCGGG</td>
      <td>RshSTT182</td>
      <td>A7G C8G A9T</td>
      <td>2</td>
    </tr>
  </tbody>
</table>


Since we retain variants with substitutions but ignore those with indels, add information about substitution mutations and number of indels:


```python
consensus = alignparse.consensus.add_mut_info_cols(
                    consensus,
                    mutation_col='gene_mutations',
                    sub_str_col='substitutions',
                    n_indel_col='number_of_indels',
                    overwrite_cols=True)

display(HTML(consensus.head().to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>barcode</th>
      <th>target</th>
      <th>gene_mutations</th>
      <th>variant_call_support</th>
      <th>substitutions</th>
      <th>number_of_indels</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>lib40</td>
      <td>AAAAAAAAAAAAATGC</td>
      <td>PRD-0038</td>
      <td>T212G T213G</td>
      <td>1</td>
      <td>T212G T213G</td>
      <td>0</td>
    </tr>
    <tr>
      <td>lib40</td>
      <td>AAAAAAAAAAAGTCAT</td>
      <td>SARS-CoV-1_2693</td>
      <td>C391G T392C G393T</td>
      <td>1</td>
      <td>C391G T392C G393T</td>
      <td>0</td>
    </tr>
    <tr>
      <td>lib40</td>
      <td>AAAAAAAAAACATGAG</td>
      <td>SARS-CoV-1_2693</td>
      <td></td>
      <td>2</td>
      <td></td>
      <td>0</td>
    </tr>
    <tr>
      <td>lib40</td>
      <td>AAAAAAAAAATCTAGA</td>
      <td>SARS-CoV-1_2693</td>
      <td>A340T del337to339</td>
      <td>2</td>
      <td>A340T</td>
      <td>1</td>
    </tr>
    <tr>
      <td>lib40</td>
      <td>AAAAAAAAAATTCGGG</td>
      <td>RshSTT182</td>
      <td>A7G C8G A9T</td>
      <td>2</td>
      <td>A7G C8G A9T</td>
      <td>0</td>
    </tr>
  </tbody>
</table>


Plot distribution of number of CCSs supporting each variant call (consensus), indicating whether or not there is an indel:


```python
max_variant_call_support = 6  # group variants with >= this much support

_ = (
 ggplot(consensus
        .assign(variant_call_support=lambda x: numpy.clip(x['variant_call_support'],
                                                          None,
                                                          max_variant_call_support),
                indel_state=lambda x: numpy.where(x['number_of_indels'] > 0,
                                                  'has indel', 'no indel')
                ),
        aes('variant_call_support')) +
 geom_bar() +
 ylab('number of variants') +
 facet_grid('indel_state ~ library') +
 theme(figure_size=(1.75 * nlibs, 3.5),
       panel_grid_major_x=element_blank(),
       ) 
 ).draw()
```


    
![png](process_ccs_files/process_ccs_84_0.png)
    


We see that most variant consensus sequences do **not** have indels, especially if we limit to the more "accurate" ones that have multiple CCSs supporting them.

We will ignore all consensus sequences with indels in the variant-barcode lookup table. 
We do this for two reasons:
 1. When there is just one CCS supporting a consensus, it is less likely to be accurate as indels are the main mode of PacBio error.
 2. For the purposes of our studies, we are interested in point mutations rather than indels anyway.
 
Here are number of valid consensus sequence (no indels) for each library and target:


```python
consensus = consensus.query('number_of_indels == 0')

lib_target_counts = (
    consensus
    .groupby(['library', 'target'])
    .size()
    .rename('consensus sequences')
    .reset_index()
    )

display(HTML(lib_target_counts.to_html(index=False)))

p = (ggplot(lib_target_counts.assign(xlabel=lambda x: x['target'] + ', ' + x['library']),
            aes('xlabel', 'consensus sequences')) +
     geom_point(size=3) +
     theme(figure_size=(0.5 * nlibs * ntargets, 1.75),
           axis_text_x=element_text(angle=90)) +
     xlab('') +
     scale_y_log10()
     )

_ = p.draw()
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>target</th>
      <th>consensus sequences</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>lib40</td>
      <td>PRD-0038</td>
      <td>33714</td>
    </tr>
    <tr>
      <td>lib40</td>
      <td>RmYN02</td>
      <td>24980</td>
    </tr>
    <tr>
      <td>lib40</td>
      <td>RsYN04</td>
      <td>35029</td>
    </tr>
    <tr>
      <td>lib40</td>
      <td>RshSTT182</td>
      <td>36190</td>
    </tr>
    <tr>
      <td>lib40</td>
      <td>SARS-CoV-1_2693</td>
      <td>37000</td>
    </tr>
    <tr>
      <td>lib40</td>
      <td>SARS-CoV-2_2649</td>
      <td>1</td>
    </tr>
    <tr>
      <td>lib41</td>
      <td>PRD-0038</td>
      <td>38282</td>
    </tr>
    <tr>
      <td>lib41</td>
      <td>RmYN02</td>
      <td>31689</td>
    </tr>
    <tr>
      <td>lib41</td>
      <td>RsYN04</td>
      <td>39039</td>
    </tr>
    <tr>
      <td>lib41</td>
      <td>RshSTT182</td>
      <td>37190</td>
    </tr>
    <tr>
      <td>lib41</td>
      <td>SARS-CoV-1_2693</td>
      <td>33829</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>273-2005</td>
      <td>368</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>279-2005</td>
      <td>403</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>AncAsia_MAP</td>
      <td>5572</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>AncAsia_alt</td>
      <td>291</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>AncAsia_tree1</td>
      <td>178</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>AncAsia_tree2</td>
      <td>325</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>AncClade2_MAP</td>
      <td>5013</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>AncClade2_alt</td>
      <td>447</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>AncClade2_alt1_subs-only</td>
      <td>397</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>AncClade2_alt2_del1-only</td>
      <td>385</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>AncClade2_alt3_del2-only</td>
      <td>354</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>AncClade2_alt4_dels-only</td>
      <td>371</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>AncClade2_tree2</td>
      <td>440</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>AncEurAf_alt</td>
      <td>396</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>AncEurAf_tree1</td>
      <td>336</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>AncSARS-CoV-1_MAP</td>
      <td>192</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>AncSARS-CoV-1_alt</td>
      <td>361</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>AncSARS1a_MAP</td>
      <td>4277</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>AncSARS1a_alt</td>
      <td>447</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>AncSARS1a_tree1</td>
      <td>217</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>AncSARS1a_tree2</td>
      <td>394</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>AncSARS1c_MAP</td>
      <td>154</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>AncSARS2a_MAP</td>
      <td>5410</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>AncSARS2a_alt</td>
      <td>210</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>AncSARS2c_MAP</td>
      <td>3543</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>AncSarbecovirus_MAP</td>
      <td>4848</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>AncSarbecovirus_alt</td>
      <td>180</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>AncSarbecovirus_alt1_ins117ins118</td>
      <td>181</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>AncSarbecovirus_tree1</td>
      <td>383</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>As6526</td>
      <td>309</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>BANAL-20-103</td>
      <td>340</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>BANAL-20-236</td>
      <td>323</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>BB9904</td>
      <td>393</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>BM48-31</td>
      <td>4605</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>BtKY72</td>
      <td>4187</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>GD-Pangolin</td>
      <td>4300</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>GX-Pangolin</td>
      <td>254</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>GX2013</td>
      <td>642</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>HKU3-1</td>
      <td>236</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>HKU3-8</td>
      <td>579</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>HeB2013</td>
      <td>427</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>HuB2013</td>
      <td>452</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>JL2012</td>
      <td>544</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>Khosta-1</td>
      <td>280</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>Khosta-2</td>
      <td>410</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>LYRa11</td>
      <td>365</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>Longquan-140</td>
      <td>406</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>PDF-2370</td>
      <td>293</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>PRD-0038</td>
      <td>395</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>RaTG13</td>
      <td>4068</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>Rc-o319</td>
      <td>239</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>Rf1</td>
      <td>321</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>Rf4092</td>
      <td>514</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>RhGB01</td>
      <td>395</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>RmYN02</td>
      <td>427</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>Rp3</td>
      <td>584</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>Rs4081</td>
      <td>411</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>Rs4084</td>
      <td>299</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>Rs4231</td>
      <td>486</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>Rs4237</td>
      <td>430</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>Rs4247</td>
      <td>638</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>Rs7327</td>
      <td>3717</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>RsSHC014</td>
      <td>244</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>RsYN04</td>
      <td>281</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>RshSTT182</td>
      <td>259</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>SARS-CoV-1_2693</td>
      <td>5345</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>SARS-CoV-1_GD01_HP03L</td>
      <td>224</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>SARS-CoV-1_GD03T0013_HP04</td>
      <td>438</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>SARS-CoV-1_GZ-C_HP03L</td>
      <td>160</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>SARS-CoV-1_GZ0402_HP04</td>
      <td>255</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>SARS-CoV-1_HGZ8L1-A_HP03E</td>
      <td>310</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>SARS-CoV-1_PC4-127_PC04</td>
      <td>442</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>SARS-CoV-1_PC4-137_PC04</td>
      <td>3174</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>SARS-CoV-1_PC4-13_PC04</td>
      <td>404</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>SARS-CoV-1_SZ1_PC03</td>
      <td>331</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>SARS-CoV-1_Sin852_HP03L</td>
      <td>353</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>SARS-CoV-1_Sino1-11_HP03L</td>
      <td>451</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>SARS-CoV-1_Urbani_HP03L</td>
      <td>4053</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>SARS-CoV-2_2649</td>
      <td>803</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>SARS-CoV-2_Alpha</td>
      <td>463</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>SARS-CoV-2_Beta</td>
      <td>357</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>SARS-CoV-2_Delta</td>
      <td>362</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>SARS-CoV-2_WH1</td>
      <td>2403</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>Shaanxi2011</td>
      <td>462</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>WIV1</td>
      <td>376</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>YN2013</td>
      <td>378</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>YN2020B</td>
      <td>387</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>Yunnan2011</td>
      <td>518</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>ZC45</td>
      <td>460</td>
    </tr>
    <tr>
      <td>lib46</td>
      <td>ZXC21</td>
      <td>765</td>
    </tr>
  </tbody>
</table>



    
![png](process_ccs_files/process_ccs_86_1.png)
    


Below we write the retained consensus sequences to a CSV file that links the nucleotide mutations to the barcodes:


```python
print(f"Writing nucleotide variants to {config['nt_variant_table_file']}")
      
(consensus
 [['target', 'library', 'barcode', 'substitutions', 'variant_call_support','number_of_indels']]
 .to_csv(config['nt_variant_table_file'], index=False)
 )
      
print('Here are the first few lines of this file:')
display(HTML(
    pd.read_csv(config['nt_variant_table_file'], na_filter=None)
    .head()
    .to_html(index=False)
    ))
```

    Writing nucleotide variants to results/variants/nucleotide_variant_table.csv
    Here are the first few lines of this file:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>target</th>
      <th>library</th>
      <th>barcode</th>
      <th>substitutions</th>
      <th>variant_call_support</th>
      <th>number_of_indels</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>PRD-0038</td>
      <td>lib40</td>
      <td>AAAAAAAAAAAAATGC</td>
      <td>T212G T213G</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <td>SARS-CoV-1_2693</td>
      <td>lib40</td>
      <td>AAAAAAAAAAAGTCAT</td>
      <td>C391G T392C G393T</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <td>SARS-CoV-1_2693</td>
      <td>lib40</td>
      <td>AAAAAAAAAACATGAG</td>
      <td></td>
      <td>2</td>
      <td>0</td>
    </tr>
    <tr>
      <td>RshSTT182</td>
      <td>lib40</td>
      <td>AAAAAAAAAATTCGGG</td>
      <td>A7G C8G A9T</td>
      <td>2</td>
      <td>0</td>
    </tr>
    <tr>
      <td>RshSTT182</td>
      <td>lib40</td>
      <td>AAAAAAAAACAAACTA</td>
      <td>G151A G152A C153A</td>
      <td>1</td>
      <td>0</td>
    </tr>
  </tbody>
</table>


What happened to the barcodes that we "dropped" because we could not construct a reliable consensus?
The `dropped` data frame from [alignparse.consensus.simple_mutconsensus](https://jbloomlab.github.io/alignparse/alignparse.consensus.html?highlight=simple_mutconsensus#alignparse.consensus.simple_mutconsensus) has this information:


```python
display(HTML(dropped.head().to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>barcode</th>
      <th>target</th>
      <th>drop_reason</th>
      <th>nseqs</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>lib40</td>
      <td>AAAAAAAATATCAGTG</td>
      <td>PRD-0038</td>
      <td>subs diff too large</td>
      <td>2</td>
    </tr>
    <tr>
      <td>lib40</td>
      <td>AAAAAAAATATTATAA</td>
      <td>PRD-0038</td>
      <td>subs diff too large</td>
      <td>4</td>
    </tr>
    <tr>
      <td>lib40</td>
      <td>AAAAAAACACCTGACA</td>
      <td>RshSTT182</td>
      <td>minor subs too frequent</td>
      <td>7</td>
    </tr>
    <tr>
      <td>lib40</td>
      <td>AAAAAAACTTGCTGAC</td>
      <td>RsYN04</td>
      <td>subs diff too large</td>
      <td>7</td>
    </tr>
    <tr>
      <td>lib40</td>
      <td>AAAAAAAGGTTTTCAA</td>
      <td>SARS-CoV-1_2693</td>
      <td>subs diff too large</td>
      <td>4</td>
    </tr>
  </tbody>
</table>


Summarize the information in this data frame on dropped barcodes with the plot below.
This plot shows several things.
First, we see that the total number of barcodes dropped is modest (just a few thousand per library) relative to the total number of barcodes per library (seen above to be on the order of hundreds of thousands).
Second, the main reason that barcodes are dropped is that there are CCSs within the same barcode with suspiciously large numbers of mutations relative to the consensus---which we use as a filter to discard the entire barcode as it could indicate strand exchange or some other issue.
In any case, the modest number of dropped barcodes indicates that there probably isn't much of a need to worry: 


```python
max_nseqs = 8  # plot together all barcodes with >= this many sequences

_ = (
 ggplot(
    dropped.assign(nseqs=lambda x: numpy.clip(x['nseqs'], None, max_nseqs)),
    aes('nseqs')) + 
 geom_bar() + 
 scale_x_continuous(limits=(1, None)) +
 xlab('number of sequences for barcode') +
 ylab('number of barcodes') +
 facet_grid('library ~ drop_reason') +
 theme(figure_size=(10, 1.5 * nlibs),
       panel_grid_major_x=element_blank(),
       )
 ).draw()
```


    
![png](process_ccs_files/process_ccs_92_0.png)
    



```python

```
