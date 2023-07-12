"""``snakemake`` file that runs entire analysis."""

# Imports ---------------------------------------------------------------------
import glob
import itertools
import os.path
import os
import textwrap
import urllib.request

import Bio.SeqIO

import dms_variants.codonvarianttable
import dms_variants.illuminabarcodeparser

import pandas as pd

# Configuration  --------------------------------------------------------------
configfile: 'config.yaml'

# run "quick" rules locally:
localrules: make_dag,
            make_summary

# Functions -------------------------------------------------------------------
def nb_markdown(nb):
    """Return path to Markdown results of notebook `nb`."""
    return os.path.join(config['summary_dir'],
                        os.path.basename(os.path.splitext(nb)[0]) + '.md')

# Global variables extracted from config --------------------------------------
pacbio_runs = (pd.read_csv(config['pacbio_runs'], dtype = str)
               .assign(pacbioRun=lambda x: x['library'] + '_' + x['run'])
               )
assert len(pacbio_runs['pacbioRun'].unique()) == len(pacbio_runs['pacbioRun'])

# Information on samples and barcode runs -------------------------------------
barcode_runs = pd.read_csv(config['barcode_runs'])

# Rules -----------------------------------------------------------------------

# making this summary is the target rule (in place of `all`) since it
# is first rule listed.
rule make_summary:
    """Create Markdown summary of analysis."""
    input:
        dag=os.path.join(config['summary_dir'], 'dag.svg'),
        WH1_get_bc_variant_lookup=config['WH1_bc_variant_lookup'],
        WH1_get_mut_antibody_escape=config['WH1_mut_antibody_escape'],
        process_ccs=nb_markdown('process_ccs.ipynb'),
        build_variants_RshSTT182=nb_markdown('build_variants_RshSTT182.ipynb'),
        build_variants_PRD0038=nb_markdown('build_variants_PRD0038.ipynb'),
        build_variants_SARS1=nb_markdown('build_variants_SARS1.ipynb'),
        build_variants_RsYN04=nb_markdown('build_variants_RsYN04.ipynb'),
        build_variants_RmYN02=nb_markdown('build_variants_RmYN02.ipynb'),
        build_variants_lib46='results/summary/build_variants_lib46.md',
        barcode_variant_table_RshSTT182=config['codon_variant_table_file_RshSTT182'],
        barcode_variant_table_PRD0038=config['codon_variant_table_file_PRD-0038'],
        barcode_variant_table_SARS1=config['codon_variant_table_file_SARS-CoV-1_2693'],
        barcode_variant_table_RsYN04=config['codon_variant_table_file_RsYN04'],
        barcode_variant_table_RmYN02=config['codon_variant_table_file_RmYN02'], 
        barcode_variant_table_lib46=config['codon_variant_table_file_lib46'], 
        nt_variant_table=config['nt_variant_table_file'],
        variant_counts_file=config['variant_counts_file'],
        count_variants=nb_markdown('count_variants.ipynb'),
        fit_titrations_huACE2='results/summary/compute_binding_Kd_huACE2.md',
        variant_Kds_file_huACE2=config['Titeseq_Kds_file_huACE2'],
        fit_titrations_RshACE2='results/summary/compute_binding_Kd_RshACE2.md',
        variant_Kds_file_RshACE2=config['Titeseq_Kds_file_RshACE2'],
        fit_titrations_Ra787ACE2='results/summary/compute_binding_Kd_Ra787ACE2.md',
        variant_Kds_file_Ra787ACE2=config['Titeseq_Kds_file_Ra787ACE2'],
        fit_titrations_Ra9479ACE2='results/summary/compute_binding_Kd_Ra9479ACE2.md',
        variant_Kds_file_Ra9479ACE2=config['Titeseq_Kds_file_Ra9479ACE2'],
        fit_titrations_RlanACE2='results/summary/compute_binding_Kd_RlanACE2.md',
        variant_Kds_file_RlanACE2=config['Titeseq_Kds_file_RlanACE2'],
        fit_titrations_RalcACE2='results/summary/compute_binding_Kd_RalcACE2.md',
        variant_Kds_file_RalcACE2=config['Titeseq_Kds_file_RalcACE2'],
        fit_titrations_RpearACE2='results/summary/compute_binding_Kd_RpearACE2.md',
        variant_Kds_file_RpearACE2=config['Titeseq_Kds_file_RpearACE2'],
        calculate_expression='results/summary/compute_expression_meanF.md',
        variant_expression_file=config['expression_sortseq_file'],
        collapse_bc_lib40_41='results/summary/collapse_barcodes_lib40_41.md',
        collapse_bc_lib40_41_file=config['final_variant_scores_lib40_41_file'],
        collapse_bc_lib46='results/summary/collapse_barcodes_lib46.md',
        collapse_bc_lib46_muts_file=config['final_variant_scores_lib46_muts_file'],
        collapse_bc_lib46_wts_file=config['final_variant_scores_lib46_wts_file'],
        epistatic_shifts='results/summary/epistatic_shifts.md',
        heatmap_viz_delta=os.path.join(config['visualization_dir'], "heatmaps_delta.html"),
        heatmap_viz_absolute=os.path.join(config['visualization_dir'], "heatmaps_absolute.html"),
    output:
        summary = os.path.join(config['summary_dir'], 'summary.md')
    run:
        def path(f):
            """Get path relative to `summary_dir`."""
            return os.path.relpath(f, config['summary_dir'])
        with open(output.summary, 'w') as f:
            f.write(textwrap.dedent(f"""
            # Summary

            Analysis run by [Snakefile]({path(workflow.snakefile)})
            using [this config file]({path(workflow.configfiles[0])}).
            See the [README in the top directory]({path('README.md')})
            for details.

            Here is the DAG of the computational workflow:
            ![{path(input.dag)}]({path(input.dag)})

            Here is the Markdown output of each Jupyter notebook in the
            workflow:



            1. Get prior Wuhan-1 RBD barcode-variant lookup table from the [SARS-CoV-2-RBD_DMS_variants repository](https://github.com/jbloomlab/SARS-CoV-2-RBD_DMS_variants) and antibody escape data from the [escape map aggregator repository](https://github.com/jbloomlab/SARS2_RBD_Ab_escape_maps/blob/main/processed_data/escape_data.csv). 
            
            2. [Process PacBio CCSs]({path(input.process_ccs)}) for each SARSr background and generate barcode consensus sequences.
            
            3. Build barcode-variant lookup tables for the [RshSTT182]({path(input.build_variants_RshSTT182)}), [PRD-0038]({path(input.build_variants_PRD0038)}), [SARS-CoV-1 Urbani]({path(input.build_variants_SARS1)}), [RsYN04]({path(input.build_variants_RsYN04)}), and [RmYN02]({path(input.build_variants_RmYN02)}) backgrounds, as well as the [lib46 pan-sarbecoviruses pool]({path(input.build_variants_lib46)}). Barcode-variant lookup tables are saved for each background: [RshSTT182]({path(input.barcode_variant_table_RshSTT182)}), [PRD-0038]({path(input.barcode_variant_table_PRD0038)}), [SARS-CoV-1]({path(input.barcode_variant_table_SARS1)}), [RsYN04]({path(input.barcode_variant_table_RsYN04)}), and [RmYN02]({path(input.barcode_variant_table_RmYN02)}), as well as the [lib46 pan-sarbecovirus pool]({path(input.barcode_variant_table_lib46)}).

            4. [Count variants by barcode]({path(input.count_variants)}).
               Creates a [variant counts file]({path(input.variant_counts_file)})
               giving counts of each barcoded variant in each condition.

            5. Fit titration curves for binding to 
            [human ACE2]({path(input.fit_titrations_huACE2)}), 
            [R. shameli ACE2]({path(input.fit_titrations_RshACE2)}), 
            [R. affinis 787 ACE2]({path(input.fit_titrations_Ra787ACE2)}), 
            [R. affinis 9479 ACE2]({path(input.fit_titrations_Ra9479ACE2)}), 
            [R. landeri ACE2]({path(input.fit_titrations_RlanACE2)}), 
            [R. alcyone ACE2]({path(input.fit_titrations_RalcACE2)}), and 
            [R. pearsonii ACE2]({path(input.fit_titrations_RpearACE2)}) (RpearACE2 is MFI not Kd), to calculate per-barcode K<sub>D</sub>, recorded in these files for 
            [huACE2]({path(input.variant_Kds_file_huACE2)}), 
            [RshACE2]({path(input.variant_Kds_file_RshACE2)}), 
            [Ra787ACE2]({path(input.variant_Kds_file_Ra787ACE2)}), 
            [Ra9479ACE2]({path(input.variant_Kds_file_Ra9479ACE2)}), 
            [RlanACE2]({path(input.variant_Kds_file_RlanACE2)}), 
            [RalcACE2]({path(input.variant_Kds_file_RalcACE2)}), and
            [RpearACE2]({path(input.variant_Kds_file_RpearACE2)}) (RpearACE2 is MFI not Kd).
          
            6. [Analyze Sort-seq]({path(input.calculate_expression)}) to calculate per-barcode RBD expression, recorded in [this file]({path(input.variant_expression_file)}).
            
            7. Collapse internal replicate barcodes of each variant to final variant phenotypes for the pan-sarbecovirus lib46 pool. Analysis [here]({path(input.collapse_bc_lib46)}) and final output file [here for SSM mutants]({path(input.collapse_bc_lib46_muts_file)}) and [here for wildtype variants]({path(input.collapse_bc_lib46_wts_file)}).
            
            8. Collapse internal replicate barcodes of each variant to final variant phenotypes for the sarbecovirus DMS pools. Analysis [here]({path(input.collapse_bc_lib40_41)}) and final output file [here]({path(input.collapse_bc_lib40_41_file)}).
            
            9. [Analyze patterns of epistasis in the DMS data]({path(input.epistatic_shifts)}).
            
            10. Make interactive data visualizations, available [here](https://jbloomlab.github.io/SARSr-CoV-RBD_DMS/)
                        
            """
            ).strip())

rule make_dag:
    # error message, but works: https://github.com/sequana/sequana/issues/115
    input:
        workflow.snakefile
    output:
        os.path.join(config['summary_dir'], 'dag.svg')
    shell:
        "snakemake --forceall --dag | dot -Tsvg > {output}"

rule interactive_heatmap_absolute:
    """ Make the interactive heatmaps for absolute expression and binding phenotypes.
    """
    input: 
        scores=config['final_variant_scores_lib40_41_file']
    params:
        annotations=config['RBD_sites']
    output:
        html=os.path.join(config['visualization_dir'], "heatmaps_absolute.html")
    notebook: "RBD-Heatmaps-Interactive-Visualization_absolute.ipynb"

rule interactive_heatmap_delta:
    """ Make the interactive heatmaps for delta expression and binding phenotypes.
    """
    input: 
        scores=config['final_variant_scores_lib40_41_file']
    params:
        annotations=config['RBD_sites']
    output:
        html=os.path.join(config['visualization_dir'], "heatmaps_delta.html")
    notebook: "RBD-Heatmaps-Interactive-Visualization_delta.ipynb"

rule epistatic_shifts:
    input:
        config['final_variant_scores_lib40_41_file'],
        config['WH1_mut_antibody_escape'],
    output:
        config['JSD_huACE2_file'],
        config['JSD_Ra787ACE2_file'],
        config['JSD_Ra9479ACE2_file'],
        config['JSD_RlanACE2_file'],
        config['JSD_RalcACE2_file'],
        config['JSD_RshACE2_file'],
        config['JSD_expr_file'],
        md='results/summary/epistatic_shifts.md',
        md_files=directory('results/summary/epistatic_shifts_files')
    envmodules:
        'R/3.6.2-foss-2019b'
    params:
        nb='epistatic_shifts.Rmd',
        md='epistatic_shifts.md',
        md_files='epistatic_shifts_files'
    shell:
        """
        R -e \"rmarkdown::render(input=\'{params.nb}\')\";
        mv {params.md} {output.md};
        mv {params.md_files} {output.md_files}
        """

rule collapse_bcs_lib40_41_SARSr_DMS:
    input:
        config['Titeseq_Kds_file_huACE2'],
        config['Titeseq_Kds_file_RshACE2'],
        config['Titeseq_Kds_file_Ra787ACE2'],
        config['Titeseq_Kds_file_Ra9479ACE2'],
        config['Titeseq_Kds_file_RlanACE2'],
        config['Titeseq_Kds_file_RalcACE2'],
        config['Titeseq_Kds_file_RpearACE2'],
        config['expression_sortseq_file'],
        config['RBD_annotation_file'],
        config['final_variant_scores_lib46_muts_file'],
    output:
        config['final_variant_scores_lib40_41_file'],
        md='results/summary/collapse_barcodes_lib40_41.md',
        md_files=directory('results/summary/collapse_barcodes_lib40_41_files')
    envmodules:
        'R/3.6.2-foss-2019b'
    params:
        nb='collapse_barcodes_lib40_41.Rmd',
        md='collapse_barcodes_lib40_41.md',
        md_files='collapse_barcodes_lib40_41_files'
    shell:
        """
        R -e \"rmarkdown::render(input=\'{params.nb}\')\";
        mv {params.md} {output.md};
        mv {params.md_files} {output.md_files}
        """

rule collapse_bcs_lib46_SARSr_wts:
    input:
        config['Titeseq_Kds_file_huACE2'],
        config['Titeseq_Kds_file_RshACE2'],
        config['Titeseq_Kds_file_Ra787ACE2'],
        config['Titeseq_Kds_file_Ra9479ACE2'],
        config['Titeseq_Kds_file_RlanACE2'],
        config['Titeseq_Kds_file_RalcACE2'],
        config['Titeseq_Kds_file_RpearACE2'],
        config['SARSr_lib46_mut_bind_expr'],
        config['SARSr_lib47_mut_bind_expr'],
    output:
        config['final_variant_scores_lib46_muts_file'],
        config['final_variant_scores_lib46_wts_file'],
        md='results/summary/collapse_barcodes_lib46.md',
        md_files=directory('results/summary/collapse_barcodes_lib46_files')
    envmodules:
        'R/3.6.2-foss-2019b'
    params:
        nb='collapse_barcodes_lib46.Rmd',
        md='collapse_barcodes_lib46.md',
        md_files='collapse_barcodes_lib46_files'
    shell:
        """
        R -e \"rmarkdown::render(input=\'{params.nb}\')\";
        mv {params.md} {output.md};
        mv {params.md_files} {output.md_files}
        """

rule fit_titrations_huACE2:
    input:
        config['codon_variant_table_file_RshSTT182'],
        config['codon_variant_table_file_PRD-0038'],
        config['codon_variant_table_file_SARS-CoV-1_2693'],
        config['codon_variant_table_file_RsYN04'],
        config['codon_variant_table_file_RmYN02'],
        config['WH1_bc_variant_lookup'],
        config['codon_variant_table_file_lib46'],
        config['variant_counts_file']
    output:
        config['Titeseq_Kds_file_huACE2'],
        md='results/summary/compute_binding_Kd_huACE2.md',
        md_files=directory('results/summary/compute_binding_Kd_huACE2_files')
    envmodules:
        'R/3.6.2-foss-2019b'
    params:
        nb='compute_binding_Kd_huACE2.Rmd',
        md='compute_binding_Kd_huACE2.md',
        md_files='compute_binding_Kd_huACE2_files'
    shell:
        """
        R -e \"rmarkdown::render(input=\'{params.nb}\')\";
        mv {params.md} {output.md};
        mv {params.md_files} {output.md_files}
        """

rule fit_titrations_RshACE2:
    input:
        config['codon_variant_table_file_RshSTT182'],
        config['codon_variant_table_file_PRD-0038'],
        config['codon_variant_table_file_SARS-CoV-1_2693'],
        config['codon_variant_table_file_RsYN04'],
        config['codon_variant_table_file_RmYN02'],
        config['WH1_bc_variant_lookup'],
        config['codon_variant_table_file_lib46'],
        config['variant_counts_file']
    output:
        config['Titeseq_Kds_file_RshACE2'],
        md='results/summary/compute_binding_Kd_RshACE2.md',
        md_files=directory('results/summary/compute_binding_Kd_RshACE2_files')
    envmodules:
        'R/3.6.2-foss-2019b'
    params:
        nb='compute_binding_Kd_RshACE2.Rmd',
        md='compute_binding_Kd_RshACE2.md',
        md_files='compute_binding_Kd_RshACE2_files'
    shell:
        """
        R -e \"rmarkdown::render(input=\'{params.nb}\')\";
        mv {params.md} {output.md};
        mv {params.md_files} {output.md_files}
        """

rule fit_titrations_Ra787ACE2:
    input:
        config['codon_variant_table_file_RshSTT182'],
        config['codon_variant_table_file_PRD-0038'],
        config['codon_variant_table_file_SARS-CoV-1_2693'],
        config['codon_variant_table_file_RsYN04'],
        config['codon_variant_table_file_RmYN02'],
        config['WH1_bc_variant_lookup'],
        config['codon_variant_table_file_lib46'],
        config['variant_counts_file']
    output:
        config['Titeseq_Kds_file_Ra787ACE2'],
        md='results/summary/compute_binding_Kd_Ra787ACE2.md',
        md_files=directory('results/summary/compute_binding_Kd_Ra787ACE2_files')
    envmodules:
        'R/3.6.2-foss-2019b'
    params:
        nb='compute_binding_Kd_Ra787ACE2.Rmd',
        md='compute_binding_Kd_Ra787ACE2.md',
        md_files='compute_binding_Kd_Ra787ACE2_files'
    shell:
        """
        R -e \"rmarkdown::render(input=\'{params.nb}\')\";
        mv {params.md} {output.md};
        mv {params.md_files} {output.md_files}
        """
        
rule fit_titrations_Ra9479ACE2:
    input:
        config['codon_variant_table_file_RshSTT182'],
        config['codon_variant_table_file_PRD-0038'],
        config['codon_variant_table_file_SARS-CoV-1_2693'],
        config['codon_variant_table_file_RsYN04'],
        config['codon_variant_table_file_RmYN02'],
        config['WH1_bc_variant_lookup'],
        config['codon_variant_table_file_lib46'],
        config['variant_counts_file']
    output:
        config['Titeseq_Kds_file_Ra9479ACE2'],
        md='results/summary/compute_binding_Kd_Ra9479ACE2.md',
        md_files=directory('results/summary/compute_binding_Kd_Ra9479ACE2_files')
    envmodules:
        'R/3.6.2-foss-2019b'
    params:
        nb='compute_binding_Kd_Ra9479ACE2.Rmd',
        md='compute_binding_Kd_Ra9479ACE2.md',
        md_files='compute_binding_Kd_Ra9479ACE2_files'
    shell:
        """
        R -e \"rmarkdown::render(input=\'{params.nb}\')\";
        mv {params.md} {output.md};
        mv {params.md_files} {output.md_files}
        """
        
rule fit_titrations_RlanACE2:
    input:
        config['codon_variant_table_file_RshSTT182'],
        config['codon_variant_table_file_PRD-0038'],
        config['codon_variant_table_file_SARS-CoV-1_2693'],
        config['codon_variant_table_file_RsYN04'],
        config['codon_variant_table_file_RmYN02'],
        config['WH1_bc_variant_lookup'],
        config['codon_variant_table_file_lib46'],
        config['variant_counts_file']
    output:
        config['Titeseq_Kds_file_RlanACE2'],
        md='results/summary/compute_binding_Kd_RlanACE2.md',
        md_files=directory('results/summary/compute_binding_Kd_RlanACE2_files')
    envmodules:
        'R/3.6.2-foss-2019b'
    params:
        nb='compute_binding_Kd_RlanACE2.Rmd',
        md='compute_binding_Kd_RlanACE2.md',
        md_files='compute_binding_Kd_RlanACE2_files'
    shell:
        """
        R -e \"rmarkdown::render(input=\'{params.nb}\')\";
        mv {params.md} {output.md};
        mv {params.md_files} {output.md_files}
        """
        
rule fit_titrations_RalcACE2:
    input:
        config['codon_variant_table_file_RshSTT182'],
        config['codon_variant_table_file_PRD-0038'],
        config['codon_variant_table_file_SARS-CoV-1_2693'],
        config['codon_variant_table_file_RsYN04'],
        config['codon_variant_table_file_RmYN02'],
        config['WH1_bc_variant_lookup'],
        config['codon_variant_table_file_lib46'],
        config['variant_counts_file']
    output:
        config['Titeseq_Kds_file_RalcACE2'],
        md='results/summary/compute_binding_Kd_RalcACE2.md',
        md_files=directory('results/summary/compute_binding_Kd_RalcACE2_files')
    envmodules:
        'R/3.6.2-foss-2019b'
    params:
        nb='compute_binding_Kd_RalcACE2.Rmd',
        md='compute_binding_Kd_RalcACE2.md',
        md_files='compute_binding_Kd_RalcACE2_files'
    shell:
        """
        R -e \"rmarkdown::render(input=\'{params.nb}\')\";
        mv {params.md} {output.md};
        mv {params.md_files} {output.md_files}
        """
        
rule fit_titrations_RpearACE2:
    input:
        config['codon_variant_table_file_RshSTT182'],
        config['codon_variant_table_file_PRD-0038'],
        config['codon_variant_table_file_SARS-CoV-1_2693'],
        config['codon_variant_table_file_RsYN04'],
        config['codon_variant_table_file_RmYN02'],
        config['WH1_bc_variant_lookup'],
        config['codon_variant_table_file_lib46'],
        config['variant_counts_file']
    output:
        config['Titeseq_Kds_file_RpearACE2'],
        md='results/summary/compute_binding_Kd_RpearACE2.md',
        md_files=directory('results/summary/compute_binding_Kd_RpearACE2_files')
    envmodules:
        'R/3.6.2-foss-2019b'
    params:
        nb='compute_binding_Kd_RpearACE2.Rmd',
        md='compute_binding_Kd_RpearACE2.md',
        md_files='compute_binding_Kd_RpearACE2_files'
    shell:
        """
        R -e \"rmarkdown::render(input=\'{params.nb}\')\";
        mv {params.md} {output.md};
        mv {params.md_files} {output.md_files}
        """

rule calculate_expression:
    input:
        config['codon_variant_table_file_RshSTT182'],
        config['codon_variant_table_file_PRD-0038'],
        config['codon_variant_table_file_SARS-CoV-1_2693'],
        config['codon_variant_table_file_RsYN04'],
        config['codon_variant_table_file_RmYN02'],
        config['WH1_bc_variant_lookup'],
        config['codon_variant_table_file_lib46'],
        config['variant_counts_file']
    output:
        config['expression_sortseq_file'],
        md='results/summary/compute_expression_meanF.md',
        md_files=directory('results/summary/compute_expression_meanF_files')
    envmodules:
        'R/3.6.2-foss-2019b'
    params:
        nb='compute_expression_meanF.Rmd',
        md='compute_expression_meanF.md',
        md_files='compute_expression_meanF_files'
    shell:
        """
        R -e \"rmarkdown::render(input=\'{params.nb}\')\";
        mv {params.md} {output.md};
        mv {params.md_files} {output.md_files}
        """


rule count_variants:
    """Count codon variants from Illumina barcode runs."""
    input:
        config['barcode_runs'],
        config['codon_variant_table_file_RshSTT182'],
        config['codon_variant_table_file_PRD-0038'],
        config['codon_variant_table_file_SARS-CoV-1_2693'],
        config['codon_variant_table_file_RsYN04'],
        config['codon_variant_table_file_RmYN02'],
        config['WH1_bc_variant_lookup'],
        config['codon_variant_table_file_lib46']
    output:
        config['variant_counts_file'],
        nb_markdown=nb_markdown('count_variants.ipynb')
    params:
        nb='count_variants.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"


rule build_variants_lib46_wildtypes:
    input:
        config['nt_variant_table_file'],
    output:
        config['codon_variant_table_file_lib46'],
        md='results/summary/build_variants_lib46.md',
        md_files=directory('results/summary/build_variants_lib46_files')
    envmodules:
        'R/3.6.2-foss-2019b'
    params:
        nb='build_variants_lib46.Rmd',
        md='build_variants_lib46.md',
        md_files='build_variants_lib46_files'
    shell:
        """
        R -e \"rmarkdown::render(input=\'{params.nb}\')\";
        mv {params.md} {output.md};
        mv {params.md_files} {output.md_files}
        """

rule build_variants_RmYN02:
    """Build barcode-variant lookup table for RmYN02 background."""
    input:
        config['nt_variant_table_file'],
    output:
        config['codon_variant_table_file_RmYN02'],
        nb_markdown=nb_markdown('build_variants_RmYN02.ipynb')
    params:
        nb='build_variants_RmYN02.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

rule build_variants_RsYN04:
    """Build barcode-variant lookup table for RsYN04 background."""
    input:
        config['nt_variant_table_file'],
    output:
        config['codon_variant_table_file_RsYN04'],
        nb_markdown=nb_markdown('build_variants_RsYN04.ipynb')
    params:
        nb='build_variants_RsYN04.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

rule build_variants_SARS1:
    """Build barcode-variant lookup table for SARS1 background."""
    input:
        config['nt_variant_table_file'],
    output:
        config['codon_variant_table_file_SARS-CoV-1_2693'],
        nb_markdown=nb_markdown('build_variants_SARS1.ipynb')
    params:
        nb='build_variants_SARS1.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

rule build_variants_PRD0038:
    """Build barcode-variant lookup table for PRD0038 background."""
    input:
        config['nt_variant_table_file'],
    output:
        config['codon_variant_table_file_PRD-0038'],
        nb_markdown=nb_markdown('build_variants_PRD0038.ipynb')
    params:
        nb='build_variants_PRD0038.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

rule build_variants_RshSTT182:
    """Build barcode-variant lookup table for RshSTT182 background."""
    input:
        config['nt_variant_table_file'],
    output:
        config['codon_variant_table_file_RshSTT182'],
        nb_markdown=nb_markdown('build_variants_RshSTT182.ipynb')
    params:
        nb='build_variants_RshSTT182.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

rule process_ccs:
    """Process the PacBio CCSs."""
    input:
        expand(os.path.join(config['ccs_dir'], "{pacbioRun}_ccs.fastq.gz"),
               pacbioRun=pacbio_runs['pacbioRun']),
        config['amplicons'],
    output:
        config['processed_ccs_file'],
        config['nt_variant_table_file'],
        nb_markdown=nb_markdown('process_ccs.ipynb')
    params:
        nb='process_ccs.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

rule get_SARSr_mut_bind_expr:
    """Download SARSr mini-SSM mutant library ACE2-binding and expression scores from URL."""
    output:
        file=config['SARSr_lib46_mut_bind_expr']
    run:
        urllib.request.urlretrieve(config['SARSr_lib46_mut_bind_expr_url'], output.file)

rule get_SARSr_wts_bind_expr:
    """Download SARSr wts library ACE2-binding and expression scores from URL."""
    output:
        file=config['SARSr_lib47_mut_bind_expr']
    run:
        urllib.request.urlretrieve(config['SARSr_lib47_mut_bind_expr_url'], output.file)

rule get_WH1_bc_variant_lookup:
    """Download SARS-CoV-2 Wuhan-Hu-1 bc-variant lookup table from URL."""
    output:
        file=config['WH1_bc_variant_lookup']
    run:
        urllib.request.urlretrieve(config['WH1_bc_variant_lookup_url'], output.file)

rule get_WH1_mut_antibody_escape:
    """Download SARS-CoV-2 mutation antibody-escape data from URL."""
    output:
        file=config['WH1_mut_antibody_escape']
    run:
        urllib.request.urlretrieve(config['WH1_mut_antibody_escape_url'], output.file)
        
rule build_ccs:
    """Run PacBio ``ccs`` program to build CCSs from subreads or get ccs for samples already processed on-instrument."""
    input:
        input_file=lambda wildcards: (pacbio_runs
                                    .set_index('pacbioRun')
                                    .at[wildcards.pacbioRun, 'input_file']
                                    )
    output:
        ccs_fastq=os.path.join(config['ccs_dir'], "{pacbioRun}_ccs.fastq.gz")
    params:
        min_ccs_length=config['min_ccs_length'],
        max_ccs_length=config['max_ccs_length'],
        min_ccs_passes=config['min_ccs_passes'],
        min_ccs_accuracy=config['min_ccs_accuracy'],
        fastq_prebuilt=lambda wc: (pacbio_runs
                                    .set_index('pacbioRun')
                                    .at[wc.pacbioRun, 'fastq_prebuilt']
                                    )
    threads: config['max_cpus']
    shell:
        """
        if [ {params.fastq_prebuilt} == "yes" ]; then
        	ln -s {input.input_file} {output.ccs_fastq}
        else
        	ccs \
            --min-length {params.min_ccs_length} \
            --max-length {params.max_ccs_length} \
            --min-passes {params.min_ccs_passes} \
            --min-rq {params.min_ccs_accuracy} \
            --num-threads {threads} \
            {input.input_file} \
            {output.ccs_fastq}
        fi
        
        """

