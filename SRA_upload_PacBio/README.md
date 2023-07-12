# Uploading FASTQ files to the SRA

Details of how the raw PacBio sequencing files were uploaded to the NIH [Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra).
The submitted files are in BioProject [PRJNA770094](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA770094).

We want to merge multiple PacBio read sets into a single `fastq.gz` file, since multiples files with different directory paths but the same file name (replicate PacBio lanes) can cause issues for the upload.

To create our merged SARSr-DMS PacBio read sets, we ran:
```
cat ../results/ccs/lib40_210917_ccs.fastq.gz /fh/fast/bloom_j/SR/ngs/pacbio/211118_TylerStarr/211108_pool_Starr-Cell1/demultiplex.bc1003_BAK8A_OA--bc1003_BAK8A_OA.hifi_reads.fastq.gz /fh/fast/bloom_j/SR/ngs/pacbio/211118_TylerStarr/211108_pool_Starr-Cell2/demultiplex.bc1003_BAK8A_OA--bc1003_BAK8A_OA.hifi_reads.fastq.gz /fh/fast/bloom_j/SR/ngs/pacbio/211118_TylerStarr/211108_pool_Starr-Cell1/demultiplex.bc1008_BAK8A_OA--bc1008_BAK8A_OA.hifi_reads.fastq.gz /fh/fast/bloom_j/SR/ngs/pacbio/211118_TylerStarr/211108_pool_Starr-Cell2/demultiplex.bc1008_BAK8A_OA--bc1008_BAK8A_OA.hifi_reads.fastq.gz > ./lib40_ccs.fastq.gz 

cat ../results/ccs/lib41_210917_ccs.fastq.gz > ./lib41_ccs.fastq.gz

cat /fh/fast/bloom_j/computational_notebooks/tstarr/2020/SARSr-CoV_homolog_survey/results/ccs/lib1_200916_1_ccs.fastq.gz /fh/fast/bloom_j/computational_notebooks/tstarr/2020/SARSr-CoV_homolog_survey/results/ccs/lib1_200916_2_ccs.fastq.gz /fh/fast/bloom_j/computational_notebooks/tstarr/2020/SARSr-CoV_homolog_survey/results/ccs/lib1_201110_ccs.fastq.gz /fh/fast/bloom_j/SR/ngs/pacbio/211118_TylerStarr/211108_pool_Starr-Cell1/demultiplex.bc1009_BAK8A_OA--bc1009_BAK8A_OA.hifi_reads.fastq.gz /fh/fast/bloom_j/SR/ngs/pacbio/211118_TylerStarr/211108_pool_Starr-Cell2/demultiplex.bc1009_BAK8A_OA--bc1009_BAK8A_OA.hifi_reads.fastq.gz > lib46_ccs.fastq.gz
```

Remember to delete these fastq.gz files along with the .tar archives after upload is complete.

The Python Jupyter notebooks [upload_to_SRA_SARSr-DMS.ipynb](upload_to_SRA_SARSr-DMS.ipynb) and [upload_to_SRA_pan-sarbeco.ipynb](upload_to_SRA_pan-sarbeco.ipynb) have instructions and does the uploading.

Because the FTP upload takes a while, you may want to run the Jupyter notebook using `slurm` so there is no timeout with::

    sbatch --wrap="jupyter nbconvert --to notebook --execute --inplace --ExecutePreprocessor.timeout=-1 upload_to_SRA.ipynb" --time 2-0
