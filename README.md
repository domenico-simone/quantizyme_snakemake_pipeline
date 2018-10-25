# Quantizyme on the Mykopat grid

## Installation of conda environment

The grid implementation of the Mykopat pipeline is deployed in a conda environment, i.e. a virtual environment with all the needed tools/modules (including hmmer and clustal-omega). Installing conda is therefore essential, before installing the pipeline.

### Installation of Anaconda

Follow instructions at https://conda.io/docs/user-guide/install/linux.html

### Installation of the pipeline

```bash
# update git
conda install git

# fetch repo
git clone -b ultralight https://github.com/domenico-simone/quantizyme_snakemake_pipeline.git

# install environment
conda env create -n quantizyme_model -f environment.yaml
```

## Running the pipeline

Since the pipeline is not interactive, it can be useful to get the transcript length distribution first, display the plot and set up the most suitable parameters for

- `remove_seqs` (default: FALSE)
- `remove_lower_t` (default: 0)
- `remove higher_t` (default: 0)


### Get transcript distribution

```bash
export projectID=ref_AA2

snakemake -rp \
--config projectID=${projectID} \
-s quantizyme_model.snake \
reference_transcripts_length_distribution/${projectID}_transcript_length_distribution.pdf
```

### Run the whole MODEL pipeline

```bash
nohup snakemake -rpF \
--config projectID=${projectID} remove_seqs=TRUE remove_lower_t=200 remove_higher_t=1500 \
-s quantizyme_model.snake \
-j 50 --cluster-config cluster.json \
--cluster 'qsub -V -l h_rt=1:00:00 -pe smp 20 -cwd -j y' &
```
