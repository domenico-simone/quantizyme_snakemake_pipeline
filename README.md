# Quantizyme on the Mykopat grid

## Installation of conda environment

The grid implementation of the Mykopat pipeline is deployed in a conda environment, i.e. a virtual environment with all the needed tools/modules (including hmmer and clustal-omega). Installing conda is therefore essential, before installing the pipeline.

### Installation of Anaconda

Follow instructions at https://conda.io/docs/user-guide/install/linux.html

### Installation of the pipeline

```bash
### update git
conda install git

### fetch repo
# shouldn't be harder than this but git gets stuck
# git clone -b ultralight https://github.com/domenico-simone/quantizyme_snakemake_pipeline.git
# so use this more complicated procedure
curl -OL https://github.com/domenico-simone/quantizyme_snakemake_pipeline/archive/ultralight.zip
unzip ultralight.zip
mv quantizyme_snakemake_pipeline-ultralight quantizyme_snakemake_pipeline

# install environment
cd quantizyme_snakemake_pipeline
conda env create -n quantizyme_model -f environment.yaml

# activate environment
source activate quantizyme_model
```

### Setup working directory

Setup a working directory:

- create input data folder (`wdir/reference_transcripts`)
- symlink some pipeline files

```bash
export PFOLDER=/path/to/pipeline/folder
cd /path/to/working/directory

# create input data folder
mkdir -p reference_transcripts

# create symlinks
ln -s ${PFOLDER}/quantizyme_model.2.snakefile .
ln -s ${PFOLDER}/cluster.yaml .
ln -s ${PFOLDER}/config.2.yaml .
```

## Running the pipeline

To run the pipeline, you need to use R/3.2.0.

```bash
# load appropriate R version
module load R/3.2.0
```

Since the pipeline is not interactive, it can be useful to get the transcript length distribution first, display the plot and set up the most suitable parameters for

- `remove_seqs` (default: FALSE)
- `remove_lower_t` (default: 0)
- `remove higher_t` (default: 0)


### Get transcript distribution

```bash
snakemake -prF --until transcript_length_distribution -s quantizyme_model.3.snakefile
```

### Run the whole MODEL pipeline

```bash
nohup snakemake -prF \
-j 100 --cluster-config cluster.yaml --cluster 'qsub -V -l h_rt=3:00:00 -pe smp {cluster.threads} -cwd -j y' \
-s quantizyme_model.3.snakefile &> quantizyme_model.log &
```

### Get the DAG representation of the pipeline

```bash
snakemake --dag \
-s quantizyme_model.3.snakefile | dot -Tsvg > quantizyme_model.dag.svg
```
