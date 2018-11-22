# Quantizyme on the Mykopat grid

## Installation of conda environment

The grid implementation of the quantizyme workflow is deployed in a conda environment, *i.e.* a virtual environment with all the needed tools/modules (including hmmer and clustal-omega). Installing conda is therefore essential, before installing the workflow.

### Installation of Anaconda

Follow instructions at https://conda.io/docs/user-guide/install/linux.html

### Installation of the workflow

#### Standard installation

Clone repo and create working directory

```bash
git clone -b analysis_tab https://github.com/domenico-simone/quantizyme_snakemake_pipeline.git
cd quantizyme_snakemake_pipeline

# create conda env based on environment.yaml file
conda env create -n quantizyme_model -f environment.yaml

# activate environment
source activate quantizyme_model
```

#### Git clone gets stuck?

```bash
### update git
conda install git

curl -OL https://github.com/domenico-simone/quantizyme_snakemake_pipeline/archive/analysis_tab.zip
unzip analysis_tab.zip
mv quantizyme_snakemake_pipeline-analysis_tab quantizyme_snakemake_pipeline

# install environment
cd quantizyme_snakemake_pipeline
conda env create -n quantizyme_model -f environment.yaml
```

### Setup a working directory

Create working directory `path/to/working/directory` and copy/symlink files needed to run the workflow. This is automatically performed by the script `wdir_setup.sh`, which takes as single argument the path of the directory where you want to run the workflow. If the directory name is not provided, the working directory `quantizyme_runs` will be created as subdirectory of the current one.

```bash
# setup working directory
bash wdir_setup.sh /path/to/working/directory
```

Assuming a directory name like `test2`, the structure of the newly created directory should be

```bash
test2
├── analysis.tab
├── cluster.yaml
├── config.yaml
├── quantizyme_model.3.snakefile -> /home/adm2/Research/slu/quantizyme_snakemake_pipeline/quantizyme_model.3.snakefile
└── reference_transcripts

```

## Running the workflow

Go to the working directory you have created and activate the conda environment

```bash
cd /path/to/working/directory
source activate quantizyme_model
```

To run the workflow, you need to use R/3.2.0.

```bash
# load appropriate R version
module load R/3.2.0
```

Since the workflow is not interactive, it can be useful to get the transcript length distribution first, display the plot and set up the most suitable parameters for

- `remove_seqs` (default: FALSE)
- `remove_lower_t` (default: 0)
- `remove higher_t` (default: 0)

### Tell the workflow what to do: compile the `analysis.tab` datasheet

The file `analysis.tab` will be used by the workflow to infer all the parameters used in the generation of gene models. For each model, one row defines the parameters used. A sample of the structure of the `analysis.tab` file is provided as follows:

| projectID | remove_seqs | remove_lower_t | remove_higher_t | subtrees | nr_trials_random_picking | subgroup_percent |
|:---------:|:-----------:|:--------------:|:---------------:|:--------:|:------------------------:|:----------------:|
| ref_AA2   | TRUE        | 200            | 1500            | 3        | 10                       | 30               |
| ref_AA3   | TRUE        | 200            | 1600            | 3        | 10                       | 30               |
| ref_AA2   | TRUE        | 300            | 1500            | 3        | 10                       | 30               |
| ref_AA2   | TRUE        | 300            | 1500            | 4        | 10                       | 30               |
| ref_AA2   | TRUE        | 300            | 1500            | 4        | 10                       | 40               |
| ref_AA3   | FALSE       | 200            | 1600            | 3        | 10                       | 30               |

In this case, the workflow will compute five models, based on the specified parameters. You can edit the `analysis.tab` directly as text file, or copy it in a spreadsheet (eg in Excel) and exporting it back as text file (with **tab** separator).

**Note on parameters**: if you set the `remove_seqs` option as FALSE, the `remove_lower_t` and `remove_higher_t` options will make no difference in the setup.

### Get the DAG representation of the workflow

Before executing the workflow, it could be a good idea to display a diagram of it. The following command creates a plot of your “directed acyclic graph” (namely DAG, a plot of all of the rules Snakemake will execute), which you can view using any picture viewing program.

```bash
snakemake --dag \
-s quantizyme_model.3.snakefile | dot -Tsvg > quantizyme_model.dag.svg
```

### Get transcript distribution

```bash
snakemake -pr --until transcript_length_distribution -s quantizyme_model.3.snakefile
```

### Run the whole MODEL workflow

```bash
nohup snakemake -pr \
-j 100 --cluster-config cluster.yaml --cluster 'qsub -V -l h_rt=3:00:00 -pe smp {cluster.threads} -cwd -j y' \
-s quantizyme_model.3.snakefile &> quantizyme_model.log &
```

#### Notes on the execution and re-execution of the workflow...

... and, as a consequence, on the beauty of Snakemake :-) .

- It is not possible to have the workflow running more than once in the same directory. If you wish to do so, please set up another working directory (by following the instructions)
- If you wish to compute a model on the same gene but with different parameters, you can just add it to a pre-existing `analysis.tab` file. Assuming you are keeping all the output files of previous runs, the workflow will automatically understand which files need to be generated. This means that, if you are computing a model for which only the number of subtrees is changed (compared to another model already computed), the workflow won't need to perform the first multiple sequence alignment again (which is usually the most time consuming step!!!), since it will re-use the one previously performed.