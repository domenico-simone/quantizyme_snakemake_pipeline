#/bin/bash

if [[ $# -ne 1 ]]; then
    export wdir="quantizyme_runs"
    echo "The working directory \"${wdir}\" will be created as subdirectory of the current one."
else
    wdir=$1
fi

# MODEL
mkdir -p ${wdir}/reference_transcripts

ln -sf $(pwd -P)/quantizyme_model.3.snakefile ${wdir}
cp cluster.yaml ${wdir}
cp analysis.tab ${wdir}

echo "pipeline_dir: \"`pwd`\"" | cat config.yaml - > ${wdir}/config.yaml

# MATCH
mkdir -p ${wdir}/match
mkdir -p ${wdir}/read_datasets

ln -sf $(pwd -P)/quantizyme_match.snakefile ${wdir}
cp analysis.match.tab ${wdir}

echo "pipeline_dir: \"`pwd`\"" | cat config.match.yaml - > ${wdir}/config.match.yaml

echo "${wdir} created with this structure:"
tree ${wdir}

echo ""
echo "### MODEL workflow"
echo "Please put your input reference transcripts in the `reference_transcripts` dir."
echo "Distribution of lengths of reference transcripts will be located in the `reference_transcripts_length_distribution` folder."
echo ""
echo "### MATCH workflow"
echo "Please put your read datasets in fasta format in the `read_datasets` dir."
echo "Results of the MATCH workflow will be located in the `match` folder."