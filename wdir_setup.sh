#/bin/bash

if [[ $# -ne 1 ]]; then
    export wdir="quantizyme_runs"
    echo "The working directory \"${wdir}\" will be created as subdirectory of the current one."
else
    wdir=$1
fi

mkdir -p ${wdir}/reference_transcripts

ln -sf $(pwd -P)/quantizyme_model.3.snakefile ${wdir}
cp cluster.yaml ${wdir}
cp analysis.tab ${wdir}

echo "pipeline_dir: \"`pwd`\"" | cat config.yaml - > ${wdir}/config.yaml

echo "${wdir} created with this structure:"
tree ${wdir}