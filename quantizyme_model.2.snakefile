configfile: "config.2.yaml"

pipeline_dir = config['pipeline_dir']
scriptsDir = pipeline_dir + "/scripts"
reference_transcripts_dir = config['reference_transcripts_dir']
model_dir = config['model_dir']
projectID = config['projectID']
subtrees = range(1,config['subtrees']+1)
remove_seqs = config['remove_seqs']
remove_lower_t = config['remove_lower_t']
remove_higher_t = config['remove_higher_t']

localrules: all, transcript_length_distribution, transcript_filtering, run_MSA, clustering, subtreeing1, compress_out_folder

rule all:
    input:
        expand("reference_transcripts_length_distribution/{projectID}_transcript_length_distribution.pdf", projectID=projectID),
        #model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/{projectID}_MODEL.zip"
        expand(model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/{projectID}_annot_subtree_{n}.fasta", projectID=projectID, remove_lower_t=remove_lower_t, remove_higher_t=remove_higher_t, n=subtrees),
        expand(model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/{projectID}_MODEL.zip", projectID=projectID, remove_lower_t=remove_lower_t, remove_higher_t=remove_higher_t)


rule transcript_length_distribution:
    input:
        file=reference_transcripts_dir + "/{projectID}.fasta"
    output:
        file="reference_transcripts_length_distribution/{projectID}_transcript_length_distribution.pdf"
    log:
        "logs/quantizyme_model_figure1/{projectID}_quantizyme_model_figure1.log"
    shell:
        """
        Rscript --vanilla {scriptsDir}/quantizyme_model_figure1.R {output.file} {input.file} 2> {log}
        """

rule transcript_filtering:
    input:
        file = reference_transcripts_dir + "/{projectID}.fasta"
    output:
        file= model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/{projectID}_start.fasta",
        outplot= model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/{projectID}_thresholded_distribution.pdf"
    params:
        filter_seqs = {remove_seqs},
        filter_lower_t = {remove_lower_t},
        filter_higher_t = {remove_higher_t}
    log:
        "logs/quantizyme_model_transcript_filtering/{projectID}_filtering_{remove_lower_t}_{remove_higher_t}.log"
    shell:
        """
        #!/bin/bash

        export outdir=$(dirname {output.file})
        echo ${{outdir}}
        Rscript --vanilla {scriptsDir}/quantizyme_model_ref_sequence_filtering.R -d ${{outdir}} -p {projectID} -i {input.file} -o {output.file} -f {params.filter_seqs} -l {params.filter_lower_t} -u {params.filter_higher_t} --outplot={output.outplot} 2> {log}
        """

rule run_MSA:
    input:
        file=model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/{projectID}_start.fasta"
    output:
        alignment= model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/align1.phy"
    threads: 40
    log:
        "logs/quantizyme_model_MSA/{projectID}_filtering_{remove_lower_t}_{remove_higher_t}.log"
    shell:
        """
        clustalo \
        --threads {threads} \
        -i {input} \
        -o {output.alignment} \
        --outfmt=phy  > {log} 2>&1
        """

rule clustering:
    input:
        alignment= model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/align1.phy"
    output:
        env_cluster_file= model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/envcluster.RData"
    params:
        projectID = "{projectID}",
        fasta = model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/{projectID}_start.fasta"
    shell:
        """
        outdir=$(dirname {input.alignment})
        Rscript --vanilla {scriptsDir}/quantizyme_model_clustering.R -a {input.alignment} -i {params.fasta} -d ${{outdir}} -p {params.projectID} -e {output.env_cluster_file}
        """

rule subtreeing1:
    input:
        alignment = model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/align1.phy",
        env_cluster_file = model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/envcluster.RData"
    output:
        subtree_fasta = model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/{projectID}_annot_subtree_{n}.fasta",
        #model_params = model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/quantizyme_model_params.RData",
    params:
        fasta = model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/{projectID}_start.fasta",
        projectID = "{projectID}",
        subtrees = "{n}"
    shell:
        """
        outdir=$(dirname {input.alignment})
        Rscript --vanilla {scriptsDir}/quantizyme_model_cut_subtree.R -a {input.alignment} -i {params.fasta} -n {params.subtrees} -e {input.env_cluster_file}
        """

rule compress_out_folder:
    input:
        model_params = expand(model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/{projectID}_annot_subtree_{n}.fasta", projectID=projectID, remove_lower_t=remove_lower_t, remove_higher_t=remove_higher_t, n=subtrees)
        #model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/{projectID}_annot_subtree_{n}.fasta"
        #expand(model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/{projectID}_annot_subtree_{n}.fasta", projectID=projectID, n=subtrees, remove_lower_t=remove_lower_t, remove_higher_t=remove_higher_t)
    output:
        #expand(model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/{projectID}_MODEL.zip", projectID=projectID, remove_lower_t=remove_lower_t, remove_higher_t=remove_higher_t)
        model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/{projectID}_MODEL.zip"
        #expand(model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/{projectID}_MODEL.zip"
        #"out.zip"
    shell:
        """
        gzip -c {input} > {output}
        """
