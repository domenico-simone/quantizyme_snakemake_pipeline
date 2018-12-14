import pandas as pd

configfile: "config.yaml"

pipeline_dir = config['pipeline_dir']
#print(pipeline_dir)
scriptsDir = pipeline_dir + "/scripts"
reference_transcripts_dir = config['reference_transcripts_dir']
model_dir = config['model_dir']

def get_out_files(df, model_dir="models"):
    outpaths = []
    for row in df.itertuples():
        if str(getattr(row, "remove_seqs")) == "True":
            outpaths.append("{}/{}_MODEL_{}_{}_subtrees_{}_trials_{}_subgroup_{}_MODEL.tar.gz".format(model_dir, \
                                                                                getattr(row, "projectID"), \
                                                                                getattr(row, "remove_lower_t"), \
                                                                                getattr(row, "remove_higher_t"), \
                                                                                getattr(row, "subtrees"), \
                                                                                getattr(row, "nr_trials_random_picking"), \
                                                                                getattr(row, "subgroup_percent")))
        elif str(getattr(row, "remove_seqs")) == "False":
            outpaths.append("{}/{}_MODEL_0_0_subtrees_{}_trials_{}_subgroup_{}_MODEL.tar.gz".format(model_dir, \
                                                                                getattr(row, "projectID"), \
                                                                                getattr(row, "subtrees"), \
                                                                                getattr(row, "nr_trials_random_picking"), \
                                                                                getattr(row, "subgroup_percent")))
    return outpaths

def get_transcript_length_distribution_plots(df, plot_dir="reference_transcripts_length_distribution"):
    outplots = []
    q = list(set(df['projectID']))
    for ref in q:
        outplots.append("{}/{}_transcript_length_distribution.pdf".format(plot_dir, ref))
    return outplots

def get_projectIDs(df):
    return list(set(df['projectID']))

def get_overall_trees(df):
    trees = []
    for row in df.itertuples():
        if str(getattr(row, "remove_seqs")) == "True":
            trees.append("{}/{}_MODEL_{}_{}/align1.tree.jpeg".format(model_dir, \
                                                                    getattr(row, "projectID"), \
                                                                    getattr(row, "remove_lower_t"), \
                                                                    getattr(row, "remove_higher_t")))
        elif str(getattr(row, "remove_seqs")) == "False":
            trees.append("{}/{}_MODEL_0_0/align1.tree.jpeg".format(model_dir, \
                                                                    getattr(row, "projectID")))
    return trees

analysis_tab = pd.read_table("analysis.tab", sep = "\t", comment='#')
#print(analysis_tab)

localrules: all, transcript_length_distribution, transcript_filtering, subtreeing1, compress_out_folder

rule all:
    input:
        transcript_length_distribution_plots = get_transcript_length_distribution_plots(analysis_tab),
        overall_trees = get_overall_trees(analysis_tab),
        model_archives = get_out_files(analysis_tab, model_dir=model_dir)

rule transcript_length_distribution:
    input:
        #file=expand(reference_transcripts_dir + "/{projectID}.fasta", projectID = get_projectIDs(analysis_tab))
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
        file = reference_transcripts_dir + "/{projectID}.fasta",
        inplot="reference_transcripts_length_distribution/{projectID}_transcript_length_distribution.pdf"
    output:
        file= model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/{projectID}_start.fasta",
        outplot= model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/{projectID}_thresholded_distribution.pdf"
    params:
        #filter_seqs = {remove_seqs},
        projectID = "{projectID}",
        filter_seqs = "TRUE",
        filter_lower_t = "{remove_lower_t}",
        filter_higher_t = "{remove_higher_t}"
    log:
        "logs/quantizyme_model_transcript_filtering/{projectID}_filtering_{remove_lower_t}_{remove_higher_t}.log"
    shell:
        """
        #!/bin/bash

        export outdir=$(dirname {output.file})
        echo ${{outdir}}
        
        if [[ {wildcards.remove_lower_t} -eq 0 && {wildcards.remove_higher_t} -eq 0 ]] 
        then
            cat {input.file} > {output.file}
            cp {input.inplot} {output.outplot}
        else
            Rscript --vanilla {scriptsDir}/quantizyme_model_ref_sequence_filtering.R -d ${{outdir}} -p {params.projectID} -i {input.file} -o {output.file} -f {params.filter_seqs} -l {params.filter_lower_t} -u {params.filter_higher_t} --outplot={output.outplot} 2> {log}
        fi
        """

rule run_MSA:
    input:
        file=model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/{projectID}_start.fasta"
    output:
        alignment= model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/align1.fasta"
    # threads: cluster_config['run_MSA']['threads']
    threads: 20
    log:
        "logs/quantizyme_model_MSA/{projectID}_filtering_{remove_lower_t}_{remove_higher_t}.log"
    shell:
        """
        clustalo \
        --threads {threads} \
        -i {input} \
        -o {output.alignment} \
        --outfmt=fa  > {log} 2>&1
        """

rule clustering:
    input:
        alignment= model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/align1.fasta"
    output:
        env_cluster_file= model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/envcluster.RData"
    params:
        projectID = "{projectID}",
        fasta = model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/{projectID}_start.fasta"
        #subtrees = "{n}"
    shell:
        """
        outdir=$(dirname {input.alignment})
        Rscript --vanilla {scriptsDir}/quantizyme_model_clustering.R -a {input.alignment} -i {params.fasta} -d ${{outdir}} -p {params.projectID} -e {output.env_cluster_file}
        """

rule overall_tree:
    input:
        alignment = model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/align1.fasta"
    output:
        tree = model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/align1.tree"
    #threads: 3
    shell:
        """
        export OMP_NUM_THREADS={threads}
        FastTreeMP -noboot -nt -gtr \
            {input.alignment} > {output.tree}
        """

rule plot_overall_tree:
    input:
        tree = model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/align1.tree"
    output:
        plot_tree = model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/align1.tree.jpeg"
    params:
        projectID = lambda wildcards: wildcards.projectID
    shell:
        """
        quantizyme_model_plot_overall_tree.R --tree {input.tree} --tree_plot {output.plot_tree} --projectID {params.projectID}
        """

rule subtreeing1:
    input:
        #alignment = model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/align1.phy",
        env_cluster_file = model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/envcluster.RData",
    output:
        subtree_fasta = model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/subtrees_{subtrees}/{projectID}_subtree_{n}.fasta",
    params:
        fasta = model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/{projectID}_start.fasta",
        projectID = "{projectID}",
        subtrees = "{n}"
    shell:
        """
        outdir=$(dirname {output.subtree_fasta})
        Rscript --vanilla {scriptsDir}/quantizyme_model_cut_subtree.R -i {params.fasta} -n {params.subtrees} -p {params.projectID} -e {input.env_cluster_file} -d ${{outdir}}
        """

rule run_MSA_subtree:
    input:
        subtree_fasta = model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/subtrees_{subtrees}/{projectID}_subtree_{n}.fasta"
    output:
        alignment= model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/subtrees_{subtrees}/{projectID}_subtree_{n}.phy"
    threads: 10
    log:
        "logs/quantizyme_model_MSA/{projectID}_filtering_{remove_lower_t}_{remove_higher_t}_subtrees_{subtrees}_{n}.log"
    shell:
        """
        clustalo \
        --threads {threads} \
        -i {input} \
        -o {output.alignment} \
        --outfmt=phy  > {log} 2>&1
        """

rule subtreeing2:
    input:
        alignment = model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/subtrees_{subtrees}/{projectID}_subtree_{n}.phy",
        subtree_fasta = model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/subtrees_{subtrees}/{projectID}_subtree_{n}.fasta",
    output:
        flag_file = model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/subtrees_{subtrees}/trials_{nr_trials_random_picking}_subgroup_{subgroup_percent}/{projectID}_subtree_{n}_trials_{nr_trials_random_picking}_subgroup_{subgroup_percent}_hmm.done"
    params:
        projectID = lambda wildcards: wildcards.projectID,
        subtrees = "{n}",
        trials = "{nr_trials_random_picking}",
        subgroup_percent = "{subgroup_percent}"
    threads: 15
    shell:
        """
        outdir=$(dirname {output.flag_file})
        Rscript --vanilla {scriptsDir}/quantizyme_model_subtree2.R -a {input.alignment} -i {input.subtree_fasta} -n {params.subtrees} -t {params.trials} -p {params.projectID} -d ${{outdir}} -s {params.subgroup_percent} -z {output.flag_file} --clustalo_threads {threads} --hmmbuild_threads {threads}
        """

rule compress_out_folder:
    input:
        flag_file = lambda wildcards: expand(model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/subtrees_{subtrees}/trials_{nr_trials_random_picking}_subgroup_{subgroup_percent}/{projectID}_subtree_{n}_trials_{nr_trials_random_picking}_subgroup_{subgroup_percent}_hmm.done", \
                                projectID = wildcards.projectID, remove_lower_t = wildcards.remove_lower_t, remove_higher_t = wildcards.remove_higher_t, \
                                subtrees = wildcards.subtrees, nr_trials_random_picking = wildcards.nr_trials_random_picking, n = range(1,int(wildcards.subtrees)+1), subgroup_percent = wildcards.subgroup_percent),
    output:
        model_archives = model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}_subtrees_{subtrees}_trials_{nr_trials_random_picking}_subgroup_{subgroup_percent}_MODEL.tar.gz"
    shell:
        """
        destDir=`pwd`
        fileDir=$(dirname {input.flag_file[0]})
        cd ${{fileDir}}
        tar -cvzf ${{destDir}}/{output} *.hmm
        """
