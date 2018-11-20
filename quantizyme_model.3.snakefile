import pandas as pd

configfile: "config.2.yaml"

pipeline_dir = config['pipeline_dir']
scriptsDir = pipeline_dir + "/scripts"
reference_transcripts_dir = config['reference_transcripts_dir']
model_dir = config['model_dir']
#projectID = config['projectID']
#nSubtrees = config['subtrees']
# subtrees = range(1,config['subtrees']+1)
# nTrials = config['nr.trials.random.picking']
# trials = range(1, config['nr.trials.random.picking']+1)
# remove_seqs = config['remove_seqs']
# remove_lower_t = config['remove_lower_t']
# remove_higher_t = config['remove_higher_t']
# subgroup_percent = config['subgroup.percent']

def get_other_fields(df, ref_genome_mt, field):
    return list(df.loc[df['ref_genome_mt'] == ref_genome_mt, field])

analysis_tab = pd.read_table("analysis.tab", sep = "\t", comment='#')

def get_out_files(df, model_dir="models"):
    outpaths = []
    for row in df.itertuples():
        outpaths.append("{}/{}_MODEL_{}_{}_subtrees_{}_trials_{}_subgroup_{}_MODEL.tar.gz".format(model_dir, \
                                                                                getattr(row, "projectID"), \
                                                                                getattr(row, "remove_lower_t"), \
                                                                                getattr(row, "remove_higher_t"), \
                                                                                getattr(row, "subtrees"), \
                                                                                getattr(row, "nr_trials_random_picking"), \
                                                                                getattr(row, "subgroup_percent")))
        #outpaths.append("{}/OUT_{}_{}_{}/{}/OUT.sam".format(res_dir, getattr(row, "sample"), getattr(row, "ref_genome_mt"), getattr(row, "ref_genome_n"), map_dir))
    return outpaths

localrules: all, transcript_length_distribution, transcript_filtering, clustering, subtreeing1, compress_out_folder

rule all:
    input:
        #expand("reference_transcripts_length_distribution/{projectID}_transcript_length_distribution.pdf", projectID=projectID)
        model_archives = get_out_files(analysis_tab, model_dir=model_dir)

# rule fake:
#     input:
#         reference_transcripts = reference_transcripts_dir + "/{projectID}.fasta"
#     output:
#         model_archives = model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}_subtrees_{subtrees}_trials_{nr_trials_random_picking}_subgroup_{subgroup_percent}_MODEL.tar.gz"
#     shell:
#         """
#         # do sntg
#         """

# # rule all:
# #     input:
# #         expand("reference_transcripts_length_distribution/{projectID}_transcript_length_distribution.pdf", projectID=projectID),
# #         expand(model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}_subtrees_" + str(nSubtrees) + "_trials_" + str(nTrials) + "_MODEL.tar.gz", projectID=projectID, remove_lower_t=remove_lower_t, remove_higher_t=remove_higher_t)
# #         #model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/{projectID}_MODEL.zip"
# #         #expand(model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/subtrees_" + str(nSubtrees) + "/{projectID}_subtree_{n}.fasta", projectID=projectID, remove_lower_t=remove_lower_t, remove_higher_t=remove_higher_t, n=subtrees),
# #         #expand(model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/{projectID}_MODEL.zip", projectID=projectID, remove_lower_t=remove_lower_t, remove_higher_t=remove_higher_t)
#
# rule transcript_length_distribution:
#     input:
#         file=reference_transcripts_dir + "/{projectID}.fasta"
#     output:
#         file="reference_transcripts_length_distribution/{projectID}_transcript_length_distribution.pdf"
#     log:
#         "logs/quantizyme_model_figure1/{projectID}_quantizyme_model_figure1.log"
#     shell:
#         """
#         Rscript --vanilla {scriptsDir}/quantizyme_model_figure1.R {output.file} {input.file} 2> {log}
#         """
#
# rule transcript_filtering:
#     input:
#         file = reference_transcripts_dir + "/{projectID}.fasta"
#     output:
#         file= model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/{projectID}_start.fasta",
#         outplot= model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/{projectID}_thresholded_distribution.pdf"
#     params:
#         filter_seqs = {remove_seqs},
#         filter_lower_t = {remove_lower_t},
#         filter_higher_t = {remove_higher_t}
#     log:
#         "logs/quantizyme_model_transcript_filtering/{projectID}_filtering_{remove_lower_t}_{remove_higher_t}.log"
#     shell:
#         """
#         #!/bin/bash
#
#         export outdir=$(dirname {output.file})
#         echo ${{outdir}}
#         Rscript --vanilla {scriptsDir}/quantizyme_model_ref_sequence_filtering.R -d ${{outdir}} -p {projectID} -i {input.file} -o {output.file} -f {params.filter_seqs} -l {params.filter_lower_t} -u {params.filter_higher_t} --outplot={output.outplot} 2> {log}
#         """
#
# rule run_MSA:
#     input:
#         file=model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/{projectID}_start.fasta"
#     output:
#         alignment= model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/align1.phy"
#     # threads: cluster_config['run_MSA']['threads']
#     threads: 20
#     log:
#         "logs/quantizyme_model_MSA/{projectID}_filtering_{remove_lower_t}_{remove_higher_t}.log"
#     shell:
#         """
#         clustalo \
#         --threads {threads} \
#         -i {input} \
#         -o {output.alignment} \
#         --outfmt=phy  > {log} 2>&1
#         """
#
# rule clustering:
#     input:
#         alignment= model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/align1.phy"
#     output:
#         env_cluster_file= model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/envcluster.RData"
#     params:
#         projectID = "{projectID}",
#         fasta = model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/{projectID}_start.fasta"
#         #subtrees = "{n}"
#     shell:
#         """
#         outdir=$(dirname {input.alignment})
#         Rscript --vanilla {scriptsDir}/quantizyme_model_clustering.R -a {input.alignment} -i {params.fasta} -d ${{outdir}} -p {params.projectID} -e {output.env_cluster_file}
#         """
#
# rule subtreeing1:
#     input:
#         alignment = model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/align1.phy",
#         env_cluster_file = model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/envcluster.RData",
#         # alignment = expand(model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/align1.phy", projectID=projectID, remove_lower_t=remove_lower_t, remove_higher_t=remove_higher_t),
#         # env_cluster_file = expand(model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/envcluster.RData", projectID=projectID, remove_lower_t=remove_lower_t, remove_higher_t=remove_higher_t)
#     output:
#         subtree_fasta = model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/subtrees_" + str(nSubtrees) + "/{projectID}_subtree_{n}.fasta",
#         #subtree_fasta = expand(model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/subtrees_" + str(nSubtrees) + "/{projectID}_subtree_{n}.fasta", projectID=projectID, remove_lower_t=remove_lower_t, remove_higher_t=remove_higher_t, n=subtrees)
#         #model_params = model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/quantizyme_model_params.RData",
#     params:
#         fasta = model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/{projectID}_start.fasta",
#         projectID = "{projectID}",
#         subtrees = "{n}"
#         # fasta = expand(model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/{projectID}_start.fasta", projectID=projectID, remove_lower_t=remove_lower_t, remove_higher_t=remove_higher_t),
#         # projectID = expand("{projectID}", projectID=projectID),
#         # subtrees = expand("{n}", n=subtrees)
#     shell:
#         """
#         outdir=$(dirname {output.subtree_fasta})
#         Rscript --vanilla {scriptsDir}/quantizyme_model_cut_subtree.R -a {input.alignment} -i {params.fasta} -n {params.subtrees} -p {params.projectID} -e {input.env_cluster_file} -d ${{outdir}}
#         """
#
# rule run_MSA_subtree:
#     input:
#         subtree_fasta = model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/subtrees_" + str(nSubtrees) + "/{projectID}_subtree_{n}.fasta"
#         #file=model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/{projectID}_start.fasta"
#     output:
#         alignment= model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/subtrees_" + str(nSubtrees) + "/{projectID}_subtree_{n}.phy"
#     # threads: cluster_config['run_MSA_subtree']['threads']
#     threads: 10
#     log:
#         "logs/quantizyme_model_MSA/{projectID}_filtering_{remove_lower_t}_{remove_higher_t}_subtrees_" + str(nSubtrees) + "_{n}.log"
#     shell:
#         """
#         clustalo \
#         --threads {threads} \
#         -i {input} \
#         -o {output.alignment} \
#         --outfmt=phy  > {log} 2>&1
#         """
#
# rule subtreeing2:
#     input:
#         alignment= model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/subtrees_" + str(nSubtrees) + "/{projectID}_subtree_{n}.phy",
#         subtree_fasta = model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/subtrees_" + str(nSubtrees) + "/{projectID}_subtree_{n}.fasta",
#         #env_cluster_file = model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/envcluster.RData",
#         # alignment = expand(model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/align1.phy", projectID=projectID, remove_lower_t=remove_lower_t, remove_higher_t=remove_higher_t),
#         # env_cluster_file = expand(model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/envcluster.RData", projectID=projectID, remove_lower_t=remove_lower_t, remove_higher_t=remove_higher_t)
#     output:
#         flag_file = model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/subtrees_" + str(nSubtrees) + "/trials_" + str(nTrials) + "/{projectID}_subtree_{n}_" + str(nTrials) + "_hmm.done",
#         #  projectID=projectID, remove_lower_t=remove_lower_t, remove_higher_t=remove_higher_t, n=subtrees, t=trials)
#         #subtree_fasta = model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/subtrees_" + str(nSubtrees) + "/{projectID}_subtree_{n}.fasta",
#         #subtree_fasta = expand(model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/subtrees_" + str(nSubtrees) + "/{projectID}_subtree_{n}.fasta", projectID=projectID, remove_lower_t=remove_lower_t, remove_higher_t=remove_higher_t, n=subtrees)
#         #model_params = model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/quantizyme_model_params.RData",
#     params:
#         # fasta = model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/subtrees_" + str(nSubtrees) + "{projectID}_subtree_{n}.fasta",
#         projectID = "{projectID}",
#         #projectID = lambda wildcards: {projectID},
#         subtrees = "{n}",
#         trials = {nTrials},
#         subgroup_percent = {subgroup_percent}
#         # fasta = expand(model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/{projectID}_start.fasta", projectID=projectID, remove_lower_t=remove_lower_t, remove_higher_t=remove_higher_t),
#         # projectID = expand("{projectID}", projectID=projectID),
#         # subtrees = expand("{n}", n=subtrees)
#     # threads: cluster_config['subtreeing2']['threads']
#     threads: 15
#     shell:
#         """
#         outdir=$(dirname {output.flag_file})
#         Rscript --vanilla {scriptsDir}/quantizyme_model_subtree2.R -a {input.alignment} -i {input.subtree_fasta} -n {params.subtrees} -t {params.trials} -p {params.projectID} -d ${{outdir}} -s {params.subgroup_percent} -z {output.flag_file} --clustalo_threads {threads} --hmmbuild_threads {threads}
#         """
#
# # rule compress_out_folder:
# #     input:
# #         hmm_models = expand(model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/subtrees_" + str(nSubtrees) + "/{projectID}_subtree_{n}_{t}.hmm", projectID=projectID, remove_lower_t=remove_lower_t, remove_higher_t=remove_higher_t, n=subtrees, t=trials)
# #         #model_params = expand(model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/subtrees_" + str(nSubtrees) + "/{projectID}_subtree_{n}_{t}.hmm", projectID=projectID, remove_lower_t=remove_lower_t, remove_higher_t=remove_higher_t, n=subtrees, t=trials)
# #         #model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/{projectID}_annot_subtree_{n}.fasta"
# #         #expand(model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/{projectID}_annot_subtree_{n}.fasta", projectID=projectID, n=subtrees, remove_lower_t=remove_lower_t, remove_higher_t=remove_higher_t)
# #     output:
# #         #expand(model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/{projectID}_MODEL.zip", projectID=projectID, remove_lower_t=remove_lower_t, remove_higher_t=remove_higher_t)
# #         model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/{projectID}_MODEL.zip"
# #         #expand(model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/{projectID}_MODEL.zip"
# #         #"out.zip"
# #     shell:
# #         """
# #         gzip -c {input} > {output}
# #         """
#
rule compress_out_folder:
    input:
        lambda wildcards: expand(model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/subtrees_{nSubtrees}/trials_{nTrials}/{projectID}_subtree_{n}_{nTrials}_{subgroup_percent}_hmm.done", \
                                projectID = wildcards.projectID, remove_lower_t = wildcards.remove_lower_t, remove_higher_t = wildcards.remove_higher_t,
                                nSubtrees = nSubtrees, nTrials = nTrials, n = subtrees, subgroup_percent = subgroup_percent),
        #model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/subtrees_" + str(nSubtrees) + "/trials_" + str(nTrials) + "/{projectID}_subtree_{n}_" + str(nTrials) + "_hmm.done"
        #flag_file = expand(model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/subtrees_" + str(nSubtrees) + "/trials_" + str(nTrials) + "/{projectID}_subtree_{n}_" + str(nTrials) + "_hmm.done", projectID=projectID, remove_lower_t=remove_lower_t, remove_higher_t=remove_higher_t, n=subtrees)
        #model_params = expand(model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/subtrees_" + str(nSubtrees) + "/{projectID}_subtree_{n}_{t}.hmm", projectID=projectID, remove_lower_t=remove_lower_t, remove_higher_t=remove_higher_t, n=subtrees, t=trials)
        #model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/{projectID}_annot_subtree_{n}.fasta"
        #expand(model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/{projectID}_annot_subtree_{n}.fasta", projectID=projectID, n=subtrees, remove_lower_t=remove_lower_t, remove_higher_t=remove_higher_t)
    output:
        model_archives = model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}_subtrees_{subtrees}_trials_{nr_trials_random_picking}_subgroup_{subgroup_percent}_MODEL.tar.gz"
        #expand(model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/{projectID}_MODEL.zip", projectID=projectID, remove_lower_t=remove_lower_t, remove_higher_t=remove_higher_t)
        #model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}_subtrees_" + str(nSubtrees) + "_trials_" + str(nTrials) + "_MODEL.tar.gz"
        #expand(model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/{projectID}_MODEL.zip"
        #"out.zip"
    # params:
    #     hmm_models = expand(model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/subtrees_" + str(nSubtrees) + "/{projectID}_subtree_{n}_{t}.hmm", projectID=projectID, remove_lower_t=remove_lower_t, remove_higher_t=remove_higher_t, n=subtrees, t=trials)
    shell:
        """
        destDir=`pwd`
        fileDir=$(dirname {input.flag_file[0]})
        cd ${{fileDir}}
        tar -cvzf ${{destDir}}/{output} *.hmm
        """
