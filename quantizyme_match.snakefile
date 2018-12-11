import pandas as pd

configfile: "config.match.yaml"
model_dir = config['model_dir']
res_dir = config['results_dir']
pipeline_dir = config['pipeline_dir']
scriptsDir = pipeline_dir + "/scripts"

# Steps:
# - ...
# - dataset_model_Venn_diagram.png

def get_venn_files(df, res_dir="match"):
    outpaths = []
    for row in df.itertuples():
        dataset = getattr(row, "dataset")
        model_base = getattr(row, "model").replace("_MODEL.tar.gz", "")
        outpaths.append("{}/{}-{}/Venn_diagram.jpg".format(res_dir, dataset, model_base))
    return outpaths

def get_report_files(df, res_dir="match"):
    outpaths = []
    for row in df.itertuples():
        dataset = getattr(row, "dataset")
        model_base = getattr(row, "model").replace("_MODEL.tar.gz", "")
        outpaths.append("{}/{}-{}/{}-{}_match_report.html".format(res_dir, dataset, model_base, dataset, model_base))
    return outpaths

# def get_model_params(model_name):
#     params = {}
#     model_name_split = model_name.split('_MODEL_')
#     params['projectID'] = model_name_split[0]
#     params['remove_lower_t'] = model_name_split[1].split('_')[0]
#     params['remove_higher_t'] = model_name_split[1].split('_')[1]
#     params['subtrees'] = model_name_split[1].split('_')[3]
#     params['nr_trials_random_picking'] = model_name_split[1].split('_')[5]
#     params['subgroup_percent'] = model_name_split[1].split('_')[7]
#     return params

analysis_match_tab = pd.read_table("analysis.match.tab", sep = "\t", comment='#')

localrules: process_nhmmer_results

rule all:
    input:
        venn = get_venn_files(analysis_match_tab, res_dir = res_dir),
        html_report = get_report_files(analysis_match_tab, res_dir = res_dir)

rule nhmmer:
    input:
        dataset_file = "read_datasets/{dataset}.fasta",
        model_file = model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/subtrees_{subtrees}/trials_{nr_trials_random_picking}_subgroup_{subgroup_percent}/subtree_{n}_trial_{t}_profile.hmm"
    output:
        nhmmer_outputs = res_dir + "/{dataset}-{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}_subtrees_{subtrees}_trials_{nr_trials_random_picking}_subgroup_{subgroup_percent}/subtree_{n}_trial_{t}_subgroup_{subgroup_percent}_matches.tbl"
    params:
        nhmmerEvalue = 1,
        nhmmer_alignment = res_dir + "/{dataset}-{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}_subtrees_{subtrees}_trials_{nr_trials_random_picking}_subgroup_{subgroup_percent}/subtree_{n}_trial_{t}_subgroup_{subgroup_percent}_nhmmer.msa"
    threads: 15
    log:
        "logs/{dataset}-{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}_subtrees_{subtrees}_trials_{nr_trials_random_picking}_subgroup_{subgroup_percent}/subtree_{n}_trial_{t}_subgroup_{subgroup_percent}_matches.log"
    shell:
        """
        nhmmer --cpu {threads} -o {log} -A {params.nhmmer_alignment} -E {params.nhmmerEvalue} --tblout {output.nhmmer_outputs} --tformat fasta {input.model_file} {input.dataset_file}
        """

rule process_nhmmer_results:
    input:
        nhmmer_outputs = lambda wildcards: expand(res_dir + "/{dataset}-{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}_subtrees_{subtrees}_trials_{nr_trials_random_picking}_subgroup_{subgroup_percent}/subtree_{n}_trial_{t}_subgroup_{subgroup_percent}_matches.tbl", \
                                                        dataset = wildcards.dataset, projectID = wildcards.projectID, remove_lower_t = wildcards.remove_lower_t, remove_higher_t = wildcards.remove_higher_t, \
                                                        t = range(1,int(wildcards.nr_trials_random_picking)+1), \
                                                        subtrees = wildcards.subtrees, nr_trials_random_picking = wildcards.nr_trials_random_picking, \
                                                        n = range(1,int(wildcards.subtrees)+1), subgroup_percent = wildcards.subgroup_percent)
    output:
        venn = res_dir + "/{dataset}-{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}_subtrees_{subtrees}_trials_{nr_trials_random_picking}_subgroup_{subgroup_percent}/Venn_diagram.jpg"
    params:
        subtrees = lambda wildcards: wildcards.subtrees,
        nr_trials_random_picking = lambda wildcards: wildcards.nr_trials_random_picking,
        subgroup_percent = lambda wildcards: wildcards.subgroup_percent
    shell:
        """
        outdir=$(dirname $(echo "{input.nhmmer_outputs}" | awk '{{print $1}}'))
        Rscript --vanilla {scriptsDir}/quantizyme_match_process_nhmmer_results.R -s {params.subtrees} -t {params.nr_trials_random_picking} -d ${{outdir}} -u {params.subgroup_percent} -o {output.venn}
        #echo ${{outdir}}
        #touch {output.venn}
        """

rule report:
    input:
        venn = res_dir + "/{dataset}-{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}_subtrees_{subtrees}_trials_{nr_trials_random_picking}_subgroup_{subgroup_percent}/Venn_diagram.jpg"
    output:
        html_report = res_dir + "/{dataset}-{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}_subtrees_{subtrees}_trials_{nr_trials_random_picking}_subgroup_{subgroup_percent}/{dataset}-{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}_subtrees_{subtrees}_trials_{nr_trials_random_picking}_subgroup_{subgroup_percent}_match_report.html"
    params:
        subtrees = lambda wildcards: wildcards.subtrees,
        trials = lambda wildcards: wildcards.nr_trials_random_picking,
        css = "{}/templates/github.css".format(pipeline_dir),
        report_template = "{}/templates/report_template.md".format(pipeline_dir)
    shell:
        """
        {scriptsDir}/match_report.py --html-report={output.html_report} \
            -s {params.subtrees} \
            -t {params.trials} \
            --venn={input.venn} \
            --css={params.css} \
            --report-template={params.report_template}
        """
