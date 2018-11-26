import pandas as pd

configfile: "config.match.yaml"
res_dir = config['results_dir']

# Steps:
# - ...
# - dataset_model_Venn_diagram.png

def get_venn_files(df, res_dir="match"):
    outpaths = []
    for row in df.itertuples():
        dataset = getattr(row, "dataset")
        model_base = getattr(row, "model").replace("_MODEL.tar.gz", "")
        outpaths.append("{}/{}_{}/Venn_diagram.png".format(res_dir, dataset, model_base))
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

rule all:
    input:
        venn = get_venn_files(analysis_match_tab, res_dir = res_dir)

rule nhmmer:
    input:
        dataset_file = "read_datasets/{dataset}.fasta"
        model_file = model_dir + "/{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}/subtrees_{subtrees}/trials_{nr_trials_random_picking}_subgroup_{subgroup_percent}/{projectID}_subtree_{n}_trials_{nr_trials_random_picking}_profile.hmm"
    output:
        nhmmer_outputs = res_dir + "/{dataset}_{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}_subtrees_{subtrees}_trials_{nr_trials_random_picking}_subgroup_{subgroup_percent}/subtree_{n}_trials_{nr_trials_random_picking}_subgroup_{subgroup_percent}_matches.tbl"

rule process_nhmmer_results:
    input:
        nhmmer_outputs = lambda wildcards: expand(res_dir + "/{dataset}_{projectID}_MODEL_{remove_lower_t}_{remove_higher_t}_subtrees_{subtrees}_trials_{nr_trials_random_picking}_subgroup_{subgroup_percent}/subtree_{n}_trials_{nr_trials_random_picking}_subgroup_{subgroup_percent}_matches.tbl", \
                                                        dataset = wildcards.dataset, projectID = wildcards.projectID, remove_lower_t = wildcards.remove_lower_t, remove_higher_t = wildcards.remove_higher_t, \
                                                        subtrees = wildcards.subtrees, nr_trials_random_picking = wildcards.nr_trials_random_picking, \
                                                        n = range(1,int(wildcards.subtrees)+1), subgroup_percent = wildcards.subgroup_percent)
    output:
        venn = get_venn_files(analysis_match_tab, res_dir = res_dir)
    shell: "touch bwe"
