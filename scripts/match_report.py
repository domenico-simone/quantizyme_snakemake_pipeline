#!/usr/bin/env python

import datetime, random, sys, os
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--html-report",            help="HTML output")
#parser.add_argument("--md-report",              help="Markdown output")
parser.add_argument("--css",                    help="CSS stylesheet for HTML output")
parser.add_argument("--report-template",        help="Report template")
parser.add_argument("-s", "--subtrees",         help="Number of subtrees in the used model")
parser.add_argument("-t", "--trials",           help="Number of subsamplings in the used model")
#parser.add_argument("--results-directory",      help="Directory with results")
parser.add_argument("--venn",                   help="Venn diagram of results overlap")

args = parser.parse_args()

html_report = args.html_report
md_report = args.html_report.replace(".html", ".md")
css = args.css
report_template = args.report_template
subtrees = int(args.subtrees)
trials = int(args.trials)
#results_dir = args.results_directory
venn = args.venn
results_dir = os.path.split(venn)[0]
# print(args)
# print(html_report)

### markdown hit table template data
# The report template is in a separate file.
# The only exception is the hit table which is built in this script
# and attached afterwards.
simple_cell = "| {}"
header_cell = "| trial{trial} "
cell = "| [{n_hits}][{subtree}_{trial}] "
header_line = "|:---:|"

table_ref_element = "[{subtree}_{trial}]: subtree_{subtree}_trial_{trial}_matched_reads.txt"
###

# get today's date
date = datetime.datetime.now()
today = "{year} - {month} - {day}\n".format(year = date.year, month = date.month, day = date.day)

# # table generation
# # subtrees are rows, trials are cols
# t = np.empty(shape = (subtrees, trials))
# 
# for i in range(trials):
#     for j in range(subtrees):
#         t[j,i] = random.randint(3000,5000)
# 
# print(t)

table_string = ""
table_refs = []

# header
for j in range(trials + 1):
    if j == 0:
        table_string += simple_cell.format("")
    else:
        table_string += header_cell.format(trial = j)

for k in ['mean', 'sd', 'cv']:
    table_string += simple_cell.format(k)

table_string += "|\n"

# header line
table_string += '|' + '|'.join([':---:']*(trials+4)) + '|\n'

for j in range(subtrees):
    n_hits_list = []
    for i in range(trials+1):
        if i == 0:
            table_string += simple_cell.format("**subtree{}**".format(j+1))
        else:
            result_file = "{}/subtree_{}_trial_{}_matched_reads.txt".format(results_dir, j+1, i)
            n_hits = len(open(result_file, 'r').readlines())
            table_string += cell.format(n_hits = n_hits, subtree = j+1, trial = i)
            n_hits_list.append(n_hits)
            table_refs.append(table_ref_element.format(subtree = j+1, trial = i))
    table_string += simple_cell.format(np.mean(n_hits_list))
    table_string += simple_cell.format(np.std(n_hits_list))
    table_string += simple_cell.format(np.std(n_hits_list)/np.mean(n_hits_list))
    table_string += "|\n"

# # get number of hits
# n_hits = len(open("results/hits.txt", 'r').readlines())

# report file (Markdown format)
report_template_txt = open(report_template, 'r').read()
r = open(md_report, 'w')
r.write(report_template_txt.format(date = date, subtrees = subtrees, trials = trials, table = table_string, table_refs = '\n'.join(table_refs)))
r.close()

# convert Markdown in HTML with pandoc
print("pandoc -f markdown -t html -o {output} -s --css={css} {input}".format(output = html_report, css = css, input = md_report))
os.system("pandoc -f markdown -t html -o {output} -s --css={css} -H {css} {input}".format(output = html_report, css = css, input = md_report))