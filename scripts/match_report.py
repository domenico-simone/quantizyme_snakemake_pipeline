#!/usr/bin/env python

import datetime, random, sys
import numpy as np
import argparse

def usage():
    print("""match_report.py <subtrees> <trials> <table_ref_element> """)
    sys.exit("\nYeah.")

parser = argparse.ArgumentParser()
parser.add_argument("--html-report",            help="HTML output")
parser.add_argument("--md-report",              help="Markdown output")
parser.add_argument("--css",                    help="CSS stylesheet for HTML output")
parser.add_argument("--report-template",        help="Report template")
parser.add_argument("-s", "--subtrees",         help="Number of subtrees in the used model")
parser.add_argument("-t", "--trials",           help="Number of subsamplings in the used model")
parser.add_argument("--results-directory",      help="Directory with results")
parser.parse_args()

usage()

# markdown template data
simple_cell = "| {}"
header_cell = "| trial{trial} "
cell = "| [{n_hits}][{subtree}_{trial}] "
header_line = "|:---:|"

table_ref_element = "[{subtree}_{trial}]: results/subtree{subtree}_trial{trial}.txt"

# read report template
report_template = open("sample.md", 'r').read()

# get today's date
date = datetime.datetime.now()
today = "{year} - {month} - {day}\n".format(year = date.year, month = date.month, day = date.day)

# get subtrees and trials
(subtrees, trials) = [int(i.split()[1]) for i in open("results/parameters.txt", 'r').readlines()]

# table generation
# subtrees are rows, trials are cols
t = np.empty(shape = (subtrees, trials))

for i in range(trials):
    for j in range(subtrees):
        t[j,i] = random.randint(3000,5000)

print(t)

table_string = ""
table_refs = []

# header
for j in range(trials + 1):
    if j == 0:
        table_string += simple_cell.format("")
    else:
        table_string += header_cell.format(trial = j)
table_string += "|\n"

# header line
table_string += '|' + '|'.join([':---:']*(trials+1)) + '|\n'

for j in range(subtrees):
    for i in range(trials+1):
        if i == 0:
            table_string += simple_cell.format("**subtree{}**".format(j+1))
        else:
            table_string += cell.format(n_hits = t[j,i-1], subtree = j+1, trial = i)
            table_refs.append(table_ref_element.format(subtree = j+1, trial = i))
    table_string += "|\n"

# get number of hits
n_hits = len(open("results/hits.txt", 'r').readlines())

# open report file
r = open('sample_compiled.md', 'w')
r.write(report_template.format(date = date, subtrees = subtrees, trials = trials, hits = n_hits, table = table_string, table_refs = '\n'.join(table_refs)))

r.close()
