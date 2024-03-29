# vim: syntax=snakemake

import os
from os.path import join

configfile: "config.yaml"

if 'output' not in config:
    config['output'] = 'output'

cancer_types = config['cancer_types'].split(' ')

include: 'rules/preprocess.rules'
include: 'rules/diffExpr.rules'
include: 'rules/rabit.rules'
include: 'rules/misc.rules'

rule all:
    input:
        join(config['output'], 'rabit_formatted_results.txt')

#########################
# Preprocess mut/expression/copy number data
#########################
rule preprocess:
    input:
        expand(join(config['data_dir'], "mutation_flags/{ctype}.txt"), ctype=cancer_types),
        expand(join(config['data_dir'], "expr/{ctype}.txt"), ctype=cancer_types)

########################
# run diff expr analysis
########################
rule diffExpr:
    input:
        expand(join(config['output'], 'gene_correlations/{ctype}.txt'), ctype=cancer_types)

#########################
# Rabit analysis
#########################
rule rabit:
    input:
        expand(join(config['output'], 'rabit/output/{ctype}.txt'), ctype=cancer_types)

#########################
# Re-run rabit with correction for top tfs
#########################
rule rabit_rerun:
    input:
        dynamic(join(config['output'], 'rabit_rerun/output/{gene}.{ctype2}.txt'))
    params:
        output_dir=config['output']
    output:
        join(config['output'], 'rabit_formatted_results.txt')
    shell:
        "python scripts/format_rabit_results.py "
        "   -i {params.output_dir}/rabit_rerun/output "
        "   -b {params.output_dir}/rabit/output "
        "   -o {output} "
