#################################
#       nCovPopDyn pipeline     #
#################################

import os


report: "report/workflow.rst"

basefilename = os.path.basename(config["samples"])[:-6]

report: "report/workflow.rst"

if config["R0"]:
	rule all:
		input:
			"results/r0/r0.csv"
else:
	rule all:
		input:
			"results/interpolation/interpolation.csv"


rule strip_whitespaces:
	input:
		config["samples"]
	output:
		expand("results/raw/{sample}_fixed1.fasta", sample=basefilename)
	conda:
		"env/env.yml"
	log:
		expand("logs/ws_{sample}.log", sample=basefilename)
	shell:
		"reformat.sh in={input} out={output} underscore ignorejunk overwrite=true 2> {log}"


rule samtools_faidx:
	input:
		expand("{reference}", reference=config["consensus"])
	output:
		expand("{reference}.fai", reference=config["consensus"])
	shell:
		"samtools faidx {input}"

rule replace_dashes:
	input:
		"results/raw/{sample}_fixed1.fasta"
	output:
		"results/raw/{sample}_fixed12.fasta"
	conda:
		"env/env.yml"
	log:
		"logs/dashes_{sample}.log"
	shell:
		"seqkit replace -s -p '-' -r 'N' {input} -o {output} 2> {log}"

rule samtools_dict:
	input:
		expand("{reference}", reference=config["consensus"])
	output:
		expand("{reference}.dict", reference=config["consensus"][:-5])
	shell:
		"samtools dict {input} -o {output}"

rule minimap_index_ref:
	input:
		expand("{reference}", reference=config["consensus"])
	output:
		expand("{reference}.mmi", reference=config["consensus"][:-5])
	conda:
		"env/env.yml"
	shell:
		"minimap2 -d {output} {input}"

rule minimap:
	input:
		ref = expand("{reference}.mmi", reference=config["consensus"][:-5]),
		s = "results/raw/{sample}_fixed12.fasta"
	output:
		"results/bam/{sample}.bam"
	conda:
		"env/env.yml"
	log:
		"logs/map_{sample}.log"
	shell:
		"minimap2 -a --eqx {input.ref} {input.s} | samtools view -Sb -F 0x900 > {output} 2> {log}"

rule sort_bam:
	input:
		"results/bam/{sample}.bam"
	output:
		"results/bam/{sample}-sorted.bam"
	log:
		"logs/sort_{sample}.log"
	conda:
		"env/env.yml"
	shell:
		"samtools sort {input} > {output} 2> {log}"

rule index_bam:
	input:
		"results/bam/{sample}-sorted.bam"
	output:
		"results/bam/{sample}-sorted.bam.bai"
	conda:
		"env/env.yml"
	log:
		"logs/index_{sample}.log"
	shell:
		"samtools index {input} 2> {log}"

rule samtools_idxstats:
	input:
		"results/bam/{sample}-sorted.bam"
	output:
		"results/bam/{sample}-sorted-stats.txt"
	conda:
		"env/env.yml"
	log:
		"logs/idxstats_{sample}.log"
	shell:
		"samtools idxstats {input} > {output}"

rule run_binning:
	input:
		bam = expand("results/bam/{sample}-sorted.bam", sample=basefilename),
		bai = expand("results/bam/{sample}-sorted.bam.bai", sample=basefilename),
		stats = expand("results/bam/{sample}-sorted-stats.txt", sample=basefilename)
	output:
		files_list = "results/bins/list_of_binnings.tsv",
		meta = "results/meta/meta_dates.tsv"
	params:
		eq_num = config["number_per_bin"],
		eq_days = config["days_per_bin"],
		bin_dir = "results/bins",
		reference = config["consensus"]
	conda:
		"env/env.yml"
	script:
		"scripts/binning/run_binning_count.py"

rule theta_estimates:
	input:
		"results/bins/list_of_binnings.tsv"
	output:
		"results/bins_incidence/table_merged_phi_estimates_var_from_size.tsv"
	params:
		ref = config["consensus"],
		cutoff = config["freq_cutoff"],
		min_bin_size = config["min_bin_size"],
		min_days_span = config["min_days_span"],
		max_days_span = config["max_days_span"],
		group = config["group"]
	conda:
		"env/env.yml"
	script:
		"scripts/metrics/run_fp.py"

rule interpolation:
	input:
		infile = "results/bins_incidence/table_merged_phi_estimates_var_from_size.tsv"
	params:
		rep_cases = config["reported_cases"],
		group = config["group"]
	conda:
		"env/env.yml"
	output:
		result = "results/incidence/incidence.csv",
		#abs_path = os.path.join(workflow.basedir,"results/splines/out_spline.pdf"),
	shell:
		"Rscript {workflow.basedir}/scripts/Rscripts/splines/computeInterpolation.R {input.infile} \
		{params.group} {output.result} {params.rep_cases}"

if config["R0"]:
	rule r0:
		input:
			infile = "results/incidence/incidence.csv"
		output:
			outfile = "results/r0/r0.csv"
		conda:
			"env/env.yml"
		shell:
			"Rscript {workflow.basedir}/scripts/Rscripts/splines/computeR0.R {input.infile} \
			{output.outfile}"
