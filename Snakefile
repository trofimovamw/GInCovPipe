#################################
#       nCovPopDyn pipeline     #
#################################

import os

import numpy as np

configfile: "config.yaml"
working_dir = os.getcwd()

# List of bin files (BAM)
bin_ind = np.arange(1000)
bin_ind_list = bin_ind.tolist()
bin_ind_padded = []
for ind in bin_ind_list:
	inds = str(ind)
	inds0 = inds.zfill(4)
	bin_ind_padded.append(inds0)

# List of binning modes for wildcards
eq_size = "eq_size_names_"
binnings_names = []
for name in config["number_per_bin"]:
	s = eq_size+str(name)
	binnings_names.append(s)
eq_days = "eq_days_"
fuzzy_days = "fuzzy_days_"
for name in config["days_per_bin"]:
	s1 = eq_days+str(name)
	s2 = fuzzy_days+str(name)
	binnings_names.append(s1)
	binnings_names.append(s2)
binnings_names.append("cal_week")


rule all:
	input:
		"results/splines/out_spline.pdf"

if config["samples_meta"] != "":
	rule add_metadata:
		input:
			fasta = "raw/{sm}.fasta",
			meta = "raw/{sm}.tsv"
		output:
			"raw/{sm}_reform.fasta"
		conda:
			"env/env.yml"
		shell:
			"python scripts/add_date_from_metadata.py {input.fasta} {input.meta} {output}"

	rule strip_whitespaces2:
		input:
			"raw/{sm}_reform.fasta"
		output:
			"results/raw/{sm}_fixed1.fasta"
		conda:
			"env/env.yml"
		log:
			"logs/ws_{sm}.log"
		shell:
			"reformat.sh in={input} out={output} underscore ignorejunk overwrite=true 2> {log}"

if config["samples"] != "":
	rule strip_whitespaces1:
		input:
			expand("raw/{sample}.fasta", sample=config["samples"])
		output:
			expand("results/raw/{sample}_fixed1.fasta", sample=config["samples"])
		conda:
			"env/env.yml"
		log:
			expand("logs/ws_{sample}.log", sample=config["samples"])
		shell:
			"reformat.sh in={input} out={output} underscore ignorejunk overwrite=true 2> {log}"

rule bwa_index:
	input:
		expand("consensus/{reference}.fasta", reference=config["consensus"])
	output:
		expand("consensus/{reference}.fasta.ann", reference=config["consensus"])
	shell:
		"bwa index {input}"

rule samtools_dict:
	input:
		expand("consensus/{reference}.fasta", reference=config["consensus"])
	output:
		expand("consensus/{reference}.dict", reference=config["consensus"])
	shell:
		"/Users/mariatrofimova/Downloads/bin/samtools dict {input} -o {output}"

rule samtools_faidx:
	input:
		expand("consensus/{reference}.fasta", reference=config["consensus"])
	output:
		expand("consensus/{reference}.fasta.fai", reference=config["consensus"])
	shell:
		"samtools faidx {input}"

rule merge_fasta:
	input:
		fasta1 = expand("results/raw/{sm}_fixed1.fasta", sm=config["samples_meta"]),
		fasta2 = expand("results/raw/{sample}_fixed1.fasta", sample=config["samples"])
	output:
		expand("results/raw/{sample}{sm}_fixed1.fasta", sample=config["samples"], sm=config["samples_meta"])
	shell:
		"cat {input.fasta1} {input.fasta2} > {output}"

rule replace_dashes:
	input:
		"results/raw/{sample}{sm}_fixed1.fasta"
	output:
		"results/raw/{sample}{sm}_fixed12.fasta"
	conda:
		"env/env.yml"
	log:
		"logs/dashes_{sample}{sm}.log"
	shell:
		"seqkit replace -s -p '-' -r 'N' {input} -o {output} 2> {log}"

rule map:
	input:
		ref = expand("consensus/{reference}.fasta", reference=config["consensus"]),
		index = expand("consensus/{reference}.fasta.ann", reference=config["consensus"]),
		s = "results/raw/{sample}{sm}_fixed12.fasta"
	output:
		sam = "results/bam/{sample}{sm}.bam"
	conda:
		"env/env.yml"
	log:
		"logs/map_{sample}{sm}.log"
	shell:
		#"bwa mem -t 4  {input.ref} {input.s} > {output.sam} 2> {log}"
		"bwa mem -t 4 -O 20 -E 0 {input.ref} {input.s} | samtools view -Sb > {output} 2> {log}"

rule sort_bam:
	input:
		"results/bam/{sample}{sm}.bam"
	output:
		"results/bam/{sample}{sm}-sorted.bam"
	log:
		"logs/sort_{sample}{sm}.log"
	conda:
		"env/env.yml"
	shell:
		"/Users/mariatrofimova/Downloads/bin/samtools sort {input} > {output} 2> {log}"

rule index_bam:
	input:
		"results/bam/{sample}{sm}-sorted.bam"
	output:
		"results/bam/{sample}{sm}-sorted.bam.bai"
	conda:
		"env/env.yml"
	log:
		"logs/index_{sample}{sm}.log"
	shell:
		"samtools index {input} 2> {log}"

rule run_binning:
	input:
		bam = expand("results/bam/{sample}{sm}-sorted.bam", sample=config["samples"], sm=config["samples_meta"]),
		bai = expand("results/bam/{sample}{sm}-sorted.bam.bai", sample=config["samples"], sm=config["samples_meta"])
	output:
		files_list = "results/bins/list_of_binnings.tsv"
	params:
		eq_num = config["number_per_bin"],
		eq_days = config["days_per_bin"],
		bin_dir = "results/bins",
		reference = config["consensus"]
	conda:
		"env/env.yml"
	script:
		"scripts/binning/run_binning_count.py"

rule fix_cigars_subprocess:
	input:
		binnings = "results/bins/list_of_binnings.tsv",
		dict = expand("consensus/{reference}.dict", reference=config["consensus"]),
		fai = expand("consensus/{reference}.fasta.fai", reference=config["consensus"])
	params:
		java_tool = config["samfixcigars"],
		samtools = config["samtools"],
		ref = config["consensus"]
	output:
		"results/fixed_cigars_bins/list_of_binnings.tsv"
	script:
		"scripts/binning/fix_cigars_subprocess.py"

rule theta_estimates:
	input:
		"results/fixed_cigars_bins/list_of_binnings.tsv"
	output:
		"results/plots/table_merged_thetas_var_from_size.tsv"
	params:
		ref = config["consensus"],
		rep_cases = config["reported_cases"],
		min_bin_size = config["min_bin_size"],
		smoothing = config["smoothing"],
		log_transform = config["log_transform"],
		min_days_span = config["min_days_span"],
		meta = config["samples"]
	script:
		"scripts/metrics/run_fp.py"

rule splines:
	input:
		infile = "results/plots/table_merged_thetas_var_from_size.tsv",
		meta_abs_path = os.path.join(workflow.basedir,"raw/%s.tsv" % config["samples"])
	params:
		date_m = config["date_m"]
	output:
		result = "results/splines/out_spline.pdf",
		abs_path = os.path.join(workflow.basedir,"results/splines/out_spline.pdf")
	shell:
		"Rscript scripts/Rscripts/splines/computeSpline.R {input.infile} {output.abs_path} {input.meta_abs_path} {params.date_m} trueN"
