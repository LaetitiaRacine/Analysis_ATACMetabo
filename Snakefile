from snakemake.io import expand, glob_wildcards

LINK_TAB = {}
SAMPLE = set()
CONDITION_TIME = set()

import csv
with open('link_tab.csv', mode='r') as file:
	reader = csv.reader(file)
	next(reader, None)  # skip the headers
	for row in reader:
		use, condition, time, donor, manip, raw_file = row
		if use == "false":
			continue

		LINK_TAB["_".join((condition, time, donor))] = raw_file
		SAMPLE.add("_".join((condition, time, donor)))
		CONDITION_TIME.add("_".join((condition, time)))

print(SAMPLE)
print(CONDITION_TIME)

wildcard_constraints:
	donor="D\d+",
	donors="D\d+(,D\d+)+",
	sample="\w+_\d{2}h_D\d+"

rule all :
	input :
		expand("D_results/downsampled_bam/{sample}_downsampled.bam", sample = SAMPLE),
		expand("D_results/downsampled_bam/{sample}_downsampled.bam.bai", sample = SAMPLE),
		"D_results/reports/nbreads_report.csv",
		"D_results/reports/qc_report.csv",
		# expand("D_results/macs2_output/{sample}_peaks.broadPeak", sample = SAMPLE),
		expand("D_results/genomic_ranges/static_peaks/{sample}.gr.rds", sample = SAMPLE),
		expand("D_results/genomic_ranges/static_peaks/{condition_time}_D1,D2,D3.gr.rds", condition_time = CONDITION_TIME),
		expand("D_results/readCount_matrix/static_peaks/featurecounts_{condition_time}_D1,D2,D3.txt", condition_time = CONDITION_TIME)

rule link_rename_raw :
	input : lambda wildcards : "A_raw_data/bam_files/" + LINK_TAB[wildcards.sample]
	output : "D_results/bam/{sample}.bam"
	shell : """ ln -s "$(pwd)/{input}" {output} """

# ====================
# Generic BAM analysis
# ====================

rule bam_qc :
	input : "{file}.bam"
	output : "{file}.bam.qc"
	conda : "B_environments/ATACMetabo_main_env.yaml"
	shell : """ samtools flagstat -O tsv {input} > {output} """

rule bam_nbreads :
	input : "{file}.bam"
	output : "{file}.bam.nbreads"
	conda : "B_environments/ATACMetabo_main_env.yaml"
	shell : """ samtools view -c {input} > {output} """

rule bam_indexing :
	input : "{file}.bam"
	output : "{file}.bam.bai"
	conda : "B_environments/ATACMetabo_main_env.yaml"
	shell : """ samtools index {input} > {output} """


# ================
# BAM downsampling
# ================

# Compute the downsampling value to equalize the number of reads all bam files
# params.output is used to persist if the downsampling value does not change
#   It is updated only if the value changed
#   So the downsampling is redone only if the value changed
checkpoint bam_downsampling_value :
	input : expand("D_results/bam/{sample}.bam.nbreads", sample = SAMPLE)
	output : temp("D_results/downsampled_bam/downsampling_value_new")
	params : output = "D_results/downsampled_bam/downsampling_value"
	shell : """
        min=$(cat {input} | sort -n | head -1)
        echo "Downsampling value = $min"
        if [ -f {params.output} ] && [ "$(cat {params.output})" = "$min" ]; then
          echo "Unchanged"
        else
          echo "Changed"
          echo "$min" > {params.output}
        fi
        echo "$min" > {output}
        """

def bam_downsampling_input(wildcards):
	checkpoints.bam_downsampling_value.get()
	return {
		'bam': "D_results/bam/{sample}.bam",
		'nbreads': "D_results/bam/{sample}.bam.nbreads",
		'downsampling_value': rules.bam_downsampling_value.params.output
	}

# Downsample bam file to the downsampling value in the downsampling_value file
rule bam_downsampling :
	input : unpack(bam_downsampling_input)
	output : "D_results/downsampled_bam/{sample}_downsampled.bam"
	conda : "B_environments/ATACMetabo_main_env.yaml"
	shell : """
        count=$(cat {input.nbreads})
        min=$(cat {input.downsampling_value})
        echo "Downsampling value = $min"
        downsampling_ratio=$(echo "$min/$count" | bc -l)
        samtools view -s 1$downsampling_ratio -b {input.bam} > {output}
        """


# =============
# Peak analysis
# =============

rule peak_calling :
	input : rules.bam_downsampling.output
	output : "D_results/macs2_output/{sample}_peaks.broadPeak"
	params :
		prefix = "{sample}",
		macs2_output_dir = "D_results/macs2_output"
	conda : "B_environments/ATACMetabo_main_env.yaml"
	shell : """ macs2 callpeak -t {input} -n {params.prefix} --outdir {params.macs2_output_dir} -f BAMPE -g hs -B --broad --broad-cutoff 0.1 """

rule broadPeak_to_grange :
	input : rules.peak_calling.output
	output : "D_results/genomic_ranges/static_peaks/{sample}.gr.rds"
	conda : "B_environments/ATACMetabo_main_env.yaml"
	shell : """ Rscript C_scripts/GRanges.R from_broadPeak -o {output} {input} """

rule static_peaks_intersection :
	input : lambda wildcards : expand("D_results/genomic_ranges/static_peaks/{{condition}}_{{time}}_{donor}.gr.rds", donor=wildcards.donors.split(","))
	output : "D_results/genomic_ranges/static_peaks/{condition}_{time}_{donors}.gr.rds"
	conda : "B_environments/ATACMetabo_main_env.yaml"
	shell : """ Rscript C_scripts/GRanges.R intersect -o {output} {input} """

rule static_peaks_featureCounts :
	input :
		bam = lambda wildcards : expand("D_results/downsampled_bam/{{condition}}_{{time}}_{donor}_downsampled.bam", donor=wildcards.donors.split(",")),
		intersect = rules.static_peaks_intersection.output
	output :
		readcount = "D_results/readCount_matrix/static_peaks/readcount_{condition}_{time}_{donors}.rds",
		featurecounts = "D_results/readCount_matrix/static_peaks/featurecounts_{condition}_{time}_{donors}.txt"
	conda : "B_environments/ATACMetabo_main_env.yaml"
	shell : """
        Rscript C_scripts/peaks_featureCounts.R \\
            --output_rds {output.readcount} \\
            --output_txt {output.featurecounts} \\
            {input.intersect} {input.bam}
        """

rule differential_peaks_union:
	input:
		"D_results/genomic_ranges/static_peaks/{condition_time_1}.gr.rds",
		"D_results/genomic_ranges/static_peaks/{condition_time_2}.gr.rds"
	output: "D_results/genomic_ranges/differential_peaks/{condition_time_1}_vs_{condition_time_2}.gr.rds"
	conda : "B_environments/ATACMetabo_main_env.yaml"
	shell : """ Rscript C_scripts/GRanges.R union -o {output} {input} """

# =======
# Reports
# =======

# Report on quality of original bam files
rule qc_report :
	input : expand("D_results/bam/{sample}.bam.qc", sample = SAMPLE)
	output : "D_results/reports/qc_report.csv"
	params : sample = SAMPLE
	shell : """
        sample=({params.sample})
        qc_files=({input})
        for ((i=0;i<${{#sample[*]}};++i)); do
          echo -e "${{sample[i]}}\n$(cat ${{qc_files[i]}})"
        done >> {output}
        """

# Report on nbreads in all original bam files
rule nbreads_prereport :
	input : expand("D_results/bam/{sample}.bam.nbreads", sample = SAMPLE),
	output : "D_results/reports/nbreads_prereport.csv"
	params : sample = SAMPLE
	shell : """
        sample=({params.sample})
        nbreads=($(cat {input}))
        echo "SAMPLE;NBREADS BEFORE DOWNSAMPLING" > {output}
        for ((i=0;i<${{#sample[*]}};++i)); do
          echo "${{sample[i]}};${{nbreads[i]}}"
        done >> {output}
        """

# Report nbreads before and after nbreads for all bam files
rule nbreads_report :
	input :
		before_downsampling = expand("D_results/bam/{sample}.bam.nbreads", sample = SAMPLE),
		after_downsampling = expand("D_results/downsampled_bam/{sample}_downsampled.bam.nbreads", sample = SAMPLE)
	output : "D_results/reports/nbreads_report.csv"
	params : sample = SAMPLE
	shell : """
        sample=({params.sample})
        nbreads_before=($(cat {input.before_downsampling}))
        nbreads_after=($(cat {input.after_downsampling}))
        echo "SAMPLE;NBREADS BEFORE DOWNSAMPLING;NBREADS AFTER DOWNSAMPLING" > {output}
        for ((i=0;i<${{#sample[*]}};++i)); do
          echo "${{sample[i]}};${{nbreads_before[i]}};${{nbreads_after[i]}}"
        done >> {output}
        """
