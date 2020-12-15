from snakemake.io import expand, glob_wildcards
from collections import defaultdict
from functools import partial

LINK_TAB = {}
SAMPLE = defaultdict(partial(defaultdict, dict))

import csv
with open('link_tab.csv', mode='r') as file:
	reader = csv.reader(file)
	next(reader, None)  # skip the headers
	for row in reader:
		use, condition, time, donor, manip, raw_file = row
		if use == "false":
			continue

		#LINK_TAB["_".join((condition, time, donor))] = raw_file
		#SAMPLE.add("_".join((condition, time, donor)))

		SAMPLE[condition][time][donor] = raw_file

def list_sample():
	result = []
	for condition, times in SAMPLE.items():
		for time, donors in times.items():
			for donor, file in donors.items():
				result.append("_".join((condition, time, donor)))
	return result

# Fonction qui liste les condition_time disponibles
# Arguments optionnels, servent à restreindre la liste
def list_condition_time(conditions_name=None, times_name=None):
	result = []

	conditions = SAMPLE.items()
	if conditions_name is not None:
		conditions = filter(lambda item : item[0] in conditions_name, conditions)
	for condition, times in conditions:
		times = times.items()
		if times_name is not None:
			times = filter(lambda item : item[0] in times_name, times)
		for time, donors in times:
			result.append(condition + "_" + time + "_" + ",".join(donors.keys()))
	return result

def list_unions():
	result = []

	# Comparaison de Xvivo avec tous les premiers points de temps des autres conditions
	xvivo_condition_time = list_condition_time(["Xvivo"], ["00h"])[0]
	for condition in SAMPLE.keys():
		if condition == "Xvivo":
			continue
		condition_initial_time = sorted(list_condition_time([condition]))[0]
		result.append(xvivo_condition_time + "_vs_" + condition_initial_time)

	# Comparaison à un point de temps donné de MP aux autres conditions
	conditions_but_mp = list(filter(lambda item: item != "MP", SAMPLE.keys()))
	for mp_condition_time in list_condition_time(["MP"]):
		time = mp_condition_time.split("_")[1]
		for condition_time in list_condition_time(conditions_but_mp, [time]):
			result.append(mp_condition_time + "_vs_" + condition_time)

	# Comparaison entre poitns de temps au sein d'une condition
	for condition in SAMPLE.keys():
		condition_times = sorted(list_condition_time([condition]))
		for i in range(len(condition_times) - 1):
			result.append(condition_times[i] + "_vs_" + 	condition_times[i+1])

	return result

# print(list_condition_time())
# print(list_sample())
# print(list_unions())

wildcard_constraints:
	donor="D\d+",
	donors="D\d+(,D\d+)+",
	sample="\w+_\d{2}h_D\d+"


rule all :
	input :
		### Test bloc 2
		expand("D_results/downsampled_bam/{sample}_downsampled.bam.bai", sample = list_sample()),
		expand("D_results/genomic_ranges/static_peaks/{sample}.gr.rds", sample = list_sample()),
		"D_results/reports/qc_report.csv",
		"D_results/reports/nbreads_report.csv",
		"D_results/reports/nbpeaks_report.csv",
		"D_results/reports/plot_hist_donor_qc_report:QCP_nb_read.png",
		"D_results/reports/plot_hist_donor_nbreads_report:nbreads_before_downsampling.png",
		"D_results/reports/plot_hist_manip_nbreads_report:nbreads_before_downsampling.png",
		"D_results/reports/plot_hist_donor_nbreads_report:nbreads_after_downsampling.png",
		"D_results/reports/plot_hist_manip_nbreads_report:nbreads_after_downsampling.png",
		"D_results/reports/plot_hist_donor_nbpeaks_report:nbpeaks.png",
		"D_results/reports/plot_hist_manip_nbpeaks_report:nbpeaks.png",
		"D_results/reports/plot_line_cond_nbpeaks_report:nbpeaks.png"
		### Test bloc 3
		# expand("D_results/readCount_matrix/static_peaks/featurecounts_{condition_time}.txt", condition_time = list_condition_time()),
		# expand("D_results/readCount_matrix/differential_peaks/featurecounts_{union}.txt", union = list_unions())
		# expand("D_results/genomic_ranges/static_peaks/{condition_time}.gr.rds", condition_time = list_condition_time()),

#*******************************************************************************************************************************************************
#*** Bloc 1 : Alignment of FASTQ => BAM
#*******************************************************************************************************************************************************

# cf script Sophie


#*******************************************************************************************************************************************************
#*** Bloc 2 : Individual study and QC of each sample
#*******************************************************************************************************************************************************

rule link_rename_raw :
	input : lambda wildcards : "A_raw_data/bam_files/" + SAMPLE[wildcards.condition][wildcards.time][wildcards.donor]
	output : "D_results/bam/{condition}_{time}_{donor}.bam"
	shell : """ ln -s "$(pwd)/{input}" {output} """

# ====================
# Generic BAM analysis
# ====================

rule bam_qc :
	input : "{file}.bam"
	output : "{file}.bam.qc"
	conda : "B_environments/ATACMetabo_main_env.locked.yaml"
	shell : """ samtools flagstat -O tsv {input} > {output} """

rule bam_nbreads :
	input : "{file}.bam"
	output : "{file}.bam.nbreads"
	conda : "B_environments/ATACMetabo_main_env.locked.yaml"
	shell : """ samtools view -c {input} > {output} """

rule bam_indexing :
	input : "{file}.bam"
	output : "{file}.bam.bai"
	conda : "B_environments/ATACMetabo_main_env.locked.yaml"
	shell : """ samtools index {input} > {output} """


# ================
# BAM downsampling
# ================

# Compute the downsampling value to equalize the number of reads all bam files
# params.output is used to persist if the downsampling value does not change
#   It is updated only if the value changed
#   So the downsampling is redone only if the value changed
checkpoint bam_downsampling_value :
	input : expand("D_results/bam/{sample}.bam.nbreads", sample = list_sample())
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
	conda : "B_environments/ATACMetabo_main_env.locked.yaml"
	shell : """
        count=$(cat {input.nbreads})
        min=$(cat {input.downsampling_value})
        echo "Downsampling value = $min"
        downsampling_ratio=$(echo "$min/$count" | bc -l)
        samtools view -s 1$downsampling_ratio -b {input.bam} > {output}
        """

#=============
# Peak analysis
#=============

rule peak_calling :
	input : rules.bam_downsampling.output
	output : "D_results/macs2_output/{sample}_peaks.broadPeak"
	params :
		prefix = "{sample}",
		macs2_output_dir = "D_results/macs2_output"
	conda : "B_environments/ATACMetabo_main_env.locked.yaml"
	shell : """ macs2 callpeak -t {input} -n {params.prefix} --outdir {params.macs2_output_dir} -f BAMPE -g hs -B --broad --broad-cutoff 0.1 """

rule broadPeak_to_grange :
	input : rules.peak_calling.output
	output : "D_results/genomic_ranges/static_peaks/{sample}.gr.rds"
	conda : "B_environments/ATACMetabo_main_env.locked.yaml"
	shell : """ Rscript C_scripts/GRanges.R from_broadPeak -o {output} {input} """

# =======
# Reports and QC plots
# =======

# Report on quality of original bam files
rule qc_report :
	input : expand("D_results/bam/{sample}.bam.qc", sample = list_sample())
	output : "D_results/reports/qc_report.csv"
	params : sample = list_sample()
	shell : """
        sample=({params.sample})
        qc_files=({input})
        echo "condition,time,donor,QCP_nb_read,QCF_nb_read,QCP_secondary,QCF_secondary,QCP_supplementary,QCF_supplementary,QCP_duplicates,QCF_duplicates,QCP_nb_mapped,QCF_nb_mapped,QCP_pct_mapped,QCF_pct_mapped,QCP_paired_in_sequencing,QCF_paired_in_sequencing,QCP_read1,QCF_read1,QCP_read2,QCF_read2,QCP_nb_properly_paired,QCF_nb_properly_paired,QCP_pct_properly_paired,QCF_pct_properly_paired,QCP_with_itself_and_mate_mapped,QCF_with_itself_and_mate_mapped,QCP_nb_singletons,QCF_nb_singletons,QCP_pct_singletons,QCF_pct_singletons,QCP_with_mate_mapped_to_a_different_chr,QCF_with_mate_mapped_to_a_different_chr,QCP_with_mate_mapped_to_a_different_chr_mapQ_over_5,QCF_with_mate_mapped_to_a_different_chr_mapQ_over_5" > {output}
        for ((i=0;i<${{#sample[*]}};++i)); do
          cond=$(echo "${{sample[i]}}" | cut -d_ -f1)
          time=$(echo "${{sample[i]}}" | cut -d_ -f2)
          donor=$(echo "${{sample[i]}}" | cut -d_ -f3)
          qc_info=$(cat "${{qc_files[i]}}" | awk -F$'\t' '{{ print $1","$2 }}' | paste -sd ',' -)
          echo -e "${{cond}},${{time}},${{donor}},${{qc_info}}"
        done >> {output}
        """

# Report on nbreads in all original bam files
rule nbreads_prereport :
	input : expand("D_results/bam/{sample}.bam.nbreads", sample = list_sample()),
	output : "D_results/reports/nbreads_prereport.csv"
	params : sample = list_sample()
	shell : """
        sample=({params.sample})
        nbreads=($(cat {input}))
        echo "condition,time,donor,nbreads_before_downsampling" > {output}
        for ((i=0;i<${{#sample[*]}};++i)); do
          cond=$(echo "${{sample[i]}}" | cut -d_ -f1)
          time=$(echo "${{sample[i]}}" | cut -d_ -f2)
          donor=$(echo "${{sample[i]}}" | cut -d_ -f3)
          echo "$cond,$time,$donor,${{nbreads[i]}}"
        done >> {output}
        """

# Report on nbreads before and after downsampling for all bam files
rule nbreads_report :
	input :
		before_downsampling = expand("D_results/bam/{sample}.bam.nbreads", sample = list_sample()),
		after_downsampling = expand("D_results/downsampled_bam/{sample}_downsampled.bam.nbreads", sample = list_sample())
	output : "D_results/reports/nbreads_report.csv"
	params : sample = list_sample()
	shell : """
        sample=({params.sample})
        nbreads_before=($(cat {input.before_downsampling}))
        nbreads_after=($(cat {input.after_downsampling}))
        echo "condition,time,donor,nbreads_before_downsampling,nbreads_after_downsampling" > {output}
        for ((i=0;i<${{#sample[*]}};++i)); do
          cond=$(echo "${{sample[i]}}" | cut -d_ -f1)
          time=$(echo "${{sample[i]}}" | cut -d_ -f2)
          donor=$(echo "${{sample[i]}}" | cut -d_ -f3)
          echo "$cond,$time,$donor,${{nbreads_before[i]}},${{nbreads_after[i]}}"
        done >> {output}
        """

# PB : relance tout le pipeline avec les peaks calling... même avec ancient je ne comprends pas
# reason: Input files updated by another job: D_results/downsampled_bam/MP_06h_D2_downsampled.bam
# quand on utilise --touch : aucune règle n'est exécutée

# Report on nb_peaks in each sample following peak_calling rule
rule nbpeaks_report :
	input : expand("D_results/macs2_output/{sample}_peaks.broadPeak", sample = list_sample())
	output : "D_results/reports/nbpeaks_report.csv"
	params : sample = list_sample()
	shell : """
        sample=({params.sample})
        nbpeaks=($(wc -l {input} | awk '{{ print $1 }}'))
        echo "condition,time,donor,nbpeaks" > {output}
        for ((i=0;i<${{#sample[*]}};++i)); do
          cond=$(echo "${{sample[i]}}" | cut -d_ -f1)
          time=$(echo "${{sample[i]}}" | cut -d_ -f2)
          donor=$(echo "${{sample[i]}}" | cut -d_ -f3)
          echo "$cond,$time,$donor,${{nbpeaks[i]}}"
        done >> {output}
	"""

# Draw QC plots from nbreads_report.csv and nb_peaks_report.csv
rule plot_reports :
	wildcard_constraints:
		format="hist_donor|hist_manip|line_cond",
		colname="[a-zA-Z_]+",
		reportname="[a-z_]+"
	input :"D_results/reports/{reportname}.csv"
	output : "D_results/reports/plot_{format}_{reportname}:{colname}.png",
	conda : "B_environments/ATACMetabo_main_env.locked.yaml"
	shell : """ Rscript C_scripts/plot_QC_variation.R -o {output} {wildcards.format} {input} {wildcards.colname}
			"""

#*******************************************************************************************************************************************************
#*** Bloc 3 : Crossed study of samples
#*******************************************************************************************************************************************************

# rule static_peaks_intersection :
# 	input : lambda wildcards : expand("D_results/genomic_ranges/static_peaks/{{condition}}_{{time}}_{donor}.gr.rds", donor=wildcards.donors.split(","))
# 	output : "D_results/genomic_ranges/static_peaks/{condition}_{time}_{donors}.gr.rds"
# 	conda : "B_environments/ATACMetabo_main_env.locked.yaml"
# 	shell : """ Rscript C_scripts/GRanges.R intersect -o {output} {input} """
#
# rule static_peaks_featureCounts :
# 	input :
# 		bam = lambda wildcards : expand("D_results/downsampled_bam/{{condition}}_{{time}}_{donor}_downsampled.bam", donor=wildcards.donors.split(",")),
# 		intersect = rules.static_peaks_intersection.output
# 	output :
# 		readcount = "D_results/readCount_matrix/static_peaks/readcount_{condition}_{time}_{donors}.rds",
# 		featurecounts = "D_results/readCount_matrix/static_peaks/featurecounts_{condition}_{time}_{donors}.txt"
# 	conda : "B_environments/ATACMetabo_main_env.locked.yaml"
# 	shell : """
#         Rscript C_scripts/peaks_featureCounts.R \\
#             --output_rds {output.readcount} \\
#             --output_txt {output.featurecounts} \\
#             {input.intersect} {input.bam}
#         """

# rule differential_peaks_union:
# 	input:
# 		"D_results/genomic_ranges/static_peaks/{condition_time_1}_{donors_1}.gr.rds",
# 		"D_results/genomic_ranges/static_peaks/{condition_time_2}_{donors_2}.gr.rds"
# 	output: "D_results/genomic_ranges/differential_peaks/{condition_time_1}_{donors_1}_vs_{condition_time_2}_{donors_2}.gr.rds"
# 	conda : "B_environments/ATACMetabo_main_env.locked.yaml"
# 	shell : """ Rscript C_scripts/GRanges.R union -o {output} {input} """
#
# rule differential_peaks_featureCounts :
# 	input :
# 		bam1 = lambda wildcards : expand("D_results/downsampled_bam/{{condition_time_1}}_{donor}_downsampled.bam", donor=wildcards.donors_1.split(",")),
# 		bam2 = lambda wildcards : expand("D_results/downsampled_bam/{{condition_time_2}}_{donor}_downsampled.bam", donor=wildcards.donors_2.split(",")),
# 		union = rules.differential_peaks_union.output
# 	output :
# 		readcount = "D_results/readCount_matrix/differential_peaks/readcount_{condition_time_1}_{donors_1}_vs_{condition_time_2}_{donors_2}.rds",
# 		featurecounts = "D_results/readCount_matrix/differential_peaks/featurecounts_{condition_time_1}_{donors_1}_vs_{condition_time_2}_{donors_2}.txt"
# 	conda : "B_environments/ATACMetabo_main_env.locked.yaml"
# 	shell : """
#         Rscript C_scripts/peaks_featureCounts.R \\
#             --output_rds {output.readcount} \\
#             --output_txt {output.featurecounts} \\
#             {input.union} {input.bam1} {input.bam2}
#         """
#
# rule annotate_grange :
# 	input :
# 		annotations = "A_raw_data/Annotation_TSS_pm1kb_int_ex_53utr_ctcf_cpg_histo_gr.rda",
# 		grange = "D_results/genomic_ranges/{folder}/{file}.gr.rds"
# 	output :
# 		grange_annot = "D_results/genomic_ranges/{folder}_annotated/{file}_ann.gr.rds",
# 		csv_annot = "D_results/genomic_ranges/{folder}_annotated/{file}_ann.csv"
# 	conda : "B_environments/ATACMetabo_main_env.locked.yaml"
# 	shell : """
#         Rscript C_scripts/annotate_grange.R \\
#             --output_csv {output.csv_annot} \\
#             -a {input.annotations} \\
#             {input.grange} \\
#             {output.grange_annot}
#         """
