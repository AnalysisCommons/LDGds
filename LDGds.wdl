task calculate_LD {
	File gds_file
	File? sample_ids_file
	String? ref_var
	File? rsid_file
	String? interval
	Int? half_interval
	Int? min_mac
	Int? max_mac
	Float? min_maf
	Float? max_maf
	String? ld_method
	String out_pref
	String? visualization
	Int? memory
	Int disk

	

	command {
		R --vanilla --args ${gds_file} ${default="NA" sample_ids_file} ${default="NA" ref_var} ${default="NA" rsid_file} ${default="NA" interval} ${default="25000" half_interval} ${default="0" min_mac} ${default="10000000" max_mac} ${default="0.05" min_maf} ${default="1" max_maf} ${default="r" ld_method} ${out_pref} ${default="F" visualization} < /LDGds/calculate_LD.R
	}

	runtime {
		docker: "tmajarian/ldgds:v0.1"
		disks: "local-disk " + disk + " HDD"
		memory: select_first([memory,"5"]) + " GB"
	}

	output {
		File out_file = select_first(glob("${out_pref}*.csv"))
		File out_visualization = select_first(glob("${out_pref}*.png"))
	}
}

workflow LD_wf {
	File this_gds_file
	File? this_sample_ids_file
	String? this_ref_var
	File? this_rsid_file
	String? this_interval
	Int? this_half_interval
	Int? this_min_mac
	Int? this_max_mac
	Float? this_min_maf
	Float? this_max_maf
	String? this_ld_method
	String this_out_pref
	String? this_visualization

	Int? this_memory
	Int? this_disk

	Int this_def_disk = select_first([this_disk, ceil(size(this_gds_file, "GB")) + 20])

	if (defined(this_ref_var) || defined(this_interval) || defined(this_rsid_file)) {
		call calculate_LD {
			input: 
				gds_file = this_gds_file,
				sample_ids_file = this_sample_ids_file,
				ref_var = this_ref_var,
				rsid_file = this_rsid_file,
				interval = this_interval,
				half_interval = this_half_interval,
				min_mac = this_min_mac,
				max_mac = this_max_mac,
				min_maf = this_min_maf,
				max_maf = this_max_maf,
				ld_method = this_ld_method,
				out_pref = this_out_pref,
				visualization = this_visualization,
				memory = this_memory,
				disk = this_def_disk
		}
		String? pass_message = "Workflow completed"
	}

	if (!defined(this_ref_var) && !defined(this_interval) && !defined(this_rsid_file)) {
		String? fail_message = "Reference variant, interval, and rsid file were unspecified, no computation was initiated."
	}

	String? workflow_message_var = select_first([pass_message, fail_message])

	output {
        File? ld_file = calculate_LD.out_file
        File? ld_plot = calculate_LD.out_visualization
        String? workflow_message = workflow_message_var
    }

    parameter_meta {
	    this_gds_file: "[file, *.gds] GDS file of genotypes per sample."
		this_sample_ids_file: "[file, default = all samples] File of sample IDs desired for LD calculation. This file should contain one sample ID per line with no header."
		this_ref_var: "[string, chr:pos] Genetic variant for which LD should be calculated. If provided, output is a row vector with pairwise LD with this variant in each row entry. Variant format should be 'chromosome:position'. Any punctuation seperator may be used. Only the first two values separated by punctuation will be considered."
		this_rsid_file: "[file] A file with a list of rsids, one per line to calculate LD from."
		this_interval: "[string, chr:start:end] Genomic interval for whcih LD should be calculated. If provided, LD will be calculated for only those variants falling within this interval. Interval format should be 'chromosome:start:end'. Any punctuation seperators may be used and need not match. Only the first three values separated by punctuation will be considered."
		this_half_interval: "[int, default = 25kb] 1/2 of desired interval length if no interval is provided. When only a reference variant is provided, this value will be added and subtracted from the reference variant position to define the interval end and start, respectively. "
		this_min_mac: "[int, default = 0] Minimum minor allele count for variant to be included in LD calculation."
		this_max_mac: "[int, default = inf] Maximum minor allele count for variant to be included in LD calculation."
		this_min_maf: "[int, default = 5%] Minimum minor allele frequency for variant to be included in LD calculation."
		this_max_maf: "[int, default = 1] Maximum minor allele frequency for variant to be included in LD calculation."
		this_ld_method: "[string, default = 'r'] LD calculation method. This value refers to the output LD values. Refer to documentation for the SNPRelate package, specifically the function snpgdsLDMat. Possible values are: composite, correlation, r, and dprime with reasonable abbreviations accepted."
		this_out_pref: "[string] Prefix for output file."
		this_memory: "[int, default = 5GB] Amount of memory to request for computation in GB."
		this_disk: "[int, default = size(this_gds_file) + 20 GB] Amount of disk space to request for computation in GB."
		  
	}

    meta {
        author: "Tim Majarian"
        email: "tmajaria@broadinstitute.org"
        description: "Generate LD measures from genotypes in GDS format. This workflow will return LD information for a set of defined samples over a set of variants or a defined variant range. A flat file of LD values and a simple visualization are returned."
    }
}
