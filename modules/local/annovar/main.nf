process ANNOVAR {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(vcf), path(tbi)
    path humandb

    output:
    tuple val(meta), path("*_multianno.vcf"), emit: vcf          
    tuple val(meta), path("*_multianno.txt"), emit: annotations 
    tuple val(meta), path("*_final_multianno_report.tsv"), emit: tsv, optional: true  
    path "versions.yml"                     , emit: versions
    
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    if [[ "$vcf" == *.gz ]]; then
        bgzip -d -c $vcf > local_unzipped_input.vcf
    else
        cp $vcf local_unzipped_input.vcf
    fi

    _JAVA_OPTIONS="-Xlog:perf+memops=off -XX:+PerfDisableSharedMem -XX:-UsePerfData" \\
    JAVA_TOOL_OPTIONS="-Xlog:perf+memops=off -XX:+PerfDisableSharedMem -XX:-UsePerfData" \\
    table_annovar.pl \\
        local_unzipped_input.vcf \\
        $humandb \\
        -out $prefix \\
        -vcfinput \\
        $args

    # =====================================================================
    # DEPENDENCY-FREE HEADER RE-MAPPING, SNPEFF & VAF INJECT (PYTHON 3)
    # =====================================================================
    python3 - << 'ANNOVAR_PY'
    import csv

    txt_input_path = "${prefix}.hg38_multianno.txt"
    tsv_output_path = "${prefix}_final_multianno_report.tsv"

    clean_headers = {
        "Chr": "Chrom",
        "Func.refGene": "Gene_Region",
        "Gene.refGene": "Gene_Symbol",
        "GeneDetail.refGene": "Transcript_Detail",
        "ExonicFunc.refGene": "Exonic_Consequence",
        "AAChange.refGene": "Amino_Acid_Change",
        "avsnp150": "dbSNP_ID",
        "clinvar_20250715": "ClinVar_Significance",
        "gnomad41_exome": "gnomAD_Exome_Freq",
        "dbnsfp42c": "Functional_Predictions",
        "intervar_20250721": "InterVar_Automated"
    }

    standard_vcf_labels = ["VCF_ID", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE_DATA"]

    try:
        with open(txt_input_path, mode='r', newline='', encoding='utf-8') as infile, \
             open(tsv_output_path, mode='w', newline='', encoding='utf-8') as outfile:
            
            reader = csv.reader(infile, delimiter='\t')
            writer = csv.writer(outfile, delimiter='\t')
            
            raw_headers = next(reader)
            
            # Anchor positions using the absolute end of the raw row layout
            sample_data_idx = len(raw_headers) - 1
            format_string_idx = len(raw_headers) - 2

            new_headers = []
            other_info_count = 0
            
            for header in raw_headers:
                if "Otherinfo" in header:
                    if other_info_count < len(standard_vcf_labels):
                        new_headers.append(standard_vcf_labels[other_info_count])
                    else:
                        new_headers.append(f"SAMPLE_DATA_{other_info_count - 4}")
                    other_info_count += 1
                elif header in clean_headers:
                    new_headers.append(clean_headers[header])
                else:
                    new_headers.append(header)
            
            try:
                target_index = new_headers.index("CLNALLELEID")
            except ValueError:
                target_index = len(new_headers)
            
            vaf_headers = ["GT", "AD_Ref", "AD_Alt1", "AD_Alt2", "DP", "GQ", "PL", "VAF1", "VAF2"]
            snpeff_headers = ["SnpEff_Variant_Effect", "SnpEff_Functional_Impact"]
            
            for h in reversed(vaf_headers):
                new_headers.insert(target_index, h)
            for h in reversed(snpeff_headers):
                new_headers.insert(target_index, h)
                
            writer.writerow(new_headers)
            
            for row in reader:
                cleaned_row = [col if col not in ["", "NA", "-"] else "." for col in row]
                
                # --- STEP 1: DYNAMIC KEY-VALUE VAF EXTRACTION ---
                gt, ad_ref, ad_alt1, ad_alt2, dp_val, gq, pl, vaf1, vaf2 = [".", ".", ".", ".", ".", ".", ".", ".", "."]
                
                if sample_data_idx < len(cleaned_row) and format_string_idx < len(cleaned_row):
                    format_field = cleaned_row[format_string_idx]
                    sample_field = cleaned_row[sample_data_idx]
                    
                    if format_field != "." and sample_field != "." and ":" in sample_field:
                        fmt_keys = format_field.split(":")
                        smpl_vals = sample_field.split(":")
                        
                        # Dynamically map the keys to values (e.g., maps DP to its real value automatically)
                        if len(fmt_keys) == len(smpl_vals):
                            vcf_data = dict(zip(fmt_keys, smpl_vals))
                            
                            gt = vcf_data.get("GT", ".")
                            gq = vcf_data.get("GQ", ".")
                            pl = vcf_data.get("PL", ".")
                            dp_val = vcf_data.get("DP", ".")
                            
                            ad_raw = vcf_data.get("AD", ".")
                            if ad_raw != "." and "," in ad_raw:
                                ad_parts = ad_raw.split(",")
                                ad_ref = ad_parts[0]
                                ad_alt1 = ad_parts[1]
                                if len(ad_parts) > 2:
                                    ad_alt2 = ad_parts[2]
                            
                            # Calculate VAF using the dynamically discovered fields
                            try:
                                total_dp = float(dp_val)
                                if total_dp > 0:
                                    if ad_alt1 != ".":
                                        vaf1 = f"{round((float(ad_alt1) / total_dp) * 100, 2)}"
                                    if ad_alt2 != ".":
                                        vaf2 = f"{round((float(ad_alt2) / total_dp) * 100, 2)}"
                            except ValueError:
                                pass

                # --- STEP 2: EXTRACT SNPEFF TEXT MUTATION ---
                snpeff_effect = "."
                snpeff_impact = "."
                row_string = "|".join(cleaned_row)
                
                if "HIGH" in row_string:
                    snpeff_impact = "HIGH"
                elif "MODERATE" in row_string:
                    snpeff_impact = "MODERATE"
                elif "LOW" in row_string:
                    snpeff_impact = "LOW"
                elif "MODIFIER" in row_string:
                    snpeff_impact = "MODIFIER"
                    
                consequences = [
                    "synonymous_variant", "missense_variant", "stop_gained", 
                    "stop_lost", "frameshift_variant", "splice_donor_variant", 
                    "splice_acceptor_variant", "intron_variant", "upstream_gene_variant"
                ]
                for cons in consequences:
                    if cons in row_string:
                        snpeff_effect = cons
                        break
                
                # --- STEP 3: ARMED ARRAY INJECTION ---
                vaf_vals = [gt, ad_ref, ad_alt1, ad_alt2, dp_val, gq, pl, vaf1, vaf2]
                snpeff_vals = [snpeff_effect, snpeff_impact]
                
                for val in reversed(vaf_vals):
                    cleaned_row.insert(target_index, val)
                for val in reversed(snpeff_vals):
                    cleaned_row.insert(target_index, val)
                    
                writer.writerow(cleaned_row)
                
        print("Successfully generated final TSV matrix with fully reactive dynamic calculations.")
    except Exception as e:
        print(f"Header conversion failed: {e}")

    ANNOVAR_PY
    # =====================================================================

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        annovar: \$(echo "custom_local_install")
    END_VERSIONS
    """
}