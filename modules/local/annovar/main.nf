process ANNOVAR {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(vcf), path(tbi)
    path humandb

    output:
    tuple val(meta), path("*_multianno.vcf"), emit: vcf          
    tuple val(meta), path("*_multianno.txt"), emit: annotations 
    tuple val(meta), path("*_final_multianno_report.tsv"), emit: tsv, optional: true  // Registers clean report
    path "versions.yml"                     , emit: versions
    
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    # 1. Manually decompress the incoming VCF to prevent nested stream contamination
    if [[ "$vcf" == *.gz ]]; then
        bgzip -d -c $vcf > local_unzipped_input.vcf
    else
        cp $vcf local_unzipped_input.vcf
    fi

    # 2. Execute ANNOVAR while hard-coding JVM suppression flags to avoid /tmp lock errors
    _JAVA_OPTIONS="-Xlog:perf+memops=off -XX:+PerfDisableSharedMem -XX:-UsePerfData" \\
    JAVA_TOOL_OPTIONS="-Xlog:perf+memops=off -XX:+PerfDisableSharedMem -XX:-UsePerfData" \\
    table_annovar.pl \\
        local_unzipped_input.vcf \\
        $humandb \\
        -out $prefix \\
        -vcfinput \\
        $args

    # =====================================================================
    # DEPENDENCY-FREE HEADER RE-MAPPING INJECT (PURE PYTHON STANDARD LIB)
    # =====================================================================
    python3 - << 'ANNOVAR_PY'
    import csv

    txt_input_path = "${prefix}.hg38_multianno.txt"
    tsv_output_path = "${prefix}_final_multianno_report.tsv"

    # Define clean clinical annotation mappings
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
            
            # Extract and parse raw headers
            raw_headers = next(reader)
            new_headers = []
            other_info_count = 0
            
            for header in raw_headers:
                # 1. Map messy Otherinfo columns to standard VCF syntax labels
                if "Otherinfo" in header:
                    if other_info_count < len(standard_vcf_labels):
                        new_headers.append(standard_vcf_labels[other_info_count])
                    else:
                        new_headers.append(f"SAMPLE_DATA_{other_info_count - 4}")
                    other_info_count += 1
                # 2. Map standard database annotations to readable system properties
                elif header in clean_headers:
                    new_headers.append(clean_headers[header])
                else:
                    new_headers.append(header)
            
            writer.writerow(new_headers)
            
            # Process remaining variant lines, formatting blank entries to dots
            for row in reader:
                cleaned_row = [col if col not in ["", "NA", "-"] else "." for col in row]
                writer.writerow(cleaned_row)
                
        print("Successfully generated final TSV matrix.")
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