process ANNOVAR {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(vcf), path(tbi)
    path humandb

    output:
    tuple val(meta), path("*_multianno.vcf"), emit: vcf          
    tuple val(meta), path("*_multianno.txt"), emit: annotations 
    path "versions.yml"                     , emit: versions
    
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    table_annovar.pl \\
        $vcf \\
        $humandb \\
        -out $prefix \\
        -vcfinput \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        annovar: \$(echo "custom_local_install")
    END_VERSIONS
    """
}