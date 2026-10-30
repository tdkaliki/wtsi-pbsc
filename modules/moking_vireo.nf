process MOKING_VIREO {
    label 'cellsnp'
    container "/software/hgi/containers/yascp/yascp.cog.sanger.ac.uk-public-yascp_qc_jan_2025.sif"

    input:
        tuple val(sample_id), path(bam), path(bam_bai), path(barcodes), path(input_vcf)
    output:
        tuple val(sample_id), path("GT_donors.vireo.vcf.gz"), emit: sample_donor_vcf//, path("vireo__${sample_id}/GT_donors.vireo.vcf.gz.csi"), emit: sample_donor_vcf
        tuple val(sample_id), path("donor_ids.tsv"), emit: sample_donor_ids

    script:
    """
        echo "donor0" > new_name.txt
        #mkdir vireo__${sample_id}
        awk 'BEGIN {
            OFS = "\t"
            print "cell", "donor_id", "prob_max", "prob_doublet", "n_vars", "best_singlet", "best_doublet", "doublet_logLikRatio"
        }
        NF { print \$1, "donor0", "", "", "", "", "", "" }' ${barcodes} > donor_ids.tsv

        bcftools reheader -s new_name.txt ${input_vcf} -o GT_donors.vireo.vcf.gz
    """
}
