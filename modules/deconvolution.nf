

process MPILEUP {
    label 'mpileup'
    container "/software/hgi/containers/yascp/yascp.cog.sanger.ac.uk-public-yascp_qc_jan_2025.sif"
    publishDir "${params.results_output}deconvolution/mpileup", mode: 'copy',
        saveAs: { filename -> (filename.endsWith('.bam') || filename.endsWith('.bai')) ? null : filename }

    input:
        tuple val(sample_id), path(bam), path(bam_bai)
        path(ref_gen)
    output:
        tuple val(sample_id), path(bam), path(bam_bai), path("${sample_id}__barcodes.txt"), path("${sample_id}__piled_up_reads.vcf")
    script:
    """
        samtools view ${bam} | awk '{ for(i=12;i<=NF;i++) if(\$i ~ /^CB:Z:/) {split(\$i,a,":"); print a[3]} }' | sort | uniq > ${sample_id}__barcodes.txt
        # Step 1: Generate VCF of observed SNPs using bcftools mpileup and call
        bcftools mpileup \
        -f ${ref_gen} \
        -q 20 -Q 20 \
        -a AD,DP \
        -Ou ${bam} | \
        bcftools call -mv -Ov -o ${sample_id}__piled_up_reads.vcf
    """
}

process SUBSET_VCF {
    label 'subset_vcf'
    container "/software/hgi/containers/yascp/yascp.cog.sanger.ac.uk-public-yascp_qc_jan_2025.sif"
    publishDir "${params.results_output}deconvolution/mpileup", mode: 'copy',
        saveAs: { filename -> (filename.endsWith('.bam') || filename.endsWith('.bai')) ? null : filename }

    input:
        tuple val(sample_id), path(bam), path(bam_bai), path(sample_id__barcodes), path(sample_id__piled_up_reads)
        path(subset_regions_bed)
    output:
        tuple val(sample_id), path(bam), path(bam_bai), path(sample_id__barcodes), path("${sample_id}__piled_up_reads__subset.vcf.gz")
    script:
    """
        bgzip -c ${sample_id__piled_up_reads} > ${sample_id}__tmp.vcf.gz
        tabix -p vcf ${sample_id}__tmp.vcf.gz
        bcftools view -R ${subset_regions_bed} ${sample_id}__tmp.vcf.gz -Oz -o ${sample_id}__piled_up_reads__subset.vcf.gz
    """
}

process CELLSNP {
    label 'cellsnp'
    container "/software/hgi/containers/yascp/yascp.cog.sanger.ac.uk-public-yascp_qc_jan_2025.sif"
    publishDir "${params.results_output}deconvolution/cellsnp", mode: 'copy'

    input:
        tuple val(sample_id), path(bam), path(bam_bai), path(barcodes), path(piled_up_reads)
    output:
        tuple val(sample_id),path("cellsnp__${sample_id}")
    script:
    """
        # Step 2: Run CellSNP-lite to extract SNPs per cell
        cellsnp-lite \
        -s ${bam} \
        -O cellsnp__${sample_id} \
        -b ${barcodes} \
        -R ${piled_up_reads} \
        -p $task.cpus \
        --minMAF 0.01 \
        --minCOUNT 1 \
        --gzip \
        --cellTAG CB \
        --UMItag None
    """
}



process VIREO {
    label 'deconvolution'

    container "/software/hgi/containers/yascp/yascp.cog.sanger.ac.uk-public-yascp_qc_jan_2025.sif"

    publishDir "${params.results_output}deconvolution/vireo", mode: 'copy', pattern: "{vireo__*,barcodes__*.tsv}"
    

    input:
        tuple val(sample_id),path(cellsnp),val(nr_samples)
    output:
        tuple val(sample_id), path("barcodes__*.tsv"), emit: barcodes_tuple
        tuple val(sample_id), path("${sample_id}.bc_list.txt.gz"), emit: barcode_list
        tuple val(sample_id), path("vireo__${sample_id}"), emit: vireo_results
        tuple val(sample_id), path("vireo__${sample_id}/GT_donors.vireo.vcf.gz"), emit: sample_donor_vcf//, path("vireo__${sample_id}/GT_donors.vireo.vcf.gz.csi"), emit: sample_donor_vcf
        tuple val(sample_id), path("vireo__${sample_id}/donor_ids.tsv"), emit: sample_donor_ids

    script:
    """
        # Step 3: Run Vireo for donor deconvolution without genotypes
        vireo \
        -c ${cellsnp} \
        -N ${nr_samples} \
        -o vireo__${sample_id} \
        --randSeed 1 \
        --nInit 200 \
        -p $task.cpus

        #Split the donor barcodes in an independent files for next step of bam splits.
        awk 'NR > 1 && \$2 != "unassigned" && \$2 != "doublet" {print > ("barcodes__" \$2 ".tsv")}' vireo__${sample_id}/donor_ids.tsv
        cat barcodes__*.tsv | awk -F'\t' -v prefix="${sample_id}_" 'BEGIN {OFS="\t"} NF >= 2 {print \$1, prefix\$2}' > ${sample_id}.bc_list.txt
        gzip ${sample_id}.bc_list.txt
    """
}

process CHECK_VIREO {
    label 'deconvolution'

    container "/software/hgi/containers/yascp/yascp.cog.sanger.ac.uk-public-yascp_qc_jan_2025.sif"    

    input:
        tuple val(sample_id),path(vireo_dir)
    output:
        tuple val(sample_id), path("*_donors_check"), emit: check_res

    script:
    def args = []
    if (params.check_vireo_min_cells != null) args << "--min_cells ${params.check_vireo_min_cells}"
    if (params.check_vireo_min_sites != null) args << "--min_sites ${params.check_vireo_min_sites}"
    def check_args = args.join(' ')
    def min_dp = params.check_vireo_min_dp != null ? params.check_vireo_min_dp : 5
    """
        cat <<'EOF' > pairs.awk
        function norm(gt,   a, n) {            # 1|0 and 0/1 count as the same genotype
            gsub(/\\|/, "/", gt)
            n = split(gt, a, "/")
            if (n == 2 && a[1] > a[2]) gt = a[2] "/" a[1]
            return gt
        }
        BEGIN {
            FS = "\\t"; OFS = "\\t"
            ns = split(names, nm, ",")
        }
        {
            for (i = 1; i <= ns; i++) {
                gt[i] = \$(2*i - 1)
                dp[i] = \$(2*i)
                ok[i] = (gt[i] !~ /\\./ && dp[i] != "." && dp[i] + 0 >= mindp)
                if (ok[i]) gt[i] = norm(gt[i])
            }
            for (i = 1; i < ns; i++) {
                if (!ok[i]) continue
                for (j = i + 1; j <= ns; j++) {
                    if (!ok[j]) continue
                    key = i SUBSEP j
                    n[key]++
                    if (gt[i] != gt[j]) d[key]++
                }
            }
        }
        END {
            print "sample1", "sample2", "sites_dp", "sites_diff_gt"
            for (i = 1; i < ns; i++)
                for (j = i + 1; j <= ns; j++) {
                    key = i SUBSEP j
                    print nm[i], nm[j], n[key]+0, d[key]+0
                }
        }
EOF

        names=\$(bcftools query -l ${vireo_dir}/GT_donors.vireo.vcf.gz | paste -sd,)

        bcftools query -f '[%GT\\t%DP\\t]\\n' ${vireo_dir}/GT_donors.vireo.vcf.gz | awk -v names="\$names" -v mindp=${min_dp} -f pairs.awk > pairwise_gt_diff.tsv

        check_vireo.py ${vireo_dir}/donor_ids.tsv pairwise_gt_diff.tsv ${check_args}

    """
}
