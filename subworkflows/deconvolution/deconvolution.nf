include {DONOR_ASSIGNMENT} from '../donor_assignment/donor_assignment.nf'
include { MPILEUP; MPILEUP as MPILEUP_ND; CELLSNP; VIREO; SUBSET_VCF; SUBSET_VCF as SUBSET_VCF_ND; CHECK_VIREO }  from '../../modules/deconvolution.nf'
//include {GET_BARCODES; SPLIT_BAM; INDEX_SPLIT_BAM} from '../../modules/split_bam.nf'
include {SPLIT_BAM_SINTO; INDEX_SPLIT_BAM; INDEX_SPLIT_BAM as INDEX_SPLIT_BAM_SUPPLEMENTARY; INDEX_SPLIT_BAM as INDEX_SPLIT_BAM_NOSUPPLEMENTARY} from '../../modules/split_bam.nf'
include {MOKING_VIREO} from '../../modules/moking_vireo.nf'

workflow DECONVOLUTION {
  take:
    mapped_reads
  main:
    Channel.empty().set { ch_versions }
    if (mapped_reads == 'independent workflow') {
    Channel
      .fromPath(params.input_samples_path)
      .splitCsv(sep: ',', header: true)
      //.filter { it -> !(it.sample_id in excluded_samples_list) }
      .map { it ->
        def sample_id = it.sample_id
        def bam_path = "${params.results_output}qc/mapped/${sample_id}.mapped.realcells_only.bam"
        def bai_path = "${params.results_output}qc/mapped/${sample_id}.mapped.realcells_only.bam.bai"
        [sample_id, bam_path, bai_path]
      }
      .set { mapped_reads }
    }

    Channel
      .fromPath(params.input_samples_path)
      .splitCsv(sep: ',', header: true)
      .map { it -> [it.sample_id, it.nr_samples_multiplexed] } // Create a tuple with bam_path and sample_id
      .set { sampleNames_nrDons }

    // Only perform deconvolution on samples that are more than 1 donor.
    combined_ch_deconv =  mapped_reads.combine(sampleNames_nrDons, by: 0)
    combined_ch_deconv.filter { experiment, bam, bai, npooled -> npooled == '1' }.map { experiment, bam, bai, _ -> [experiment, bam, bai] }.set{not_for_deconv}
    combined_ch_deconv.filter { experiment, bam, bai, npooled -> npooled != '1' }.map { experiment, bam, bai, _ -> [experiment, bam, bai] }.set{for_deconv}
    for_deconv.view { experiment, bam, bai -> "Deconvolution will be performed on ${experiment} with bam: ${bam} and bai: ${bai}" }
    not_for_deconv.view { experiment, bam, bai -> "Deconvolution will NOT be performed on ${experiment} with bam: ${bam} and bai: ${bai}" }

    MPILEUP(for_deconv, params.genome_fasta_f)
    SUBSET_VCF(MPILEUP.out, params.subset_regions_bed)
    CELLSNP(SUBSET_VCF.out)

    sampleNames_cellsnp_nrDons = CELLSNP.out.combine(sampleNames_nrDons, by: 0)
    VIREO(sampleNames_cellsnp_nrDons)
    CHECK_VIREO(VIREO.out.vireo_results)
    CHECK_VIREO.out.check_res
    .branch { sample_id, check_file ->
        passed: check_file.name.startsWith('passed')
        failed: check_file.name.startsWith('failed')
    }
    .set { vireo_check }
    passed_ids = vireo_check.passed.map { sid, f -> sid }
    failed_ids = vireo_check.failed.map { sid, f -> sid }
    failed_ids.view { sample_id -> "Sample ${sample_id} failed vireo check and will be excluded from downstream analysis." }
    mapped_reads = mapped_reads.combine(VIREO.out.barcode_list, by: 0)
    if (params.vireo_check=="TRUE"){
      mapped_reads = mapped_reads.join(passed_ids)
    }
    mapped_reads.view { sample_id, bam, bai, barcodes -> "Sample ${sample_id} passed vireo check files ${bam}, ${bai}, ${barcodes}." }
    SPLIT_BAM_SINTO(mapped_reads)

    donor_bams = SPLIT_BAM_SINTO.out.donor_bams.flatMap { val, paths ->
      paths.collect { p -> tuple(p.simpleName, p) } }
    INDEX_SPLIT_BAM(donor_bams)
    fullBam_ch_pre = INDEX_SPLIT_BAM.out.donor_bams
    fullBam_ch = fullBam_ch_pre.mix(not_for_deconv)

    donor_supplementary_bams = SPLIT_BAM_SINTO.out.donor_supplementary_bams.flatMap { val, paths ->
      paths.collect { p -> tuple(p.simpleName, p) } }
    INDEX_SPLIT_BAM_SUPPLEMENTARY(donor_supplementary_bams)
    supplementaryBam_ch = INDEX_SPLIT_BAM_SUPPLEMENTARY.out.donor_bams

    donor_nosupplementary_bams = SPLIT_BAM_SINTO.out.donor_nosupplementary_bams.flatMap { val, paths ->
      paths.collect { p -> tuple(p.simpleName, p) } }
    INDEX_SPLIT_BAM_NOSUPPLEMENTARY(donor_nosupplementary_bams)
    nosupplementaryBam_ch = INDEX_SPLIT_BAM_NOSUPPLEMENTARY.out.donor_bams
    if (params.run_gtcheck== 'TRUE' || params.run_barcode_check == 'TRUE') {
        MPILEUP_ND(not_for_deconv, params.genome_fasta_f)
        SUBSET_VCF_ND(MPILEUP_ND.out, params.subset_regions_bed)
        MOKING_VIREO(SUBSET_VCF_ND.out)
        if (params.vireo_check=="TRUE"){
        sample_vcf= VIREO.out.sample_donor_vcf.join(passed_ids).mix(MOKING_VIREO.out.sample_donor_vcf)
        sample_donor_ids=VIREO.out.sample_donor_ids.join(passed_ids).mix(MOKING_VIREO.out.sample_donor_ids)
        }else{
        sample_vcf= VIREO.out.sample_donor_vcf.mix(MOKING_VIREO.out.sample_donor_vcf)
        sample_donor_ids=VIREO.out.sample_donor_ids.mix(MOKING_VIREO.out.sample_donor_ids)
        }
        sample_vcf.view { sample_id, vcf -> "VCF for ${sample_id} is ${vcf}" }
        sample_donor_ids.view { sample_id, file_path -> "Donor IDs for ${sample_id} are in ${file_path}" }
        sample_donor_ids.collectFile(name: 'vireo_map.tsv', newLine: true) { sample_id, file_path ->
          "${sample_id}\t${file_path}"
        }
        .set { vireo_map }

        DONOR_ASSIGNMENT(vireo_map, sample_vcf)

        //SUBSET_VCF_ND.out.map{ experiment, bam, bai, barcodes, vcf -> [experiment, vcf] }.set{nd_vcfs}
        //nd_vcfs.view { experiment, vcf -> "VCF for ${experiment} is ${vcf}" }
/*
        if (for_deconv.count()>0) {
          sample_vcf=nd_vcfs.mix(VIREO.out.sample_donor_vcf)
          nd_vcfs.mix(VIREO.out.sample_donor_ids).collectFile(name: 'vireo_map_nd.tsv', newLine: true) { sample_id, file_path ->
            "${sample_id}\t${file_path}"
          }
          .set { vireo_map }
        }
        else{
          sample_vcf=nd_vcfs
          nd_vcfs.collectFile(name: 'vireo_map_nd.tsv', newLine: true) { sample_id, file_path ->
            "${sample_id}\t${file_path}"
          }
          .set { vireo_map }
        }
        //sample_vcf=VIREO.out.sample_donor_vcf
        //sample_donor_ids=VIREO.out.sample_donor_ids
        sample_donor_ids.collectFile(name: 'vireo_map.tsv', newLine: true) { sample_id, file_path ->
          "${sample_id}\t${file_path}"
        }
        .set { vireo_map }

        DONOR_ASSIGNMENT(vireo_map, sample_vcf)
*/    }

  emit:
    fullBam_ch
    supplementaryBam_ch
    nosupplementaryBam_ch
}
