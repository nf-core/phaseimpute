//
// Subworkflow to perform genotype imputation using Minimac4
//

include { MINIMAC4_COMPRESSREF } from '../../../modules/nf-core/minimac4/compressref/main'
include { MINIMAC4_IMPUTE }      from '../../../modules/nf-core/minimac4/impute/main'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_MINIMAC4 } from '../../../modules/nf-core/bcftools/index/main'

workflow VCF_IMPUTE_MINIMAC4 {

    take:
    ch_input  // channel: [ [id, chr], vcf, tbi ]
    ch_panel  // channel: [ [id, chr], vcf, tbi ]  
    ch_map    // channel: [ [chr], map]

    main:
    
    ch_versions = Channel.empty()

    // Compress reference panel to MSAV format (spécifique à MINIMAC4)
    MINIMAC4_COMPRESSREF(
        ch_panel
    )
    ch_versions = ch_versions.mix(MINIMAC4_COMPRESSREF.out.versions.first())

    // Prepare input channels for MINIMAC4 by combining VCF, panel, and map files
    ch_minimac4_input = ch_input
        .map { meta, vcf, tbi -> [meta.chr, meta, vcf, tbi] }
        .combine(MINIMAC4_COMPRESSREF.out.msav.map { meta, msav -> [meta.chr, meta.id, msav] }, by: 0)
        .combine(ch_map.map { meta, map -> [meta.chr, map] }, by: 0)
        .map { chr, target_meta, target_vcf, target_tbi, panel_id, ref_msav, map ->
            [target_meta, target_vcf, target_tbi, ref_msav, [], [], map]
        }

    // Perform imputation 
    MINIMAC4_IMPUTE(
        ch_minimac4_input
    )
    ch_versions = ch_versions.mix(MINIMAC4_IMPUTE.out.versions.first())

    // Index the output VCF file
    BCFTOOLS_INDEX_MINIMAC4(
        MINIMAC4_IMPUTE.out.vcf
    )
    ch_versions = ch_versions.mix(BCFTOOLS_INDEX_MINIMAC4.out.versions.first())

    // Join imputed and index files (même structure que BEAGLE5)
    ch_imputed_vcf_tbi = MINIMAC4_IMPUTE.out.vcf
        .join(BCFTOOLS_INDEX_MINIMAC4.out.csi)
        .map{ meta, vcf, index -> [meta + [tools: "minimac4"], vcf, index] }

    emit:
    vcf_tbi  = ch_imputed_vcf_tbi // channel: [ [id, chr, tools], vcf, tbi ]
    versions = ch_versions        // channel: [ versions.yml ]
}