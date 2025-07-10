//
// Subworkflow to perform genotype imputation using Minimac4
//

include { MINIMAC4_COMPRESSREF } from '../../../modules/nf-core/minimac4/compressref/main'
include { MINIMAC4_IMPUTE }      from '../../../modules/nf-core/minimac4/impute/main'
include { BCFTOOLS_INDEX }       from '../../../modules/nf-core/bcftools/index/main'

workflow VCF_IMPUTE_MINIMAC4 {

    take:
    ch_target       // channel: [ val(meta), path(target_vcf), path(target_index) ]
    ch_reference    // channel: [ val(meta), path(ref_vcf), path(ref_index) ]
    ch_map          // channel: [ val(meta), path(map) ] or [ [], [] ]

    main:
    
    ch_versions = Channel.empty()

    // Compress reference panel to MSAV format
    MINIMAC4_COMPRESSREF (
        ch_reference
    )
    ch_versions = ch_versions.mix(MINIMAC4_COMPRESSREF.out.versions)

    // DEBUG: View compressed reference output
    MINIMAC4_COMPRESSREF.out.msav.view { "VCF_IMPUTE_MINIMAC4 - Compressed reference: $it" }

    // Prepare input channels for MINIMAC4 by combining VCF, panel, and map files
    ch_impute_input = ch_target
        .map { meta, vcf, index -> [meta.chr, meta, vcf, index] }
        .combine(MINIMAC4_COMPRESSREF.out.msav.map { meta, msav -> [meta.chr, msav] }, by: 0)
        .map { chr, target_meta, target_vcf, target_index, msav ->
            [target_meta, target_vcf, target_index, msav, [], [], []]
        }

    // DEBUG: View what each element represents
    ch_impute_input.view { meta, vcf, index, msav, sites, sites_idx, map ->
        "MINIMAC4_IMPUTE INPUT: meta=${meta}, vcf=${vcf}, index=${index}, msav=${msav}, sites=${sites}, sites_idx=${sites_idx}, map=${map}"
    }

    // Perform imputation 
    MINIMAC4_IMPUTE (
        ch_impute_input
    )
    ch_versions = ch_versions.mix(MINIMAC4_IMPUTE.out.versions)

    // Index the output VCF file
    BCFTOOLS_INDEX (
        MINIMAC4_IMPUTE.out.vcf
    )
    ch_versions = ch_versions.mix(BCFTOOLS_INDEX.out.versions)

    // Join imputed and index files
    ch_imputed_vcf_tbi = MINIMAC4_IMPUTE.out.vcf
        .join(BCFTOOLS_INDEX.out.csi)
        .map{ meta, vcf, index -> [meta + [tools: "minimac4"], vcf, index] }

    // DEBUG: View what goes to CONCAT - with file details
    ch_imputed_vcf_tbi.view { meta, vcf, index ->
        "OUTPUT TO CONCAT: meta=${meta}, vcf=${vcf.name} (${vcf.size()} bytes), index=${index.name} (${index.size()} bytes)"
    }

    emit:
    vcf_tbi  = ch_imputed_vcf_tbi             // channel: [ [id, chr, tools], vcf, csi ]
    versions = ch_versions                    // channel: [ path(versions.yml) ]
}