//
// Subworkflow to perform genotype imputation using Minimac4
//

include { MINIMAC4_COMPRESSREF } from '../../../modules/nf-core/minimac4/compressref/main'
include { MINIMAC4_IMPUTE }      from '../../../modules/nf-core/minimac4/impute/main'

workflow VCF_IMPUTE_MINIMAC4 {

    take:
    ch_target       // channel: [ val(meta), path(target_vcf), path(target_index) ]
    ch_reference    // channel: [ val(meta), path(ref_vcf), path(ref_index) ]
    ch_sites        // channel: [ val(meta), path(sites_vcf), path(sites_index) ] or [ [], [], [] ]
    ch_map          // channel: [ val(meta), path(map) ] or [ [], [] ]

    main:
    
    ch_versions = Channel.empty()

    //
    // MODULE: Compress reference panel to MSAV format
    //
    MINIMAC4_COMPRESSREF (
        ch_reference
    )
    ch_versions = ch_versions.mix(MINIMAC4_COMPRESSREF.out.versions)

    //
    // MODULE: Perform imputation
    //
    // Combine target data with compressed reference
    ch_impute_input = ch_target
        .combine(MINIMAC4_COMPRESSREF.out.msav.map{ meta, msav -> msav })
        .combine(ch_sites.map{ meta, sites_vcf, sites_index -> [sites_vcf, sites_index] })
        .combine(ch_map.map{ meta, map -> map })
        .map{ target_meta, target_vcf, target_index, msav, sites_vcf, sites_index, map -> 
            [target_meta, target_vcf, target_index, msav, sites_vcf, sites_index, map]
        }

    MINIMAC4_IMPUTE (
        ch_impute_input
    )
    ch_versions = ch_versions.mix(MINIMAC4_IMPUTE.out.versions)

    emit:
    vcf      = MINIMAC4_IMPUTE.out.vcf      // channel: [ val(meta), path(vcf) ]
    versions = ch_versions                   // channel: [ path(versions.yml) ]
}