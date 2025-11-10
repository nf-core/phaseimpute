include { MINIMAC4_COMPRESSREF } from '../../../modules/nf-core/minimac4/compressref/main'
include { MINIMAC4_IMPUTE      } from '../../../modules/nf-core/minimac4/impute/main'
include { BCFTOOLS_INDEX       } from '../../../modules/nf-core/bcftools/index/main'

workflow VCF_IMPUTE_MINIMAC4 {

    take:
    ch_input     // channel: [ [id, chr], vcf, tbi ]
    ch_panel     // channel: [ [id, chr], vcf, tbi ]
    ch_map       // channel: [ [chr], map]
    ch_posfile   // channel: [ [chr], sites_vcf, sites_index, hap, legend ]

    main:

    ch_versions = Channel.empty()

    ch_posfile_minimac4 = ch_posfile
        .map { meta, sites_vcf, sites_index, hap, legend ->
            [meta, sites_vcf, sites_index]
        }

    // Compress reference panel to MSAV format
    MINIMAC4_COMPRESSREF(ch_panel)
    ch_versions = ch_versions.mix(MINIMAC4_COMPRESSREF.out.versions.first())

    // Prepare input channels for MINIMAC4
    ch_minimac4_input = ch_input
        .map { meta, vcf, tbi -> [meta.chr, meta, vcf, tbi] }
        .combine(
            MINIMAC4_COMPRESSREF.out.msav.map { meta, msav -> [meta.chr, meta.id, msav] },
            by: 0
        )
        .combine(
            ch_map.map { meta, map -> [meta.chr, map] },
            by: 0
        )
        .combine(
            ch_posfile_minimac4.map { meta, sites_vcf, sites_index ->
                [meta.chr, sites_vcf, sites_index]
            },
            by: 0
        )
        .map { chr, target_meta, target_vcf, target_tbi, panel_id, ref_msav, map, sites_vcf, sites_index ->
            [target_meta + [panel: panel_id], target_vcf, target_tbi, ref_msav, sites_vcf, sites_index, map]
        }
    // Perform imputation
    MINIMAC4_IMPUTE(ch_minimac4_input)
    ch_versions = ch_versions.mix(MINIMAC4_IMPUTE.out.versions.first())

    // Index the output VCF file
    BCFTOOLS_INDEX(
        MINIMAC4_IMPUTE.out.vcf
    )
    ch_versions = ch_versions.mix(BCFTOOLS_INDEX.out.versions.first())

    // Join imputed and index files
    ch_vcf_index = MINIMAC4_IMPUTE.out.vcf
        .join(
            BCFTOOLS_INDEX.out.tbi
                .mix(BCFTOOLS_INDEX.out.csi)
        )
        .map{ meta, vcf, index -> [meta + [tools: "minimac4"], vcf, index] }

    emit:
    vcf_index  = ch_vcf_index // channel: [ [id, chr, tools], vcf, index ]
    versions = ch_versions        // channel: [ versions.yml ]
}
