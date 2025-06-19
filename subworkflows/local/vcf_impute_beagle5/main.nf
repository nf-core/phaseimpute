include { BEAGLE5_BEAGLE                     } from '../../../modules/nf-core/beagle5/beagle'
include { BCFTOOLS_INDEX as VCF_BCFTOOLS_INDEX } from '../../../modules/nf-core/bcftools/index'

workflow VCF_IMPUTE_BEAGLE5 {

    take:
    ch_input        // channel (mandatory): [ [id, panel, chr], vcf, tbi ]
    ch_panel        // channel (mandatory): [ [panel, chr], vcf, tbi ]
    ch_map          // channel (optional):  [ [chr], map]

    main:

    ch_versions = Channel.empty()

    // Prepare input for Beagle imputation
    ch_impute_input = ch_input
        .map{ metaIPC, vcf, index -> [metaIPC, vcf] }

    // Combine with reference panel
    ch_panel_combined = ch_input
        .map{ metaIPC, vcf, index -> [metaIPC.subMap("panel", "chr"), metaIPC, vcf] }
        .combine(ch_panel, by:0)
        .map{ _metaPC, metaIPC, vcf, panel, panel_index -> [metaIPC, vcf, panel] }

    ch_map_file = ch_map.ifEmpty([[]]).collect()
    exclsamples = Channel.of([[]]).collect()
    exclmarkers = Channel.of([[]]).collect()

    BEAGLE5_BEAGLE (
        ch_panel_combined,
        ch_map_file,
        exclsamples,
        exclmarkers,
        Channel.of([[]]).collect()
    )
    ch_versions = ch_versions.mix(BEAGLE5_BEAGLE.out.versions.first())

    VCF_BCFTOOLS_INDEX(BEAGLE5_BEAGLE.out.vcf)
    ch_versions = ch_versions.mix(VCF_BCFTOOLS_INDEX.out.versions.first())

    // Join imputed and index files
    ch_imputed_vcf_tbi = BEAGLE5_BEAGLE.out.vcf
        .join(VCF_BCFTOOLS_INDEX.out.csi)
        .map{ metaIPC, vcf, index -> [metaIPC + [tools: "beagle5"], vcf, index] }

    emit:
    vcf_tbi             = ch_imputed_vcf_tbi    // channel: [ [id, panel, chr, tools], vcf, tbi ]
    versions            = ch_versions           // channel: [ versions.yml ]
}
