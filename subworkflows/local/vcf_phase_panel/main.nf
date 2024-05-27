include { VCF_PHASE_SHAPEIT5                     } from '../../../subworkflows/nf-core/vcf_phase_shapeit5'

workflow VCF_PHASE_PANEL {
    take:
    ch_vcf          // channel: [ [id, chr], vcf, index ]
    ch_region       // channel: [ [chr, region], region ]

    main:

    ch_versions = Channel.empty()

    // Phase panel
    if (params.phased == false) {
        ch_phase_input = ch_vcf
            .map { metaICR, vcf, csi -> [metaICR.subMap("chr"), metaICR, vcf, csi] }
            .combine(
                ch_region.map{metaCR, region -> [metaCR.subMap("chr"), region]}
                , by: 0)
            .map { metaC, metaICR, vcf, csi, region -> [metaICR.subMap("id"), vcf, csi, [], region] }
            .view()
        VCF_PHASE_SHAPEIT5(
            ch_phase_input,
            Channel.of([[],[],[]]).collect(),
            Channel.of([[],[],[]]).collect(),
            Channel.of([[],[]]).collect()
        )
        ch_versions = ch_versions.mix(VCF_PHASE_SHAPEIT5.out.versions)
        ch_panel_phased = VCF_PHASE_SHAPEIT5.out.variants_phased
            .combine(VCF_PHASE_SHAPEIT5.out.variants_index, by: 0)
    } else {
        ch_panel_phased = ch_vcf
    }

    emit:
    vcf_tbi             = ch_panel_phased  // channel: [ [id, chr], vcf, index ]
    versions            = ch_versions      // channel: [ versions.yml ]
}
