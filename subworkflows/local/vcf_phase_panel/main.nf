include { VCF_PHASE_SHAPEIT5                     } from '../../../subworkflows/nf-core/vcf_phase_shapeit5'

workflow VCF_PHASE_PANEL {
    take:
    ch_vcf          // channel: [ [id, chr, region], vcf, index ]

    main:

    ch_versions = Channel.empty()

    // Phase panel
    if (params.phased == false) {
        VCF_PHASE_SHAPEIT5(ch_vcf
            .map { metaICR, vcf, csi -> [metaICR, vcf, csi, [], metaICR.region] },
        Channel.of([[],[],[]]).collect(),
        Channel.of([[],[],[]]).collect(),
        Channel.of([[],[]]).collect())
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
