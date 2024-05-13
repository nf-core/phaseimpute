include { VCF_PHASE_SHAPEIT5                     } from '../../../subworkflows/nf-core/vcf_phase_shapeit5/main'

workflow VCF_PHASE_PANEL {
    take:
    ch_vcf          // channel: [ [id, chr], vcf, index ]
    ch_panel_norm
    ch_panel_sites
    ch_panel_tsv

    main:

    ch_versions = Channel.empty()

    // Phase panel
    if (params.phased == false) {
        VCF_PHASE_SHAPEIT5(ch_vcf
            .map { meta, vcf, csi -> [meta, vcf, csi, [], meta.region] },
        Channel.of([[],[],[]]).collect(),
        Channel.of([[],[],[]]).collect(),
        Channel.of([[],[]]).collect())
        ch_versions = ch_versions.mix(VCF_PHASE_SHAPEIT5.out.versions)
        ch_panel_phased = VCF_PHASE_SHAPEIT5.out.variants_phased
            .combine(VCF_PHASE_SHAPEIT5.out.variants_index, by: 0)
    } else {
        ch_panel_phased = ch_vcf
    }

    ch_panel = ch_panel_norm
    .combine(ch_panel_sites, by: 0)
    .combine(ch_panel_tsv, by: 0)
    .combine(ch_panel_phased, by: 0)
    .map{ metaIC, norm, n_index, sites, s_index, tsv, t_index, phased, p_index
        -> [[panel:metaIC.id, chr:metaIC.chr ], norm, n_index, sites, s_index, tsv, t_index, phased, p_index]
    }

    emit:
    vcf_tbi             = ch_panel_phased
    panel               = ch_panel
    versions            = ch_versions      // channel: [ versions.yml ]
}
