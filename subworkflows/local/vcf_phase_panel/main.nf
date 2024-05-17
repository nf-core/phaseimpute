include { VCF_PHASE_SHAPEIT5                     } from '../../../subworkflows/nf-core/vcf_phase_shapeit5'

workflow VCF_PHASE_PANEL {
    take:
    ch_vcf          // channel: [ [id, chr, region], vcf, index ]
    ch_panel_norm   // channel: [ [panel, chr], norm, index ]
    ch_panel_sites  // channel: [ [panel, chr], sites, index ]
    ch_panel_tsv    // channel: [ [panel, chr], tsv, index ]

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

    ch_panel = ch_panel_norm
        .combine(ch_panel_sites, by: 0)
        .combine(ch_panel_tsv, by: 0)
        .combine(ch_panel_phased, by: 0)
        .map{ metaPC, norm, n_index, sites, s_index, tsv, t_index, phased, p_index
            -> [[panel:metaPC.id, chr:metaPC.chr ], norm, n_index, sites, s_index, tsv, t_index, phased, p_index]
        }

    emit:
    vcf_tbi             = ch_panel_phased  // channel: [ [id, chr], vcf, index ]
    panel               = ch_panel         // channel: [ [panel, chr], norm, n_index, sites, s_index, tsv, t_index, phased, p_index ]
    versions            = ch_versions      // channel: [ versions.yml ]
}
