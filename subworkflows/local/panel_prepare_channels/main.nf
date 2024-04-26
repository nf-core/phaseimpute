include { VCF_CONCATENATE_BCFTOOLS as CONCAT_PANEL   } from '../../../subworkflows/local/vcf_concatenate_bcftools'

workflow PANEL_PREPARE_CHANNELS {
    take:
    ch_panel_norm          // channel: [ [id, chr], vcf, index ]
    ch_panel_sites
    ch_panel_tsv
    ch_panel_phased

    main:

    ch_versions = Channel.empty()

    ch_panel = ch_panel_norm
        .combine(ch_panel_sites, by: 0)
        .combine(ch_panel_tsv, by: 0)
        .combine(ch_panel_phased, by: 0)
        .map{ metaIC, norm, n_index, sites, s_index, tsv, t_index, phased, p_index
            -> [[panel:metaIC.id, chr:metaIC.chr ], norm, n_index, sites, s_index, tsv, t_index, phased, p_index]
        }


        ch_panel_sites_tsv = ch_panel
            .map{ metaPC, norm, n_index, sites, s_index, tsv, t_index, phased, p_index
                -> [metaPC, sites, tsv]
            }
        CONCAT_PANEL(ch_panel
            .map{ metaPC, norm, n_index, sites, s_index, tsv, t_index, phased, p_index
                -> [[id:metaPC.panel], sites, s_index]
            }
        )
        ch_panel_sites = CONCAT_PANEL.out.vcf_tbi_join

        ch_panel_phased = ch_panel_phased
            .map{ metaPC, norm, n_index, sites, s_index, tsv, t_index, phased, p_index
                -> [metaPC, phased, p_index]
            }


    emit:
    panel      = ch_panel
    ch_panel_sites
    ch_panel_phased
    ch_panel_sites_tsv

    versions       = ch_versions      // channel: [ versions.yml ]
}
