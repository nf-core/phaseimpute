include { BCFTOOLS_QUERY } from '../../../modules/nf-core/bcftools/query'
include { GAWK           } from '../../../modules/nf-core/gawk'


workflow PREPARE_POSFILE_TSV {

    take:
    ch_panel_sites // channel:   [ [id, chr], vcf, csi ]

    main:

    ch_versions      = Channel.empty()

    // Convert position file to tab-separated file
    BCFTOOLS_QUERY(ch_panel_sites, [], [], [])
    ch_versions = ch_versions.mix(BCFTOOLS_QUERY.out.versions)
    ch_posfile  = BCFTOOLS_QUERY.out.output

    // Remove multiallelic positions from tsv
    GAWK(ch_posfile, [])
    ch_versions = ch_versions.mix(GAWK.out.versions)

    emit:
    posfile  = GAWK.out.output // channel:   [ [id, chr], tsv ]
    versions = ch_versions     // channel:   [ versions.yml ]

}
