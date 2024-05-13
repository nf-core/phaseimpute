include { BCFTOOLS_QUERY                    } from '../../../modules/nf-core/bcftools/query/main'
include { GAWK                              } from '../../../modules/nf-core/gawk'


workflow PREPARE_POSFILE_TSV {

    take:
    ch_panel_sites
    ch_fasta

    main:

    ch_versions      = Channel.empty()

    // Convert position file to tab-separated file
    BCFTOOLS_QUERY(ch_panel_sites, [], [], [])
    ch_posfile = BCFTOOLS_QUERY.out.output

    // Remove multiallelic positions from tsv
    GAWK(ch_posfile, [])

    emit:
    posfile                    = GAWK.out.output                       // channel:   [ [id, chr], txt ]
    versions                   = ch_versions                           // channel:   [ versions.yml ]

}
