include { BCFTOOLS_QUERY     } from '../../../modules/nf-core/bcftools/query/main'
include { STITCH             } from '../../../modules/nf-core/stitch/main'

workflow BAM_IMPUTE_STITCH {

    take:
    ch_input                                 //  channel: [ val(meta), bam, bai ]
    ch_panel_sites

    main:

    ch_versions      = Channel.empty()

    // Convert position file to tab-separated file
    BCFTOOLS_QUERY(ch_panel_sites)
    ch_posfile = BCFTOOLS_QUERY.out.output


    // Run STITCH
    STITCH( stitch_input, GET_READS.out, reference, seed )



    emit:
    ch_vcf_tbi                                                          // channel:  [ meta, vcf, tbi ]
    versions                   = ch_versions                           // channel:   [ versions.yml ]

}
