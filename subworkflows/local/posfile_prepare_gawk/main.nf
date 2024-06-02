include { GAWK                         } from '../../../modules/nf-core/gawk'
include { GUNZIP                      } from '../../../modules/nf-core/gunzip'

workflow POSFILE_PREPARE_GAWK {

    take:
    ch_posfile // channel:   [ [id, chr], txt ]

    main:
    ch_versions = Channel.empty()

    // Only keep the txt from the channel
    ch_posfile = ch_posfile.map{meta,vcf,txt -> tuple(meta,txt)}

    // Decompress
    GUNZIP(ch_posfile)
    ch_posfile = GUNZIP.out.gunzip

    // Convert TSV in "Glimpse format" to "Stitch format": Replace ","" to "\t"
    GAWK(ch_posfile, [])
    ch_versions = ch_versions.mix(GAWK.out.versions)

    emit:
    posfile  = GAWK.out.output          // channel:   [ [meta], txt ]

}
