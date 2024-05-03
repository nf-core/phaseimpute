include { GLIMPSE_CHUNK                     } from '../../../modules/nf-core/glimpse/chunk/main'

workflow MAKE_CHUNKS {

    take:
    ch_reference                           // channel: [ val(meta),vcf ]

    main:

    ch_versions = Channel.empty()

    // Make chunks
    ch_vcf_csi_chr = ch_reference.map{meta, vcf, csi -> [meta, vcf, csi, meta.chr]}
    GLIMPSE_CHUNK(ch_vcf_csi_chr)

    // Rearrange chunks into channel
    ch_chunks = GLIMPSE_CHUNK.out.chunk_chr
                    .splitText()
                    .map { metamap, line ->
                        def fields = line.split("\t")
                        def startEnd = fields[2].split(':')[1].split('-')
                        [metamap, metamap.chr, startEnd[0], startEnd[1]]
                    }

    emit:
    ch_chunks                 = ch_chunks                              // channel: [ chr, val(meta), start, end, number ]
    versions                  = ch_versions                           // channel:  [ versions.yml ]
}
