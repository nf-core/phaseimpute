include { GLIMPSE_CHUNK           } from '../../../modules/nf-core/glimpse/chunk'
include { GLIMPSE2_CHUNK          } from '../../../modules/nf-core/glimpse2/chunk'

workflow VCF_CHUNK_GLIMPSE {

    take:
    ch_reference  // channel (mandatory): [ [panel_id, chr], vcf, csi ]
    ch_map        // channel (optional) : [ [panel_id, chr], map ]
    chunk_model   // channel : model

    main:

    ch_versions = Channel.empty()
    // Add chromosome to channel
    ch_vcf_csi_chr = ch_reference
        .map{metaPC, vcf, csi -> [metaPC, vcf, csi, metaPC.chr]}

    // Make chunks with Glimpse1
    GLIMPSE_CHUNK(ch_vcf_csi_chr)
    ch_versions = ch_versions.mix(GLIMPSE_CHUNK.out.versions.first())

    // Rearrange chunks into channel for QUILT
    ch_chunks_quilt = GLIMPSE_CHUNK.out.chunk_chr
        .splitText()
        .map { metaPC, line ->
            def fields = line.split("\t")
            def startEnd = fields[2].split(':')[1].split('-')
            [metaPC, metaPC.chr, startEnd[0], startEnd[1]]
        }

    // Rearrange chunks into channel for GLIMPSE1 and GLIMPSE2
    ch_chunks_glimpse1 = GLIMPSE_CHUNK.out.chunk_chr
        .splitCsv(
            header: ['ID', 'Chr', 'RegionIn', 'RegionOut', 'Size1', 'Size2'],
            sep: "\t", skip: 0
        )
        .map { metaPC, it -> [metaPC, it["RegionIn"], it["RegionOut"]]}

    ch_input_glimpse2 = ch_vcf_csi_chr
        .combine(ch_map, by:0)

    GLIMPSE2_CHUNK(ch_input_glimpse2, chunk_model)
    ch_versions = ch_versions.mix(GLIMPSE2_CHUNK.out.versions.first())

    // Rearrange channels
    ch_chunks_glimpse2 = GLIMPSE2_CHUNK.out.chunk_chr
        .splitCsv(
            header: [
                'ID', 'Chr', 'RegionBuf', 'RegionCnk', 'WindowCm',
                'WindowMb', 'NbTotVariants', 'NbComVariants'
            ], sep: "\t", skip: 0
        )
        .map { metaPC, it -> [metaPC, it["RegionBuf"], it["RegionCnk"]]}

    emit:
    chunks                    = GLIMPSE_CHUNK.out.chunk_chr           // channel:  [ [panel_id, chr], txt ]
    chunks_quilt              = ch_chunks_quilt                       // channel:  [ [panel_id, chr], chr,  start, end ]
    chunks_glimpse1           = ch_chunks_glimpse1                    // channel:  [ [panel_id, chr], chr,  region1, region2 ]
    chunks_glimpse2           = ch_chunks_glimpse2                    // channel:  [ [panel_id, chr], chr,  region1, region2 ]
    versions                  = ch_versions                           // channel:  [ versions.yml ]
}
