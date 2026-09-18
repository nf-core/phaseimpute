include { GLIMPSE_CHUNK           } from '../../../modules/nf-core/glimpse/chunk'
include { GLIMPSE2_CHUNK          } from '../../../modules/nf-core/glimpse2/chunk'

workflow VCF_CHUNK_GLIMPSE {

    take:
    ch_reference  // channel (mandatory): [ [panel_id, chr], vcf, csi ]
    ch_map        // channel (optional) : [ [panel_id, chr], map ]
    chunk_model   // value : model
    chunk_version // value : glimpse version to use

    main:

    // Add chromosome to channel
    ch_vcf_csi_chr = ch_reference
        .map{metaPC, vcf, csi -> [metaPC, vcf, csi, metaPC.chr]}

    if (chunk_version == "V1") {
        // Make chunks with Glimpse1
        GLIMPSE_CHUNK(ch_vcf_csi_chr)

        ch_chunks = GLIMPSE_CHUNK.out.chunk_chr
    } else if (chunk_version == "V2") {
        ch_input_glimpse2 = ch_vcf_csi_chr
            .combine(ch_map, by:0)

        GLIMPSE2_CHUNK(ch_input_glimpse2, chunk_model)
        ch_chunks = GLIMPSE2_CHUNK.out.chunk_chr
    } else {
        error ("Parameter chunk_version should be V1 or V2, found: ${chunk_version}.")
    }

    emit:
    chunks         = ch_chunks         // channel:  [ [panel_id, chr], txt ]
}
