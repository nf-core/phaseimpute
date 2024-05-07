include { GLIMPSE_CHUNK                     } from '../../../modules/nf-core/glimpse/chunk/main'
include { GLIMPSE2_CHUNK                    } from '../../../modules/nf-core/glimpse2/chunk/main'
include { GLIMPSE2_SPLITREFERENCE           } from '../../../modules/nf-core/glimpse2/splitreference/main'

workflow MAKE_CHUNKS {

    take:
    ch_reference                           // channel: [ val(meta),vcf ]
    ch_map                                 // channel  (optional): [ meta, map ]

    main:

    ch_versions = Channel.empty()

    // Add chromosome to channel
    ch_vcf_csi_chr = ch_reference.map{meta, vcf, csi -> [meta, vcf, csi, meta.chr]}
    ch_vcf_csi_chr.dump(tag:"ch_vcf_csi_chr")

    // Make chunks with Glimpse1
    GLIMPSE_CHUNK(ch_vcf_csi_chr)
    ch_versions = ch_versions.mix(GLIMPSE_CHUNK.out.versions)

    // Rearrange chunks into channel for QUILT
    ch_chunks_glimpse1 = GLIMPSE_CHUNK.out.chunk_chr
                    .splitText()
                    .map { metamap, line ->
                        def fields = line.split("\t")
                        def startEnd = fields[2].split(':')[1].split('-')
                        [metamap, metamap.chr, startEnd[0], startEnd[1]]
                    }
    ch_chunks_glimpse1.dump(tag:"ch_chunks_glimpse1")

    // Make chunks with Glimpse2 (does not work with "sequential" mode)
    chunk_model = "recursive"

    GLIMPSE2_CHUNK ( ch_vcf_csi_chr, ch_map, chunk_model )
    ch_versions = ch_versions.mix( GLIMPSE2_CHUNK.out.versions.first() )

    // Rearrange channels
    ch_chunks_glimpse2 = GLIMPSE2_CHUNK.out.chunk_chr
                            .splitCsv(header: ['ID', 'Chr', 'RegionBuf', 'RegionCnk', 'WindowCm',
                                    'WindowMb', 'NbTotVariants', 'NbComVariants'],
                                    sep: "\t", skip: 0)
                            .map { meta, it -> [meta, it["RegionBuf"], it["RegionCnk"]]}
    ch_chunks_glimpse2.dump(tag:"ch_chunks_glimpse2")

    // Split reference panel in bin files
    // Segmentation fault occurs in small-sized panels
    // Should be run only in full-sized panels

    ch_bins = [[]]

    if (params.binaryref == true) {
    // Create channel to split reference
    split_input = ch_reference.combine(ch_chunks_glimpse2, by: 0)

    // Create a binary reference panel for quick reading time
    GLIMPSE2_SPLITREFERENCE( split_input, ch_map )
    ch_versions = ch_versions.mix( GLIMPSE2_SPLITREFERENCE.out.versions.first() )

    ch_bins = GLIMPSE2_SPLITREFERENCE.out.bin_ref
    }

    emit:
    chunks                    = ch_chunks_glimpse1                    // channel: [ chr, val(meta), start, end, number ]
    chunks_glimpse2           = ch_chunks_glimpse2
    binary                    = ch_bins                               // channel:  [ [meta], bin]
    versions                  = ch_versions                           // channel:  [ versions.yml ]
}
