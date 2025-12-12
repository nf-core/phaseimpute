include { GLIMPSE2_PHASE                     } from '../../../modules/nf-core/glimpse2/phase'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_1 } from '../../../modules/nf-core/bcftools/index'
include { GLIMPSE2_LIGATE                    } from '../../../modules/nf-core/glimpse2/ligate'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_2 } from '../../../modules/nf-core/bcftools/index'

workflow BAM_VCF_IMPUTE_GLIMPSE2 {

    take:
    ch_input        // channel (mandatory): [ [id], [file], [index], bamlist ]
    ch_panel        // channel (mandatory): [ [panel_id, chr], vcf, index ]
    ch_chunks       // channel  (optional): [ [panel_id, chr], region1, region2 ]
    ch_fasta        // channel (mandatory): [ [genome], fa, fai ]
    ch_map          // channel (mandatory): [ [chr], map ]

    main:

    ch_versions = channel.empty()

    // Impute with Glimpse2 without using binary files
    samples_file = channel.of([[]]).collect()

    // Create input channel to impute with Glimpse2

    // Join chunks, panel and map
    ch_chunks_panel_map = ch_chunks
        .combine(ch_panel, by:0)
        .map{ metaPC, regionin, regionout, panel, panel_index ->
            [
                ["chr": metaPC.chr], ["panel": metaPC.id, "chr": metaPC.chr, "chunk": regionout],
                regionin, regionout, panel, panel_index
            ]
        }
        .combine(ch_map, by:0)
        .map{ _metaC, metaPCC, regionin, regionout, panel, panel_index, gmap ->
            [metaPCC, regionin, regionout, panel, panel_index, gmap]
        }
        .ifEmpty{
            error "BAM_VCF_IMPUTE_GLIMPSE2: join operation resulted in an empty channel. Please provide a valid ch_chunks and ch_map channel as input."
        }

    // Join input and chunks reference
    ch_phase_input = ch_input
        .combine(samples_file)
        .combine(ch_chunks_panel_map)
        .map{ metaI, bam, bai, bamlist, samples, metaPCC, regionin, regionout, panel, panel_index, gmap ->
            [metaI + metaPCC,
            bam, bai, bamlist, samples, regionin, regionout, panel, panel_index, gmap]
        }

    // Impute with Glimpse2
    GLIMPSE2_PHASE(ch_phase_input, ch_fasta)
    ch_versions = ch_versions.mix( GLIMPSE2_PHASE.out.versions.first() )

    // Index phased file
    BCFTOOLS_INDEX_1(GLIMPSE2_PHASE.out.phased_variants)
    ch_versions = ch_versions.mix( BCFTOOLS_INDEX_1.out.versions.first() )

    // Ligate all phased files in one and index it
    ligate_input = GLIMPSE2_PHASE.out.phased_variants
        .join( BCFTOOLS_INDEX_1.out.csi )
        .map{ metaIPCR, vcf, index -> [metaIPCR.subMap("id", "panel_id", "chr", "batch"), vcf, index] }
        .groupTuple()

    GLIMPSE2_LIGATE( ligate_input )
    ch_versions = ch_versions.mix( GLIMPSE2_LIGATE.out.versions.first() )

    BCFTOOLS_INDEX_2( GLIMPSE2_LIGATE.out.merged_variants )
    ch_versions = ch_versions.mix( BCFTOOLS_INDEX_2.out.versions.first() )

    // Join imputed and index files
    ch_imputed_vcf_tbi = GLIMPSE2_LIGATE.out.merged_variants
        .join(BCFTOOLS_INDEX_2.out.tbi)
        .map{ metaIPC, vcf, index -> [metaIPC + [tools: "glimpse2"], vcf, index] }

    emit:
    vcf_tbi             = ch_imputed_vcf_tbi    // channel: [ [id, panel_id, chr, tool], vcf, tbi ]
    versions            = ch_versions           // channel: [ versions.yml ]
}
