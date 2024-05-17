include { GLIMPSE2_PHASE                 } from '../../../modules/nf-core/glimpse2/phase'
include { GLIMPSE2_LIGATE                } from '../../../modules/nf-core/glimpse2/ligate'
include { BCFTOOLS_INDEX as INDEX_PHASE  } from '../../../modules/nf-core/bcftools/index'
include { BCFTOOLS_INDEX as INDEX_LIGATE } from '../../../modules/nf-core/bcftools/index'

workflow VCF_IMPUTE_GLIMPSE2 {

    take:
    ch_input        // channel (mandatory): [ [id], vcf, csi, infos ]
    ch_panel        // channel (mandatory): [ [panel, chr, region], vcf, csi, region ]
    ch_chunks       // channel  (optional): [ [chr], region1, region2 ]
    ch_fasta        // channel (mandatory): [ [genome], fa, fai ]

    main:

    ch_versions = Channel.empty()

    // Impute with Glimpse2 without using binary files
    def samples_file = [[]]
    def gmap = [[]]

    // Create input channel to impute with Glimpse2
    // Add chr as key to input
    ch_input = ch_input.map{meta, bam, bai -> return[['chr': meta.chr], meta, bam, bai]}

    // Join chunks and panel
    ch_chunks_panel = ch_chunks.join(ch_panel)

    // Change key:value names
    ch_chunks_panel = ch_chunks_panel.map{meta, vcf, csi, region1, region2 -> return[['id': meta.panel, 'chr': meta.chr], vcf, csi, region1, region2]}

    // Add chr as key
    ch_chunks_panel = ch_chunks_panel.map{meta, vcf, csi, region1, region2 -> return[['chr': meta.chr], vcf, csi, region1, region2]}

    // Join input and chunks reference
    ch_input_glimpse2 = ch_input.map { it + samples_file }.join(ch_chunks_panel).map { it + gmap }

    // Remove chr key
    ch_input_glimpse2 = ch_input_glimpse2.map{ it[1..-1] }

    //Impute with Glimpse2
    GLIMPSE2_PHASE(ch_input_glimpse2, ch_fasta) // Error: AC/AN INFO fields in VCF are inconsistent with GT field, update the values in the VCF
    ch_versions = ch_versions.mix(GLIMPSE2_PHASE.out.versions)

    emit:
    versions               = ch_versions                            // channel: [ versions.yml ]
}
