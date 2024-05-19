include { GLIMPSE2_PHASE                        } from '../../../modules/nf-core/glimpse2/phase'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_1    } from '../../../modules/nf-core/bcftools/index'

workflow VCF_IMPUTE_GLIMPSE2 {

    take:
    ch_input        // channel (mandatory): [ [id], bam, bai ]
    ch_panel        // channel (mandatory): [ [panel, chr, region], vcf, tbi ]
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

    // Impute with Glimpse2
    GLIMPSE2_PHASE(ch_input_glimpse2, ch_fasta)
    ch_versions = ch_versions.mix(GLIMPSE2_PHASE.out.versions)

    // Index phased file
    BCFTOOLS_INDEX_1(GLIMPSE2_PHASE.out.phased_variants)
    ch_versions = ch_versions.mix(BCFTOOLS_INDEX_1.out.versions)

    // Join imputed and index files
    ch_imputed_vcf_tbi = GLIMPSE2_PHASE.out.phased_variants.join(BCFTOOLS_INDEX_1.out.tbi)

    emit:
    vcf_tbi             = ch_imputed_vcf_tbi    // [ [id, chr, region], vcf, tbi ]
    versions            = ch_versions           // channel: [ versions.yml ]
}
