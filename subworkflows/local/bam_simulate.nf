include { SAMTOOLS_COVERAGE            } from '../../modules/nf-core/samtools/coverage/main.nf'
include { SAMTOOLS_INDEX as INDEX1     } from '../../modules/nf-core/samtools/index/main.nf'
include { SAMTOOLS_INDEX as INDEX2     } from '../../modules/nf-core/samtools/index/main.nf'
include { SAMTOOLS_VIEW as VIEW_REGION } from '../../modules/nf-core/samtools/view/main.nf'
include { SAMTOOLS_VIEW as VIEW_DEPTH  } from '../../modules/nf-core/samtools/view/main.nf'

workflow BAM_SIMULATE {

    take:
    ch_bam    // channel: [ [id, ref], bam, bai ]
    ch_region // channel: [ [ref, region], fasta, val(region)]
    ch_depth  // channel: val(depth)

    main:

    ch_versions = Channel.empty()

    ch_input_region = ch_bam
                            .combine(ch_region)
                            .map{ metaIR, bam, index, metaRR, fasta, region ->
                                [metaIR + metaRR, bam, index, fasta, region]}
                            .combine(Channel.of([[]]))

    VIEW_REGION(ch_input_region, [])
    ch_versions = ch_versions.mix(VIEW_REGION.out.versions.first())

    INDEX1 (VIEW_REGION.out.bam)
    ch_versions = ch_versions.mix(INDEX1.out.versions.first())

    ch_coverage = VIEW_REGION.out.bam
                    .combine(INDEX1.out.bai, by:0)
                    .map{metaIRR, bam, index -> [metaIRR, bam, index, metaIRR.region]}

    SAMTOOLS_COVERAGE ( ch_coverage )
    ch_versions = ch_versions.mix(SAMTOOLS_COVERAGE.out.versions.first())

    ch_mean_depth = SAMTOOLS_COVERAGE.out.coverage
                    .splitCsv(header: true, sep:'\t')
                    .map{metaIRR, row -> [metaIRR,"${row.meandepth}" as Float]}

    ch_depth_factor = ch_mean_depth
                        .combine(ch_depth)
                        .map{metaIRR, mean, depth -> [metaIRR, metaIRR + ["depth":depth], depth as Float / mean + 1]}

    ch_input_downsample = ch_coverage
                            .map{ metaIRR, bam, index, region ->
                                    [metaIRR.subMap(["ref","region"]), metaIRR, bam, index]}
                            .combine(ch_region.map{metaRR, fasta, region ->  [metaRR, fasta, region]},
                                    by: 0)
                            .map{metaRR, metaIRR, bam, index, fasta, region ->
                                    [metaIRR, bam, index, fasta, region ] }
                            .combine(ch_depth_factor, by:0)
                            .map{ metaIRR, bam, index, fasta, region, metaIRRD, depth ->
                                    [metaIRRD, bam, index, fasta, region, depth]}

    VIEW_DEPTH(ch_input_downsample, [])
    ch_versions = ch_versions.mix(VIEW_DEPTH.out.versions.first())

    INDEX2(VIEW_DEPTH.out.bam)

    emit:
    bam_region        = VIEW_REGION.out.bam
    bam_region_index  = INDEX1.out.bai
    bam_emul          = VIEW_DEPTH.out.bam       // channel: [ val(meta), [ bam ] ]
    bam_emul_index    = INDEX2.out.bai           // channel: [ val(meta), [ bai ] ]
    versions          = ch_versions              // channel: [ versions.yml ]
}
