include { SAMTOOLS_COVERAGE            } from '../../modules/nf-core/samtools/coverage/main.nf'
include { SAMTOOLS_INDEX as INDEX     } from '../../modules/nf-core/samtools/index/main.nf'
include { SAMTOOLS_VIEW as VIEW_REGION } from '../../modules/nf-core/samtools/view/main.nf'
include { SAMTOOLS_VIEW as VIEW_DEPTH  } from '../../modules/nf-core/samtools/view/main.nf'

workflow BAM_DOWNSAMPLE {

    take:
    ch_bam    // channel: [ [id, ref], bam, bai ]
    ch_depth  // channel: [ val(depth) ]
    ch_fasta  // channel: [ fasta ]
    main:

    ch_versions = Channel.empty()

    // Add region to channel
    ch_coverage = ch_bam
        .map{ metaIR, bam, index ->
            [metaIR, bam, index, metaIR["region"]]
        }

    // Get coverage of the region
    SAMTOOLS_COVERAGE ( ch_coverage ) // meta, bam, bai, region
    ch_versions = ch_versions.mix(SAMTOOLS_COVERAGE.out.versions.first())

    // Compute mean depth of the region
    ch_mean_depth = SAMTOOLS_COVERAGE.out.coverage.view()
        .splitCsv(header: true, sep:'\t')
        .map{metaIR, row -> [metaIR,"${row.meandepth}" as Float]}

    // Compute downsampling factor
    ch_depth_factor = ch_mean_depth
        .combine(ch_depth)
        .map{metaIR, mean, depth ->
            [metaIR, metaIR + ["depth":depth], depth as Float / mean + 1]
        }

    // Add all necessary channel for downsampling
    ch_input_downsample = ch_coverage
        .combine(ch_fasta)
        .combine(ch_depth_factor, by:0)
        .map{ metaIR, bam, index, region, fasta, metaIRD, depth ->
            [metaIRD, bam, index, fasta, region, depth]}

    // Downsample
    VIEW_DEPTH(ch_input_downsample, [])
    ch_versions = ch_versions.mix(VIEW_DEPTH.out.versions.first())

    // Index result
    INDEX(VIEW_DEPTH.out.bam)
    ch_versions = ch_versions.mix(INDEX.out.versions.first())

    ch_bam_emul = VIEW_DEPTH.out.bam
        .combine(INDEX.out.bai, by:0)

    emit:
    bam_emul          = ch_bam_emul              // channel: [ metaIRD, bam, bai ]
    versions          = ch_versions              // channel: [ versions.yml ]
}
