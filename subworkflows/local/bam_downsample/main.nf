include { SAMTOOLS_COVERAGE            } from '../../../modules/nf-core/samtools/coverage/main.nf'
include { SAMTOOLS_INDEX as INDEX      } from '../../../modules/nf-core/samtools/index/main.nf'
include { SAMTOOLS_VIEW as VIEW_REGION } from '../../../modules/nf-core/samtools/view/main.nf'
include { SAMTOOLS_VIEW as VIEW_DEPTH  } from '../../../modules/nf-core/samtools/view/main.nf'

workflow BAM_DOWNSAMPLE {

    take:
    ch_bam    // channel: [ [id, ref], bam, bai ]
    ch_depth  // channel: [ val(depth) ]
    ch_fasta  // channel: [ fasta ]

    main:
    ch_versions = Channel.empty()

    // Add fasta and region to bam channel
    ch_input_region = ch_bam
        .combine(ch_fasta)
        .combine(ch_region)
        .map{ metaI, bam, index, fasta, metaR, region ->
            [ metaI + metaR, bam, index, fasta, region ]
        }
        .combine(Channel.of([[]])) // depth parameter

    // Extract region of interest
    VIEW_REGION(ch_input_region, [])
    ch_versions = ch_versions.mix(VIEW_REGION.out.versions.first())

    // Index region of interest
    INDEX1 (VIEW_REGION.out.bam)
    ch_versions = ch_versions.mix(INDEX1.out.versions.first())

    // Add region to channel
    ch_coverage = ch_bam
        .map{ metaIR, bam, index ->
            [ metaIR, bam, index, metaIR["region"] ]
        }

    // Get coverage of the region
    SAMTOOLS_COVERAGE ( ch_coverage ) // meta, bam, bai, region
    ch_versions = ch_versions.mix(SAMTOOLS_COVERAGE.out.versions.first())

    // Compute mean depth of the region
    ch_mean_depth = SAMTOOLS_COVERAGE.out.coverage
        .splitCsv(header: true, sep:'\t')
        .map{ metaIR, row ->
            [ metaIR,"${row.meandepth}" as Float ]
        }

    // Compute downsampling factor
    ch_depth_factor = ch_mean_depth
        .combine(ch_depth)
        .map{ metaIR, mean, depth ->
            [ metaIR, metaIR + ["depth":depth], depth as Float / mean ]
        }

    // Add all necessary channel for downsampling
    ch_input_downsample = ch_coverage
        .combine(ch_fasta)
        .combine(ch_depth_factor)
        .map{ metaIR, bam, index, region, fasta, metaIRD, depth ->
            [ metaIRD, bam, index, fasta, region, depth ]
        }

    // Downsample
    VIEW_DEPTH(ch_input_downsample, [])
    ch_versions = ch_versions.mix(VIEW_DEPTH.out.versions.first())

    // Index result
    INDEX2(VIEW_DEPTH.out.bam)
    ch_versions = ch_versions.mix(INDEX2.out.versions.first())

    // Aggregate bam and index
    ch_bam_region = VIEW_REGION.out.bam
        .combine(INDEX1.out.bai)
    ch_bam_emul   = VIEW_DEPTH.out.bam
        .combine(INDEX2.out.bai)

    emit:
    bam_region        = ch_bam_region     // channel: [ metaIR, bam, bai ]
    bam_emul          = ch_bam_emul       // channel: [ metaIRD, bam, bai ]
    versions          = ch_versions       // channel: [ versions.yml ]
}
