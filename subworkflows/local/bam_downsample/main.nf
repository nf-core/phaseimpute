include { SAMTOOLS_COVERAGE   } from '../../../modules/nf-core/samtools/coverage/main.nf'
include { SAMTOOLS_INDEX      } from '../../../modules/nf-core/samtools/index/main.nf'
include { SAMTOOLS_VIEW       } from '../../../modules/nf-core/samtools/view/main.nf'

workflow BAM_DOWNSAMPLE {

    take:
    ch_bam    // channel: [ [id, genome, chr, region], bam, bai ]
    ch_depth  // channel: [ [depth], depth ]
    ch_fasta  // channel: [ [genome], fasta, fai ]

    main:
    ch_versions      = Channel.empty()

    // Add region to channel
    ch_coverage = ch_bam
        .map{ metaICR, bam, index ->
            [ metaICR, bam, index, metaICR["region"] ]
        }

    // Get coverage of the region
    SAMTOOLS_COVERAGE ( ch_coverage, ch_fasta ) // [ meta, bam, bai, region], [ meta, fasta, fai ]
    ch_versions = ch_versions.mix(SAMTOOLS_COVERAGE.out.versions.first())

    // Compute mean depth of the region
    ch_mean_depth = SAMTOOLS_COVERAGE.out.coverage
        .splitCsv(header: true, sep:'\t')
        .map{ metaICR, row ->
            [ metaICR,"${row.meandepth}" as Float ]
        }

    // Compute downsampling factor
    ch_depth_factor = ch_mean_depth
        .combine(ch_depth)
        .map{ metaICR, mean, metaD, depth ->
            [ metaICR, metaICR + metaD, depth as Float / mean ]
        }

    // Add all necessary channel for downsampling
    ch_input_downsample = ch_coverage
        .combine(ch_depth_factor, by : 0)
        .map{ metaICR, bam, index, region, metaICRD, depth ->
            [ metaICRD, bam, index, region, depth ]
        }

    // Downsample
    SAMTOOLS_VIEW(
        ch_input_downsample,
        ch_fasta.map{ metaG, fasta, fai -> [metaG, fasta] },
        []
    )
    ch_versions = ch_versions.mix(SAMTOOLS_VIEW.out.versions.first())

    // Index result
    SAMTOOLS_INDEX(SAMTOOLS_VIEW.out.bam)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    // Aggregate bam and index
    ch_bam_emul = SAMTOOLS_VIEW.out.bam
        .combine(SAMTOOLS_INDEX.out.bai, by:0)

    emit:
    bam_emul          = ch_bam_emul                    // channel: [ [id, genome, chr, region, depth], bam, bai ]
    coverage          = SAMTOOLS_COVERAGE.out.coverage // channel: [ [id, genome, chr, region, depth], txt ]
    versions          = ch_versions                    // channel: [ versions.yml ]
}
