include { SAMTOOLS_DEPTH                     } from '../../../modules/nf-core/samtools/depth'
include { SAMTOOLS_VIEW                      } from '../../../modules/nf-core/samtools/view'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_1 } from '../../../modules/nf-core/samtools/index'
include { SAMTOOLS_MERGE                     } from '../../../modules/nf-core/samtools/merge'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_2 } from '../../../modules/nf-core/samtools/index'

workflow BAM_DOWNSAMPLE {

    take:
    ch_bam    // channel: [ [id, genome, chr, region], bam, bai ]
    ch_depth  // channel: [ [depth], depth ]
    ch_fasta  // channel: [ [genome], fasta, fai ]

    main:
    ch_versions      = Channel.empty()

    // Compute mean depth
    SAMTOOLS_DEPTH(ch_bam, [[], []])
    ch_mean_depth = SAMTOOLS_DEPTH.out.tsv
        .splitCsv(header: false, sep:'\t')
        .map{ metaICR, row ->
            [ metaICR, row[2] as Float ]
        }
        .groupTuple()
        .map{ metaICR, depth ->
            [ metaICR, depth.sum()/depth.size() ]
        }
    ch_versions = ch_versions.mix(SAMTOOLS_DEPTH.out.versions.first())

    // Compute downsampling factor
    ch_depth_factor = ch_mean_depth
        .combine(ch_depth)
        .map{ metaICR, mean, metaD, depth ->
            [ metaICR, metaICR + metaD, depth as Float / mean ]
        }

    // Add all necessary channel for downsampling
    ch_input_downsample = ch_bam
        .combine(ch_depth_factor, by : 0)
        .map{ metaICR, bam, index, metaICRD, depth ->
            [ metaICRD, bam, index, [], depth ]
        }

    // Downsample
    SAMTOOLS_VIEW(
        ch_input_downsample,
        ch_fasta.map{ metaG, fasta, fai -> [metaG, fasta] },
        []
    )
    ch_versions = ch_versions.mix(SAMTOOLS_VIEW.out.versions.first())

    // Index result
    SAMTOOLS_INDEX_1(SAMTOOLS_VIEW.out.bam)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX_1.out.versions.first())

    // Aggregate bam and index
    ch_bam_emul = SAMTOOLS_VIEW.out.bam
        .combine(SAMTOOLS_INDEX_1.out.bai, by:0)

    if (params.input_region) {
        SAMTOOLS_MERGE(
            ch_bam_emul
                .map{
                    metaICRD, bam, index -> [metaICRD.subMap("id", "depth"), bam, index]
                }
                .groupTuple()
                .map{ metaID, bam, index ->
                    [ metaID + ["chr": "all"], bam, index ]
                },
            ch_fasta
        )
        ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions.first())

        SAMTOOLS_INDEX_2(SAMTOOLS_MERGE.out.bam)
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX_2.out.versions.first())

        ch_bam_emul_all = SAMTOOLS_MERGE.out.bam
            .combine(SAMTOOLS_INDEX_2.out.bai, by:0)
    } else {
        ch_bam_emul_all = ch_bam_emul
    }

    emit:
    bam_emul          = ch_bam_emul_all                // channel: [ [id, chr, region, depth], bam, bai ]
    versions          = ch_versions                    // channel: [ versions.yml ]
}
