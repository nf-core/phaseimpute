include { SAMTOOLS_COVERAGE            } from '../../modules/nf-core/samtools/coverage/main.nf'
include { SAMTOOLS_INDEX as INDEX1     } from '../../modules/nf-core/samtools/index/main.nf'
include { SAMTOOLS_INDEX as INDEX2     } from '../../modules/nf-core/samtools/index/main.nf'
include { SAMTOOLS_VIEW as VIEW_REGION } from '../../modules/nf-core/samtools/view/main.nf'
include { SAMTOOLS_VIEW as VIEW_DEPTH  } from '../../modules/nf-core/samtools/view/main.nf'

workflow BAM_SIMULATE {

    take:
    ch_bam    // channel: [ [id, ref], bam, bai ]
    ch_region // channel: [ [ref, region], val(chr), val(start), val(end)]
    ch_depth  // channel: val(depth)
    fasta     // fasta file

    main:

    ch_versions = Channel.empty()

    // Set region to chr:start-end
    ch_region = ch_region.map{ chr, start, end -> [chr + ":" + start + "-" + end]}

    // Add fasta and region to bam channel
    ch_input_region = ch_bam
        .combine(Channel.fromPath(fasta).collect())
        .combine(ch_region)
        .map{ meta, bam, index, fasta, region ->
            [meta + ["region": region], bam, index, fasta, region]
        }
        .combine(Channel.of([[]])) // depth parameter

    // Extract region of interest
    VIEW_REGION(ch_input_region, [])
    ch_versions = ch_versions.mix(VIEW_REGION.out.versions.first())

    // Index region of interest
    INDEX1 (VIEW_REGION.out.bam)
    ch_versions = ch_versions.mix(INDEX1.out.versions.first())

    // Add region to channel
    ch_coverage = VIEW_REGION.out.bam
        .combine(INDEX1.out.bai, by:0)
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
        .map{ metaIR, bam, index, region ->
            [metaIR.subMap(["region"]), metaIR, bam, index, region]}
        .combine(Channel.fromPath(fasta).collect())
        .map{metaR, metaIR, bam, index, region, fasta ->
            [metaIR, bam, index, fasta, region ] }
        .combine(ch_depth_factor, by:0)
        .map{ metaIR, bam, index, fasta, region, metaIRD, depth ->
            [metaIRD, bam, index, fasta, region, depth]}

    // Downsample
    VIEW_DEPTH(ch_input_downsample, [])
    ch_versions = ch_versions.mix(VIEW_DEPTH.out.versions.first())

    // Index result
    INDEX2(VIEW_DEPTH.out.bam)
    ch_versions = ch_versions.mix(INDEX2.out.versions.first())

    emit:
    bam_region        = VIEW_REGION.out.bam
    bam_region_index  = INDEX1.out.bai
    bam_emul          = VIEW_DEPTH.out.bam       // channel: [ val(meta), [ bam ] ]
    bam_emul_index    = INDEX2.out.bai           // channel: [ val(meta), [ bai ] ]
    versions          = ch_versions              // channel: [ versions.yml ]
}
