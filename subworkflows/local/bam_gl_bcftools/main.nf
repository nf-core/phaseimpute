include { GAWK                      } from '../../../modules/nf-core/gawk'
include { TABIX_BGZIP               } from '../../../modules/nf-core/tabix/bgzip'
include { BCFTOOLS_MPILEUP          } from '../../../modules/nf-core/bcftools/mpileup'
include { BCFTOOLS_ANNOTATE         } from '../../../modules/nf-core/bcftools/annotate'

workflow BAM_GL_BCFTOOLS {

    take:
    ch_input   // channel: [ [id], bam, bai ]
    ch_posfile // channel: [ [panel, chr], legend]
    ch_fasta   // channel: [ [genome], fasta, fai]

    main:

    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // Convert legend to TSV with ','
    GAWK(ch_posfile, [])
    ch_versions = ch_versions.mix(GAWK.out.versions)

    // Compress TSV
    TABIX_BGZIP(GAWK.out.output)
    ch_versions = ch_versions.mix(TABIX_BGZIP.out.versions.first())

    ch_mpileup       = ch_input
        .combine(TABIX_BGZIP.out.output)
        .map{metaI, bam, bai, metaPC, tsv ->
                [metaI + ["panel": metaPC.id, "chr": metaPC.chr], bam, tsv]
        }

    BCFTOOLS_MPILEUP(
        ch_mpileup,
        ch_fasta,
        false
    )
    ch_versions = ch_versions.mix(BCFTOOLS_MPILEUP.out.versions.first())

    // Annotate the variants
    BCFTOOLS_ANNOTATE(BCFTOOLS_MPILEUP.out.vcf
        .join(BCFTOOLS_MPILEUP.out.tbi)
        .combine(Channel.of([[], [], [], []]))
    )
    ch_versions = ch_versions.mix(BCFTOOLS_ANNOTATE.out.versions.first())

    // Output
    ch_output = BCFTOOLS_ANNOTATE.out.vcf
        .join(BCFTOOLS_ANNOTATE.out.tbi)

    ch_multiqc_files = ch_multiqc_files.mix(BCFTOOLS_MPILEUP.out.stats.map{ it[1] })

    emit:
    vcf           = ch_output        // channel: [ [id, panel, chr], vcf, tbi ]
    versions      = ch_versions      // channel: [ versions.yml ]
    multiqc_files = ch_multiqc_files
}
