include { BCFTOOLS_MPILEUP          } from '../../../modules/nf-core/bcftools/mpileup'
include { BCFTOOLS_INDEX            } from '../../../modules/nf-core/bcftools/index'
include { BCFTOOLS_ANNOTATE         } from '../../../modules/nf-core/bcftools/annotate'

workflow BAM_GL_BCFTOOLS {

    take:
    ch_input   // channel: [ [id], bam, bai ]
    ch_target  // channel: [ [panel, chr], sites, tsv]
    ch_fasta   // channel: [ [genome], fasta, fai]

    main:

    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()

    ch_mpileup       = ch_input
        .combine(ch_target)
        .map{metaI, bam, bai, metaPC, sites, tsv ->
                [metaI + ["panel": metaPC.id, "chr": metaPC.chr], bam, sites, tsv]
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

    // Index annotated VCF
    BCFTOOLS_INDEX(BCFTOOLS_ANNOTATE.out.vcf)
    ch_versions = ch_versions.mix(BCFTOOLS_INDEX.out.versions.first())

    // Output
    ch_output = BCFTOOLS_ANNOTATE.out.vcf
        .join(BCFTOOLS_INDEX.out.tbi)

    ch_multiqc_files = ch_multiqc_files.mix(BCFTOOLS_MPILEUP.out.stats.map{ it[1] })

    emit:
    vcf           = ch_output        // channel: [ [id, panel, chr], vcf, tbi ]
    versions      = ch_versions      // channel: [ versions.yml ]
    multiqc_files = ch_multiqc_files
}
