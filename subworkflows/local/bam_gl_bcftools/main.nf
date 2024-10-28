include { GAWK                      } from '../../../modules/nf-core/gawk'
include { TABIX_BGZIP               } from '../../../modules/nf-core/tabix/bgzip'
include { BCFTOOLS_MPILEUP          } from '../../../modules/nf-core/bcftools/mpileup'
include { BCFTOOLS_MERGE            } from '../../../modules/nf-core/bcftools/merge'
include { BCFTOOLS_ANNOTATE         } from '../../../modules/nf-core/bcftools/annotate'

workflow BAM_GL_BCFTOOLS {

    take:
    ch_bam     // channel: [ [id], bam, bai ]
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
    ch_versions = ch_versions.mix(TABIX_BGZIP.out.versions)

    ch_mpileup = ch_bam
        .combine(TABIX_BGZIP.out.output)
        .map{metaI, bam, bai, metaPC, tsv ->
                [metaI + ["panel": metaPC.id, "chr": metaPC.chr], bam, tsv]
        }

    BCFTOOLS_MPILEUP(
        ch_mpileup,
        ch_fasta,
        false
    )
    ch_versions = ch_versions.mix(BCFTOOLS_MPILEUP.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(BCFTOOLS_MPILEUP.out.stats.map{ it[1] })

    // Branch depending on number of files
    ch_all_vcf = BCFTOOLS_MPILEUP.out.vcf
        .join(BCFTOOLS_MPILEUP.out.tbi)
        .map{ metaIPC, vcf, tbi -> [metaIPC.subMap("panel", "chr", "batch"), [metaIPC, vcf, tbi]] }
        .groupTuple(sort: { it1, it2 -> it1[0]["id"] <=> it2[0]["id"] }) // Sort by id
        .map{ metaPC, filestups -> [
            metaPC + [id: "all", metas: filestups.collect{it[0]}],
            filestups.collect{it[1]},
            filestups.collect{it[2]},
            filestups.collect{it[1]}.size()
        ] } // Compute number of records
        .branch{
            one: it[3] == 1
            more: it[3] > 1
        }

    // Merge VCFs
    BCFTOOLS_MERGE(
        ch_all_vcf.more.map{ [it[0], it[1], it[2], []] },
        ch_fasta
    )
    ch_versions = ch_versions.mix(BCFTOOLS_MERGE.out.versions)

    // Mix all vcfs
    ch_to_annotate = ch_all_vcf.one
        .map{ [it[0]["metas"][0], it[1][0], it[2][0]] }
        .mix(
            BCFTOOLS_MERGE.out.vcf
                .join(BCFTOOLS_MERGE.out.tbi)
        )

    // Annotate the variants
    BCFTOOLS_ANNOTATE(ch_to_annotate
        .combine(Channel.of([[], [], [], []]))
    )
    ch_versions = ch_versions.mix(BCFTOOLS_ANNOTATE.out.versions)

    // Output
    ch_output = BCFTOOLS_ANNOTATE.out.vcf
        .join(BCFTOOLS_ANNOTATE.out.tbi)
        .map{ metaIPC, vcf, tbi -> [metaIPC + [ variantcaller:'bcftools' ], vcf, tbi] }

    emit:
    vcf_tbi       = ch_output        // channel: [ [id, panel, chr], vcf, tbi ]
    versions      = ch_versions      // channel: [ versions.yml ]
    multiqc_files = ch_multiqc_files
}
