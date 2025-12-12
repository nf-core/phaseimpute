include { BCFTOOLS_MPILEUP          } from '../../../modules/nf-core/bcftools/mpileup'
include { BCFTOOLS_MERGE            } from '../../../modules/nf-core/bcftools/merge'
include { BCFTOOLS_ANNOTATE         } from '../../../modules/nf-core/bcftools/annotate'

workflow BAM_GL_BCFTOOLS {

    take:
    ch_bam     // channel: [ [id], bam, bai ]
    ch_posfile // channel: [ [panel_id, chr], posfile_comma]
    ch_fasta   // channel: [ [genome], fasta, fai]

    main:

    ch_versions      = channel.empty()
    ch_multiqc_files = channel.empty()

    ch_mpileup = ch_bam
        .combine(ch_posfile)
        .map{metaI, bam, _bai, metaPC, tsv ->
                [metaI + metaPC, bam, tsv]
        }

    BCFTOOLS_MPILEUP(
        ch_mpileup,
        ch_fasta,
        false
    )
    ch_versions = ch_versions.mix(BCFTOOLS_MPILEUP.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(BCFTOOLS_MPILEUP.out.stats.map{ it -> it[1] })

    // Branch depending on number of files
    ch_all_vcf = BCFTOOLS_MPILEUP.out.vcf
        .join(BCFTOOLS_MPILEUP.out.tbi)
        .map{ metaIPC, vcf, tbi -> [metaIPC.subMap("panel_id", "chr", "batch"), [metaIPC, vcf, tbi]] }
        .groupTuple(sort: { it1, it2 -> it1[0]["id"] <=> it2[0]["id"] }) // Sort by id
        .map{ metaPC, filestups -> [
            metaPC + [id: "all", metas: filestups.collect{it -> it[0]}],
            filestups.collect{it -> it[1]},
            filestups.collect{it -> it[2]},
            filestups.collect{it -> it[1]}.size()
        ] } // Compute number of records
        .branch{it ->
            one: it[3] == 1
            more: it[3] > 1
        }

    // Merge VCFs
    BCFTOOLS_MERGE(
        ch_all_vcf.more.map{it -> [it[0], it[1], it[2], []] },
        ch_fasta
    )
    ch_versions = ch_versions.mix(BCFTOOLS_MERGE.out.versions.first())

    // Mix all vcfs
    ch_to_annotate = ch_all_vcf.one
        .map{it -> [it[0]["metas"][0], it[1][0], it[2][0]] }
        .mix(
            BCFTOOLS_MERGE.out.vcf
                .join(BCFTOOLS_MERGE.out.tbi)
        )

    // Annotate the variants
    BCFTOOLS_ANNOTATE(ch_to_annotate
        .combine(channel.of([[], [], [], []]))
    )
    ch_versions = ch_versions.mix(BCFTOOLS_ANNOTATE.out.versions.first())

    // Output
    ch_output = BCFTOOLS_ANNOTATE.out.vcf
        .join(BCFTOOLS_ANNOTATE.out.tbi)
        .map{ metaIPC, vcf, tbi -> [metaIPC + [ variantcaller:'bcftools' ], vcf, tbi] }

    emit:
    vcf_tbi       = ch_output        // channel: [ [id, panel, chr], vcf, tbi ]
    versions      = ch_versions      // channel: [ versions.yml ]
    multiqc_files = ch_multiqc_files
}
