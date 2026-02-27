include { BCFTOOLS_MPILEUP                 } from '../../../modules/nf-core/bcftools/mpileup'
include { BCFTOOLS_MERGE                   } from '../../../modules/nf-core/bcftools/merge'
include { BCFTOOLS_ANNOTATE                } from '../../../modules/nf-core/bcftools/annotate'
include { VCF_GATHER_BCFTOOLS              } from '../../../subworkflows/nf-core/vcf_gather_bcftools'

workflow BAM_GL_BCFTOOLS {

    take:
    ch_bam     // channel: [ [id], bam, bai ]
    ch_posfile // channel: [ [panel_id, chr], posfile_comma]
    ch_fasta   // channel: [ [genome], fasta, fai]

    main:

    ch_versions      = channel.empty()
    ch_multiqc_files = channel.empty()

    def posfile_count = ch_posfile
        .map{ _meta, posfile -> posfile}
        .collect()
        .map { posfile -> posfile.size() }

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
    ch_all_samples_vcf = BCFTOOLS_MPILEUP.out.vcf
        .join(BCFTOOLS_MPILEUP.out.tbi)
        .map{ metaIPC, vcf, tbi -> [metaIPC.subMap("panel_id", "chr", "batch"), [metaIPC, vcf, tbi]] }
        .groupTuple(sort: { it1, it2 -> it1[0]["id"] <=> it2[0]["id"] }) // Sort by id
        .map{ metaPC, filestups -> [
            metaPC + [id: "all_samples", metas: filestups.collect{it -> it[0]}],
            filestups.collect{it -> it[1]},
            filestups.collect{it -> it[2]},
            filestups.collect{it -> it[1]}.size()
        ] } // Compute number of records
        .branch{ _meta, _vcf, _index, size ->
            one: size == 1
            more: size > 1
        }

    // Merge VCFs all individuals
    BCFTOOLS_MERGE(
        ch_all_samples_vcf.more.map{ meta, vcf, index, _size -> [ meta, vcf, index, [] ] },
        ch_fasta
    )
    ch_versions = ch_versions.mix(BCFTOOLS_MERGE.out.versions.first())

    // Mix all vcfs
    ch_to_concat = ch_all_samples_vcf.one
        .map{ meta, vcf_list, index_list, _size -> [ meta, vcf_list[0], index_list[0] ] }
        .mix(
            BCFTOOLS_MERGE.out.vcf
                .join(BCFTOOLS_MERGE.out.tbi.mix(
                    BCFTOOLS_MERGE.out.csi
                ))
        )

    // Merge all chromosomes
    VCF_GATHER_BCFTOOLS(
        ch_to_concat.combine(posfile_count),
        ["id", "panel_id"],
        false
    )

    // Annotate the variants
    BCFTOOLS_ANNOTATE(VCF_GATHER_BCFTOOLS.out.vcf_index
        .combine(channel.of([[], [], [], []]))
    )
    ch_versions = ch_versions.mix(BCFTOOLS_ANNOTATE.out.versions.first())

    // Output
    ch_output = BCFTOOLS_ANNOTATE.out.vcf
        .join(BCFTOOLS_ANNOTATE.out.tbi.mix(
            BCFTOOLS_ANNOTATE.out.csi
        ))
        .map{ metaIPC, vcf, index -> [metaIPC + [ variantcaller:'bcftools' ], vcf, index] }

    emit:
    vcf_index     = ch_output        // channel: [ [id, panel], vcf, index ]
    versions      = ch_versions      // channel: [ versions.yml ]
    multiqc_files = ch_multiqc_files
}
