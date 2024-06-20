include { BCFTOOLS_ANNOTATE           } from '../../../modules/nf-core/bcftools/annotate'
include { BCFTOOLS_INDEX              } from '../../../modules/nf-core/bcftools/index'
include { GAWK                        } from '../../../modules/nf-core/gawk'

workflow VCF_CHR_RENAME_BCFTOOLS {
    take:
    ch_vcf          // channel: [ [id], vcf, index, diff, prefix ]

    main:

    ch_versions = Channel.empty()

    // Generate the chromosome renaming file
    ch_rename_file = ch_vcf
        .collectFile{ meta, vcf, index, diff, prefix ->
            chr = ""
            if (prefix == "chr") {
                for (i in diff) {
                    chr += "${i} chr${i}\n"
                }
            } else if (prefix == "nochr") {
                for (i in diff) {
                    chr += "${i} ${i.replace('chr', '')}\n"
                }
            } else {
                error "Unknown prefix: ${prefix}"
            }
            ["${meta.id}.txt", chr]
        }
        .map{ file -> [[id: file.getBaseName()], file] }

    // Rename the chromosome without prefix
    BCFTOOLS_ANNOTATE(
        ch_vcf
            .map {
                meta, vcf, index, diff, prefix ->
                [[id: meta.id], meta, vcf, index]
            } // channel: [ [id], vcf, index ]
            .combine(ch_rename_file, by:0)
            .map {
                metaI, meta, vcf, index, rename_file ->
                [meta, vcf, index, [], [], [], rename_file]
            }
    )
    ch_versions = ch_versions.mix(BCFTOOLS_ANNOTATE.out.versions.first())

    BCFTOOLS_INDEX(BCFTOOLS_ANNOTATE.out.vcf)
    ch_versions = ch_versions.mix(BCFTOOLS_INDEX.out.versions.first())

    ch_vcf_renamed = BCFTOOLS_ANNOTATE.out.vcf
        .combine(BCFTOOLS_INDEX.out.csi, by:0)

    emit:
    vcf_renamed    = ch_vcf_renamed        // [ [id], vcf, csi ]
    versions       = ch_versions           // channel: [ versions.yml ]
}
