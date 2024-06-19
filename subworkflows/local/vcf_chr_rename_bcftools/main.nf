include { BCFTOOLS_ANNOTATE           } from '../../../modules/nf-core/bcftools/annotate'
include { BCFTOOLS_INDEX              } from '../../../modules/nf-core/bcftools/index'
include { GAWK                        } from '../../../modules/nf-core/gawk'

workflow VCF_CHR_RENAME_BCFTOOLS {
    take:
    ch_vcf          // channel: [ [id], vcf, index, dif, prefix ]

    main:

    ch_versions = Channel.empty()

    // Generate the chromosome renaming file
    ch_rename_file = ch_vcf
        .map{ meta, vcf, index, diff, prefix ->
            if (prefix == "chr") :
                return [
                    meta,
                    diff.collectFile{["rename_chr.txt", "${it} chr${it}\n"]}
                ]
            else :
                return [
                    meta,
                    diff.collectFile{["rename_nochr.txt", "${it} ${it.replaceFirst("chr", "")}\n"]}
                ]
        }

    // Rename the chromosome without prefix
    BCFTOOLS_ANNOTATE(
        ch_vcf // channel: [ [id], vcf, index ]
            .combine(Channel.of([[],[],[]]))
            .combine(ch_rename_file, by:0)
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
