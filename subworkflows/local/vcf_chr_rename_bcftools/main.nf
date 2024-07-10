include { BCFTOOLS_ANNOTATE           } from '../../../modules/nf-core/bcftools/annotate'
include { BCFTOOLS_INDEX              } from '../../../modules/nf-core/bcftools/index'

workflow VCF_CHR_RENAME_BCFTOOLS {
    take:
    ch_vcf          // channel: [ [id], vcf, index, diff, prefix ]

    main:

    ch_versions = Channel.empty()

    // Check that prefix is either "chr" or "nochr"
    ch_vcf = ch_vcf.map{
        meta, vcf, index, diff, prefix ->
        if (prefix != "chr" && prefix != "nochr") {
            error "Invalid chr_prefix: ${prefix}"
        }
        [meta, vcf, index, diff, prefix]
    }

    // Generate the chromosome renaming file
    ch_rename_file = ch_vcf
        .collectFile{ meta, vcf, index, diff, prefix ->
            def chr = diff.collect { i ->
                prefix == "chr" ? "${i} chr${i}" :
                "${i} ${i.replace('chr', '')}"
            }.join('\n')
            ["${meta.id}.txt", chr]
        }
        .map{ file -> [[id: file.getBaseName()], file] }

    // Add the chromosome renaming file to the input channel
    ch_annotate_input = ch_vcf.map {
        meta, vcf, index, diff, prefix ->
        [[id: meta.id], meta, vcf, index]
    } // channel: [ [id], vcf, index ]
    .combine(ch_rename_file, by:0)
    .map {
        metaI, meta, vcf, index, rename_file ->
        [meta, vcf, index, [], [], [], rename_file]
    }

    // Rename the chromosome without prefix
    BCFTOOLS_ANNOTATE(ch_annotate_input)
    ch_versions = ch_versions.mix(BCFTOOLS_ANNOTATE.out.versions.first())

    BCFTOOLS_INDEX(BCFTOOLS_ANNOTATE.out.vcf)
    ch_versions = ch_versions.mix(BCFTOOLS_INDEX.out.versions.first())

    ch_vcf_renamed = BCFTOOLS_ANNOTATE.out.vcf
        .combine(BCFTOOLS_INDEX.out.csi, by:0)

    emit:
    vcf_renamed    = ch_vcf_renamed        // [ [id], vcf, csi ]
    versions       = ch_versions           // channel: [ versions.yml ]
}
