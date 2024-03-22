include { BCFTOOLS_ANNOTATE           } from '../../../modules/nf-core/bcftools/annotate/main.nf'
include { BCFTOOLS_INDEX              } from '../../../modules/nf-core/bcftools/index/main.nf'
include { FAITOCHR                    } from '../../../modules/local/faitochr/main.nf'
include { VCFCHREXTRACT               } from '../../../modules/local/vcfchrextract/main.nf'

workflow VCF_CHR_RENAME {
    take:
    ch_vcf          // channel: [ [id], vcf, index ]
    ch_fasta        // channel: [ [id], fasta, fai ]

    main:

    ch_versions = Channel.empty()

    // Get contig names from the VCF
    VCFCHREXTRACT(ch_vcf.map{ metaV, vcf, csi -> [metaV, vcf] })

    // Check if the contig names are the same as the reference
    chr_disjoint = VCFCHREXTRACT.out.chr
        .combine(ch_vcf, by:0)
        .combine(ch_fasta)
        .map{metaI, chr, vcf, csi, metaG, fasta, fai ->
            [
                metaI, vcf, csi,
                chr.readLines()*.split(' ').collect{it[0]},
                fai.readLines()*.split('\t').collect{it[0]}
            ]
        }
        .map { meta, vcf, csi, chr, fai ->
            [meta, vcf, csi, (chr-fai).size()]
        }
        .branch{
            no_rename: it[3] == 0
            to_rename: it[3] > 0
        }

    if (chr_disjoint.to_rename.ifEmpty(true) != true){
        if (params.rename_chr == true) {
            println 'Some contig names in the VCF do not match the reference genome. Renaming the contigs by adding / removing "chr" prefix ...'
            // Generate the chromosome renaming file
            FAITOCHR(ch_fasta.map{ metaG, fasta, fai -> [metaG, fai] })
            ch_versions = ch_versions.mix(FAITOCHR.out.versions)

            // Rename the chromosome without prefix
            BCFTOOLS_ANNOTATE(chr_disjoint.to_rename.map{ meta, vcf, csi, chr -> [meta, vcf, csi] }
                .combine(Channel.of([[], [], []]))
                .combine(FAITOCHR.out.annot_chr.map{it[1]})
            )
            ch_versions = ch_versions.mix(BCFTOOLS_ANNOTATE.out.versions.first())

            BCFTOOLS_INDEX(BCFTOOLS_ANNOTATE.out.vcf)
            ch_versions = ch_versions.mix(BCFTOOLS_INDEX.out.versions.first())

            ch_vcf_renamed = BCFTOOLS_ANNOTATE.out.vcf
                .combine(BCFTOOLS_INDEX.out.csi, by:0)

            ch_vcf_out = chr_disjoint.no_rename
                .map{meta, vcf, csi, chr -> [meta, vcf, csi]}
                .mix(ch_vcf_renamed)
        } else {
            error 'Some contig names in the VCF do not match the reference genome. Please set `rename_chr` to `true` to rename the contigs.'
        }
    } else {
        ch_vcf_out = ch_vcf
    }

    emit:
    vcf            = ch_vcf_out            // [ meta, vcf, csi ]
    versions       = ch_versions           // channel: [ versions.yml ]
}
