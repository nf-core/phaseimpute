include { VCFCHREXTRACT as VCFCHRBFR  } from '../../../modules/local/vcfchrextract/main.nf'
include { VCFCHREXTRACT as VCFCHRAFT  } from '../../../modules/local/vcfchrextract/main.nf'
include { VCF_CHR_RENAME              } from '../vcf_chr_rename/main.nf'

workflow VCF_CHR_CHECK {
    take:
    ch_vcf          // channel: [ [id], vcf, index ]
    ch_fasta        // channel: [ [id], fasta, fai ]

    main:

    ch_versions = Channel.empty()

    // Get contig names from the VCF
    VCFCHRBFR(ch_vcf.map{ metaV, vcf, csi -> [metaV, vcf] })
    ch_versions = ch_versions.mix(VCFCHRBFR.out.versions.first())

    // Check if the contig names are the same as the reference
    chr_disjoint = check_chr(VCFCHRBFR.out.chr, ch_vcf, ch_fasta)

    if (params.rename_chr == true) {
        // Generate the chromosome renaming file
        VCF_CHR_RENAME(
            chr_disjoint.to_rename.map{meta, vcf, index, nb -> [meta, vcf, index]},
            ch_fasta
        )
        ch_versions = ch_versions.mix(VCF_CHR_RENAME.out.versions.first())

        // Check if modification has solved the problem
        VCFCHRAFT(VCF_CHR_RENAME.out.vcf_renamed.map{ metaV, vcf, csi -> [metaV, vcf] })
        ch_versions = ch_versions.mix(VCFCHRAFT.out.versions.first())

        chr_disjoint_after = check_chr(VCFCHRAFT.out.chr, VCF_CHR_RENAME.out.vcf_renamed, ch_fasta)

        chr_disjoint_after.to_rename.map{
            error 'Even after renaming errors are still present. Please check that contigs name in vcf and fasta file are equivalent.'
        }
        ch_vcf_renamed = VCF_CHR_RENAME.out.vcf_renamed

    } else {
        chr_disjoint.to_rename.map {
            error 'Some contig names in the VCF do not match the reference genome. Please set `rename_chr` to `true` to rename the contigs.'
        }
        ch_vcf_renamed = Channel.empty()
    }    

    ch_vcf_out = chr_disjoint.no_rename
        .map{meta, vcf, csi, chr -> [meta, vcf, csi]}
        .mix(ch_vcf_renamed)

    emit:
    vcf            = ch_vcf_out            // [ meta, vcf, csi ]
    versions       = ch_versions           // channel: [ versions.yml ]
}


def check_chr(ch_chr, ch_vcf, ch_fasta){
    chr_checked = ch_chr
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
    return chr_checked
}