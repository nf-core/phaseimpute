include { BAMCHREXTRACT as BAMCHRBFR  } from '../../../modules/local/bamchrextract'
include { GAWK                        } from '../../../modules/nf-core/gawk'
include { BAM_CHR_RENAME_SAMTOOLS     } from '../bam_chr_rename_samtools'
include { BAMCHREXTRACT as BAMCHRAFT  } from '../../../modules/local/bamchrextract'

workflow BAM_CHR_CHECK_SAMTOOLS {
    take:
    ch_bam          // channel: [ [id], bam, bai ]
    ch_fasta        // channel: [ [genome], fasta, fai ]

    main:

    ch_versions = Channel.empty()

    // Get contig names from the BAM
    BAMCHRBFR(ch_bam.map{ meta, bam, index -> [meta, bam] })
    ch_versions = ch_versions.mix(BAMCHRBFR.out.versions)

    // Check if the contig names are the same as the reference
    chr_disjoint = check_chr(BAMCHRBFR.out.chr, ch_bam, ch_fasta)

    if (params.rename_chr == true) {
        // Generate the chromosome renaming file
        BAM_CHR_RENAME_SAMTOOLS(
            chr_disjoint.to_rename.map{meta, bam, csi, dis, prefix -> [meta, bam, csi, prefix]}
        )
        ch_versions = ch_versions.mix(BAM_CHR_RENAME_SAMTOOLS.out.versions)

        // Check if modification has solved the problem
        BAMCHRAFT(BAM_CHR_RENAME_SAMTOOLS.out.bam_renamed.map{ meta, bam, csi -> [meta, bam] })
        ch_versions = ch_versions.mix(BAMCHRAFT.out.versions)

        chr_disjoint_after = check_chr(BAMCHRAFT.out.chr, BAM_CHR_RENAME_SAMTOOLS.out.bam_renamed, ch_fasta)
        chr_disjoint_after.to_rename.map{
            error "Even after renaming errors are still present. Please check the contigs name : ${it[3]} in bam and fasta file."
        }
        ch_bam_renamed = BAM_CHR_RENAME_SAMTOOLS.out.bam_renamed

    } else {
        chr_disjoint.to_rename.map {
            error "Contig names: ${it[3]} in the BAM are not present in the reference genome. Please set `rename_chr` to `true` to rename the contigs."
        }
        ch_bam_renamed = Channel.empty()
    }

    ch_bam_out = chr_disjoint.no_rename
        .mix(ch_bam_renamed)

    emit:
    bam            = ch_bam_out            // [ [id], bam, csi ]
    versions       = ch_versions           // channel: [ versions.yml ]
}


def check_chr(ch_chr, ch_bam, ch_fasta){
    chr_checked = ch_chr
        .combine(ch_bam, by:0)
        .combine(ch_fasta)
        .map{metaI, chr, bam, csi, metaG, fasta, fai ->
            [
                metaI, bam, csi,
                chr.readLines()*.split(' ').collect{it[0]},
                fai.readLines()*.split('\t').collect{it[0]}
            ]
        }
        .branch{ meta, bam, csi, chr, fai ->
            no_rename: (chr - fai).size() == 0
                return [meta, bam, csi]
            to_rename: true
                return [meta, bam, csi, chr-fai, chr-fai =~ "chr" ? "nochr" : "chr"]
        }
    return chr_checked
}
