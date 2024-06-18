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
        // Generate the chromosome renaming command based on the fasta file
        GAWK(ch_fasta.map{ metaG, fasta, fai -> [metaG, fai] }, [])
        ch_versions = ch_versions.mix(GAWK.out.versions)

        chr_prefix = GAWK.out.output
            .splitCsv()
            .map{ it[1][0] }

        // Generate the chromosome renaming file
        BAM_CHR_RENAME_SAMTOOLS(
            chr_disjoint.to_rename.map{meta, bam, index, nb -> [meta, bam, index]},
            chr_prefix
        )
        ch_versions = ch_versions.mix(BAM_CHR_RENAME_SAMTOOLS.out.versions)

        // Check if modification has solved the problem
        BAMCHRAFT(BAM_CHR_RENAME_SAMTOOLS.out.bam_renamed.map{ meta, bam, csi -> [meta, bam] })
        ch_versions = ch_versions.mix(BAMCHRAFT.out.versions)

        chr_disjoint_after = check_chr(BAMCHRAFT.out.chr, BAM_CHR_RENAME_SAMTOOLS.out.bam_renamed, ch_fasta)
        chr_disjoint_after.to_rename.view()
        chr_disjoint_after.to_rename.map{
            error 'Even after renaming errors are still present. Please check that contigs name in bam and fasta file are equivalent.'
        }
        ch_bam_renamed = BAM_CHR_RENAME_SAMTOOLS.out.bam_renamed

    } else {
        chr_disjoint.to_rename.map {
            error 'Some contig names in the BAM do not match the reference genome. Please set `rename_chr` to `true` to rename the contigs.'
        }
        ch_bam_renamed = Channel.empty()
    }

    ch_bam_out = chr_disjoint.no_rename
        .map{meta, bam, csi, chr -> [meta, bam, csi]}
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
        .map { meta, bam, csi, chr, fai ->
            [meta, bam, csi, (chr-fai).size()]
        }
        .branch{
            no_rename: it[3] == 0
            to_rename: it[3] > 0
        }
    return chr_checked
}
