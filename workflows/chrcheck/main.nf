include { VCF_CHR_EXTRACT as VCF_CHR_BF } from '../../modules/local/vcf_chr_extract'
include { BAM_CHR_EXTRACT as BAM_CHR_BF } from '../../modules/local/bam_chr_extract'
include { BAM_CHR_RENAME_SAMTOOLS       } from '../../subworkflows/local/bam_chr_rename_samtools'
include { VCF_CHR_RENAME_BCFTOOLS       } from '../../subworkflows/local/vcf_chr_rename_bcftools'
include { VCF_CHR_EXTRACT as VCF_CHRAF  } from '../../modules/local/vcf_chr_extract'
include { BAM_CHR_EXTRACT as BAM_CHRAF  } from '../../modules/local/bam_chr_extract'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CHRCHECK {
    take:
        ch_input // [[id], file, index]
        ch_chr   // [chr]

    main:
        ch_versions = Channel.empty()

        ch_input = ch_input.combine(
            ch_chr.unique().collect().toList()
        )

        // Split the input between VCF and BAM files
        ch_input = ch_input.branch{
            bam: it[1] =~ 'bam|cram|sam'
            vcf: it[1] =~ 'vcf|bcf'
        }

        // Extract the contig names from the VCF files
        VCF_CHR_BF(ch_input.vcf.map{ meta, file, index, chr -> [meta, file] })
        ch_versions = ch_versions.mix(VCF_CHR_BF.out.versions)
        chr_vcf_disjoint = checkChr(VCF_CHR_BF.out.chr, ch_input.vcf)

        // Extract the contig names from the BAM files
        BAM_CHR_BF(ch_input.bam.map{ meta, file, index, chr -> [meta, file] })
        ch_versions = ch_versions.mix(BAM_CHR_BF.out.versions)
        chr_bam_disjoint = checkChr(BAM_CHR_BF.out.chr, ch_input.bam)

        if (params.rename_chr == true) {
            // Rename the contigs in the BAM files
            BAM_CHR_RENAME_SAMTOOLS(
                chr_bam_disjoint.to_rename.map{meta, bam, csi, diff, prefix -> [meta, bam, csi, prefix]}
            )
            ch_versions = ch_versions.mix(BAM_CHR_RENAME_SAMTOOLS.out.versions)

            // Rename the contigs in the VCF files
            VCF_CHR_RENAME_BCFTOOLS(
                chr_vcf_disjoint.to_rename
            )
            ch_versions = ch_versions.mix(VCF_CHR_RENAME_BCFTOOLS.out.versions)

            ch_bam_renamed = BAM_CHR_RENAME_SAMTOOLS.out.bam_renamed
            ch_vcf_renamed = VCF_CHR_RENAME_BCFTOOLS.out.vcf_renamed
        } else {
            chr_vcf_disjoint.to_rename.map {
                error "Contig names: ${it[3]} in VCF: ${it[1]} are not present in reference genome with same writing. Please set `rename_chr` to `true` to rename the contigs."
            }
            chr_bam_disjoint.to_rename.map {
                error "Contig names: ${it[3]} in BAM: ${it[1]} are not present in reference genome with same writing. Please set `rename_chr` to `true` to rename the contigs."
            }
            ch_vcf_renamed = Channel.empty()
            ch_bam_renamed = Channel.empty()
        }

        ch_output = chr_bam_disjoint.no_rename
            .mix(chr_vcf_disjoint.no_rename)
            .mix(ch_bam_renamed)
            .mix(ch_vcf_renamed)
    emit:
        output   = ch_output             // [ [id], file, index ]
        versions = ch_versions           // channel: [ versions.yml ]
}

def checkChr(ch_chr, ch_input){
    chr_checked = ch_chr
        .combine(ch_input, by:0)
        .map{metaI, chr, file, index, lst ->
            [
                metaI, file, index,
                chr.readLines()*.split(' ').collect{it[0]},
                lst
            ]
        }
        .branch{ meta, file, index, chr, lst ->
            lst_diff = diffChr(chr, lst, file)
            diff = lst_diff[0]
            prefix = lst_diff[1]
            no_rename: diff.size() == 0
                return [meta, file, index]
            to_rename: true
                return [meta, file, index, diff, prefix]
        }
    return chr_checked
}

def diffChr(chr_target, chr_ref, file) {
    diff = chr_ref - chr_target
    prefix = (chr_ref - chr_target) =~ "chr" ? "chr" : "nochr"
    if (diff.size() != 0) {
        // Ensure that by adding/removing the prefix we can solve the problem
        new_chr = []
        to_rename = []
        if (prefix == "chr") {
            chr_target.each{ new_chr += "chr${it}" }
            diff.each{ to_rename += it.replace('chr', '') }
        } else {
            chr_target.each{ new_chr += it.replace('chr', '') }
            diff.each{ to_rename += "chr${it}" }
        }
        new_diff = diff - new_chr
        if (new_diff.size() != 0) {
            error "Contig names: ${new_diff} absent from file: ${file} and cannot be solved by adding or removing the `chr` prefix."
        }
        diff = to_rename
    }
    return [diff, prefix]
}
