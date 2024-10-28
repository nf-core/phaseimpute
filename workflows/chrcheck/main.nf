include { VCF_CHR_EXTRACT         } from '../../modules/local/vcf_chr_extract'
include { BAM_CHR_EXTRACT         } from '../../modules/local/bam_chr_extract'
include { BAM_CHR_RENAME_SAMTOOLS } from '../../subworkflows/local/bam_chr_rename_samtools'
include { VCF_CHR_RENAME_BCFTOOLS } from '../../subworkflows/local/vcf_chr_rename_bcftools'
include { checkChr                } from '../../subworkflows/local/utils_nfcore_chrcheck_pipeline'
include { diffChr                 } from '../../subworkflows/local/utils_nfcore_chrcheck_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CHRCHECK {
    take:
        ch_input // [[id], file, index, [chr]]

    main:
        ch_versions = Channel.empty()
        // Split the input between VCF and BAM files
        ch_input = ch_input.branch{
            bam: it[1] =~ 'bam|cram|sam'
            vcf: it[1] =~ 'vcf|bcf'
        }
        // Check if channel is empty
        chr_vcf_disjoint = Channel.empty()
        // Extract the contig names from the VCF files
        VCF_CHR_EXTRACT(ch_input.vcf.map{ meta, file, index, chr -> [meta, file] })
        ch_versions = ch_versions.mix(VCF_CHR_EXTRACT.out.versions)
        chr_vcf_disjoint = checkChr(VCF_CHR_EXTRACT.out.chr, ch_input.vcf)

        chr_bam_disjoint = Channel.empty()
        // Extract the contig names from the BAM files
        BAM_CHR_EXTRACT(ch_input.bam.map{ meta, file, index, chr -> [meta, file] })
        ch_versions = ch_versions.mix(BAM_CHR_EXTRACT.out.versions)
        chr_bam_disjoint = checkChr(BAM_CHR_EXTRACT.out.chr, ch_input.bam)

        if (params.rename_chr == true) {
            ch_bam_renamed = Channel.empty()
            // Rename the contigs in the BAM files
            BAM_CHR_RENAME_SAMTOOLS(
                chr_bam_disjoint.to_rename.map{meta, bam, csi, diff, prefix -> [meta, bam, csi, prefix]}
            )
            ch_versions = ch_versions.mix(BAM_CHR_RENAME_SAMTOOLS.out.versions)
            ch_bam_renamed = BAM_CHR_RENAME_SAMTOOLS.out.bam_renamed

            ch_vcf_renamed = Channel.empty()
            // Rename the contigs in the VCF files
            VCF_CHR_RENAME_BCFTOOLS(chr_vcf_disjoint.to_rename)
            ch_versions = ch_versions.mix(VCF_CHR_RENAME_BCFTOOLS.out.versions)
            ch_vcf_renamed = VCF_CHR_RENAME_BCFTOOLS.out.vcf_renamed
        } else {
            chr_vcf_disjoint.to_rename.map {
                chr_names = it[3].size() > 5 ? it[3][0..5] + ['...'] : it[3]
                error "Contig names: ${chr_names} in VCF: ${it[1]} are not present in reference genome with same writing. Please set `rename_chr` to `true` to rename the contigs."
            }
            chr_bam_disjoint.to_rename.map {
                chr_names = it[3].size() > 5 ? it[3][0..5] + ['...'] : it[3]
                error "Contig names: ${chr_names} in BAM: ${it[1]} are not present in reference genome with same writing. Please set `rename_chr` to `true` to rename the contigs."
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
