include { VCFCHREXTRACT           } from '../../modules/local/vcfchrextract'
include { BAMCHREXTRACT           } from '../../modules/local/bamchrextract'
include { BAM_CHR_RENAME_SAMTOOLS } from '../../subworkflows/local/bam_chr_rename_samtools'
include { VCF_CHR_RENAME_BCFTOOLS } from '../../subworkflows/local/vcf_chr_rename_bcftools'
include { checkChr                } from './function.nf'
include { diffChr                 } from './function.nf'

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
            bam: it[1] =~ 'bam|cram'
            vcf: it[1] =~ 'vcf|bcf'
            other: it[1].size() > 0
            empty: true
        }

        ch_input.other.map {
            error "File: ${it[1]} is not a VCF, BCFT or BAM, CRAM file."
        }

        // Check if channel is empty
        ch_vcf_split = Channel.empty()
        // Extract the contig names from the VCF files
        VCFCHREXTRACT(ch_input.vcf.map{ meta, file, _index, _chr -> [meta, file] })
        ch_versions = ch_versions.mix(VCFCHREXTRACT.out.versions.first())
        ch_vcf_split = checkChr(VCFCHREXTRACT.out.chr, ch_input.vcf)

        ch_bam_split = Channel.empty()
        // Extract the contig names from the BAM files
        BAMCHREXTRACT(ch_input.bam.map{ meta, file, _index, _chr -> [meta, file] })
        ch_versions = ch_versions.mix(BAMCHREXTRACT.out.versions.first())
        ch_bam_split = checkChr(BAMCHREXTRACT.out.chr, ch_input.bam)

        if (params.rename_chr == true) {
            ch_bam_renamed = Channel.empty()
            // Rename the contigs in the BAM files
            BAM_CHR_RENAME_SAMTOOLS(
                ch_bam_split.to_rename.map{meta, bam, csi, _diff, prefix -> [meta, bam, csi, prefix]}
            )
            ch_versions = ch_versions.mix(BAM_CHR_RENAME_SAMTOOLS.out.versions.first())
            ch_bam_renamed = BAM_CHR_RENAME_SAMTOOLS.out.bam_renamed

            ch_vcf_renamed = Channel.empty()
            // Rename the contigs in the VCF files
            VCF_CHR_RENAME_BCFTOOLS(ch_vcf_split.to_rename)
            ch_versions = ch_versions.mix(VCF_CHR_RENAME_BCFTOOLS.out.versions.first())
            ch_vcf_renamed = VCF_CHR_RENAME_BCFTOOLS.out.vcf_renamed
        } else {
            ch_vcf_split.to_rename.map {
                def chr_names = it[3].size() > params.max_chr_names ? it[3][0..params.max_chr_names - 1] + ['...'] : it[3]
                error "Contig names: ${chr_names} in VCF: ${it[1]} are not present in reference genome with same writing. Please set `rename_chr` to `true` to rename the contigs."
            }
            ch_bam_split.to_rename.map {
                def chr_names = it[3].size() > params.max_chr_names ? it[3][0..params.max_chr_names - 1] + ['...'] : it[3]
                error "Contig names: ${chr_names} in BAM: ${it[1]} are not present in reference genome with same writing. Please set `rename_chr` to `true` to rename the contigs."
            }
            ch_vcf_renamed = Channel.empty()
            ch_bam_renamed = Channel.empty()
        }

        ch_output = ch_bam_split.no_rename
            .mix(ch_vcf_split.no_rename)
            .mix(ch_bam_renamed)
            .mix(ch_vcf_renamed)
    emit:
        output   = ch_output             // [ [id], file, index ]
        versions = ch_versions           // channel: [ versions.yml ]
}
