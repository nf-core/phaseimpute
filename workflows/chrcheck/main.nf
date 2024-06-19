include { VCFCHREXTRACT as VCFCHRBF } from '../../modules/local/vcfchrextract'
include { BAMCHREXTRACT as BAMCHRBF } from '../../modules/local/bamchrextract'
include { BAM_CHR_RENAME_SAMTOOLS   } from '../../subworkflows/local/bam_chr_rename_samtools'
include { VCF_CHR_RENAME_BCFTOOLS   } from '../../subworkflows/local/vcf_chr_rename_bcftools'
include { VCFCHREXTRACT as VCFCHRAF } from '../../modules/local/vcfchrextract'
include { BAMCHREXTRACT as BAMCHRAF } from '../../modules/local/bamchrextract'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CHR_CHECK {
    take:
        ch_input // [[id], file, index]
        ch_fasta // [[id], fasta, fai]

    main:
        ch_versions = Channel.empty()

        // Split the input between VCF and BAM files
        ch_input = ch_input.branch{
            bam: it[1] =~ 'bam|cram|sam'
            vcf: it[1] =~ 'vcf|bcf'
        }
        ch_input.vcf.view()
        ch_input.bam.view()

        // Extract the contig names from the VCF files
        VCFCHRBF(ch_input.vcf)
        ch_versions = ch_versions.mix(VCFCHRBF.out.versions)
        chr_vcf_disjoint = check_chr(VCFCHRBF.out.chr, ch_input.vcf, ch_fasta)

        // Extract the contig names from the BAM files
        BAMCHRBF(ch_input.bam)
        ch_versions = ch_versions.mix(BAMCHRBF.out.versions)
        chr_bam_disjoint = check_chr(BAMCHRBF.out.chr, ch_input.bam, ch_fasta)

        if (params.rename_chr == true) {
            // Rename the contigs in the BAM files
            BAM_CHR_RENAME_SAMTOOLS(
                chr_bam_disjoint.to_rename.map{meta, bam, csi, dis, prefix -> [meta, bam, csi, prefix]}
            )
            ch_versions = ch_versions.mix(BAM_CHR_RENAME_SAMTOOLS.out.versions)

            // Rename the contigs in the VCF files
            VCF_CHR_RENAME_BCFTOOLS(
                chr_vcf_disjoint.to_rename.map{meta, vcf, csi, dis, prefix -> [meta, vcf, csi]}
            )

            // Check if modification has solved the problem
            BAMCHRAF(BAM_CHR_RENAME_SAMTOOLS.out.bam_renamed.map{ meta, bam, csi -> [meta, bam] })
            chr_bam_disjoint_after = check_chr(BAMCHRAF.out.chr, BAM_CHR_RENAME_SAMTOOLS.out.bam_renamed, ch_fasta)

            VCFCHRAF(VCF_CHR_RENAME_BCFTOOLS.out.vcf_renamed.map{ meta, vcf, csi -> [meta, vcf] })
            chr_vcf_disjoint_after = check_chr(VCFCHRAF.out.chr, VCF_CHR_RENAME_BCFTOOLS.out.vcf_renamed, ch_fasta)

            chr_bam_disjoint_after.to_rename.map{
                error "Even after renaming errors are still present. Please check the contigs name: ${it[3]} in BAM: ${it[1]} and fasta file."
            }
            chr_vcf_disjoint_after.to_rename.map{
                error "Even after renaming errors are still present. Please check the contigs name: ${it[3]} in VCF: ${it[1]} and fasta file."
            }
            // If no errors are present, we can use the renamed files
            ch_bam_renamed = BAM_CHR_RENAME_SAMTOOLS.out.bam_renamed
            ch_vcf_renamed = VCF_CHR_RENAME_BCFTOOLS.out.vcf_renamed
        } else {
            chr_vcf_disjoint.to_rename.map {
                error "Contig names: ${it[3]} in the VCF: ${it[1]} are not present in the reference genome. Please set `rename_chr` to `true` to rename the contigs."
            }
            chr_bam_disjoint.to_rename.map {
                error "Contig names: ${it[3]} in the BAM: ${it[1]} are not present in the reference genome. Please set `rename_chr` to `true` to rename the contigs."
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


def check_chr(ch_chr, ch_input, ch_fasta){
    chr_checked = ch_chr
        .combine(ch_input, by:0)
        .combine(ch_fasta)
        .map{metaI, chr, file, index, metaG, fasta, fai ->
            [
                metaI, file, index,
                chr.readLines()*.split(' ').collect{it[0]},
                fai.readLines()*.split('\t').collect{it[0]}
            ]
        }
        .branch{ meta, file, index, chr, fai ->
            no_rename: (chr - fai).size() == 0
                return [meta, file, index]
            to_rename: true
                return [meta, file, index, chr-fai, chr-fai =~ "chr" ? "nochr" : "chr"]
        }
    return chr_checked
}
