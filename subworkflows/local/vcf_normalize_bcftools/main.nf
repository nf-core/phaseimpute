include { BCFTOOLS_NORM   } from '../../../modules/nf-core/bcftools/norm'
include { BCFTOOLS_VIEW   } from '../../../modules/nf-core/bcftools/view'
include { BCFTOOLS_PLUGINFILLTAGS } from '../../../modules/nf-core/bcftools/pluginfilltags'

workflow VCF_NORMALIZE_BCFTOOLS {
    take:
    ch_vcf_index    // channel: [ [id, chr], vcf, index ]
    ch_fasta        // channel: [ [genome], fasta, [fai, gzi] ]
    normalize       // boolean
    compute_freq    // boolean

    main:

    ch_fasta = ch_fasta.map { meta, fasta, _fai -> [meta, fasta] }

    // Join duplicated biallelic sites into multiallelic records
    if (normalize) {
        BCFTOOLS_NORM(ch_vcf_index, ch_fasta)

        // Join multiallelic VCF and index
        ch_multiallelic_vcf_index = BCFTOOLS_NORM.out.vcf
            .join(BCFTOOLS_NORM.out.index)

        // Remove all multiallelic records and samples specified in the `--remove_samples` command:
        BCFTOOLS_VIEW(ch_multiallelic_vcf_index, [], [], [])

        // Join biallelic VCF and index
        ch_vcf_index = BCFTOOLS_VIEW.out.vcf
            .join(BCFTOOLS_VIEW.out.index)
    }

    // (Optional) Recompute panel AC/AN/AF/NS INFO fields from GT
    if (compute_freq) {
        BCFTOOLS_PLUGINFILLTAGS(ch_vcf_tbi, [], [], [])

        // Join fixed vcf and index
        ch_vcf_tbi = BCFTOOLS_PLUGINFILLTAGS.out.vcf
            .join(BCFTOOLS_PLUGINFILLTAGS.out.index)
    }
    emit:
    vcf_index = ch_vcf_index // channel: [ [id, chr], vcf, [tbi, csi] ]
}
