include { GAWK                 } from '../../../modules/nf-core/gawk'
include { STITCH               } from '../../../modules/nf-core/stitch'
include { BCFTOOLS_INDEX       } from '../../../modules/nf-core/bcftools/index'

workflow BAM_IMPUTE_STITCH {

    take:
    ch_input        // channel:   [ [id], bam, bai, bamlist ]
    ch_posfile      // channel:   [ [panel, chr], legend ]
    ch_region       // channel:   [ [chr, region], region ]
    ch_fasta        // channel:   [ [genome], fa, fai ]

    main:

    ch_versions      = Channel.empty()
    // Run STITCH
    seed = params.seed

    // Value channels
    def input_empty         = [[]]
    def rdata_empty         = [[]]
    k_val                   = params.k_val
    ngen                    = params.ngen

    // Transform posfile to TSV with ','
    GAWK(ch_posfile, [])
    ch_versions = ch_versions.mix(GAWK.out.versions)

    // Get chromosomes of posfile
    ch_posfile = GAWK.out.output
        .map{metaPC, posfile -> [[chr: metaPC.chr], metaPC, posfile]}

    // Get chromosomes of fasta
    ch_chromosomes = ch_region
        .map{metaCR, region -> [[chr: metaCR.chr], metaCR.chr]}

    // Make final channel with parameters
    ch_parameters = ch_posfile
        .map { it + input_empty + rdata_empty}
        .join(ch_chromosomes)
        .map { it + k_val + ngen}
        .map { metaC, metaPC, posfile, input, rdata, chr, k_val, ngen ->
            [metaPC, posfile, input, rdata, chr, k_val, ngen]
        }

    ch_bam_params = ch_input // Add chr to meta map
        .combine(ch_parameters)
        .map{
            metaI, bam, bai, bamlist, metaPC, posfile, input, rdata, chr, k_val, ngen ->
            [
                metaI + [chr: metaPC.chr, panel:metaPC.id],
                bam, bai, bamlist, posfile, input, rdata, chr, k_val, ngen
            ]
        }

    STITCH( ch_bam_params, ch_fasta, seed )
    ch_versions = ch_versions.mix(STITCH.out.versions)

    // Index imputed annotated VCF
    BCFTOOLS_INDEX(STITCH.out.vcf)
    ch_versions = ch_versions.mix(BCFTOOLS_INDEX.out.versions)

    // Join VCFs and TBIs
    ch_vcf_tbi = STITCH.out.vcf
        .join(BCFTOOLS_INDEX.out.tbi)
        .map { metaI, vcf, tbi -> [ metaI + [tools: "stitch"], vcf, tbi ] }

    emit:
    vcf_tbi  = ch_vcf_tbi                        // channel:   [ [id, chr], vcf, tbi ]
    versions = ch_versions                       // channel:   [ versions.yml ]

}
