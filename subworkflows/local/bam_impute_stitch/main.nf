include { GAWK                 } from '../../../modules/nf-core/gawk'
include { STITCH               } from '../../../modules/nf-core/stitch'
include { BCFTOOLS_INDEX       } from '../../../modules/nf-core/bcftools/index'

workflow BAM_IMPUTE_STITCH {

    take:
    ch_input        // channel:   [ [id], [bam], [bai], bampaths, bamnames ]
    ch_posfile      // channel:   [ [panel, chr], legend ]
    ch_fasta        // channel:   [ [genome], fa, fai ]
    ch_map          // channel:   [ [chr], map ]
    seed            // value:     seed for random number generator

    main:

    ch_versions      = Channel.empty()

    // Value channels
    def input_empty         = [[]]
    def rdata_empty         = [[]]
    k_val_params            = params.k_val
    ngen_params             = params.ngen

    // Transform posfile to TSV with ','
    GAWK(ch_posfile, [], false)
    ch_versions = ch_versions.mix(GAWK.out.versions.first())

    // Make final channel with parameters
    ch_parameters = GAWK.out.output
        .map{metaPC, posfile -> [[chr: metaPC.chr], metaPC, posfile]}
        .map { it + input_empty + rdata_empty + k_val_params + ngen_params}
        .combine(ch_map, by: 0)
        .map { _metaC, metaPC, posfile, input, rdata, k_val, ngen, map ->
            [metaPC, posfile, input, map, rdata, metaPC.chr, k_val, ngen]
        }

    ch_bam_params = ch_input // Add chr to meta map
        .combine(ch_parameters)
        .map{
            metaI, bam, bai, bampath, bamname, metaPC, posfile, input, map, rdata, chr, k_val, ngen ->
            [
                metaI + metaPC,
                bam, bai, bampath, bamname, posfile, input, map, rdata, chr, k_val, ngen
            ]
        }

    STITCH( ch_bam_params, ch_fasta, seed )
    ch_versions = ch_versions.mix(STITCH.out.versions.first())

    // Index imputed annotated VCF
    BCFTOOLS_INDEX(STITCH.out.vcf)
    ch_versions = ch_versions.mix(BCFTOOLS_INDEX.out.versions.first())

    // Join VCFs and TBIs
    ch_vcf_tbi = STITCH.out.vcf
        .join(BCFTOOLS_INDEX.out.tbi)
        .map { metaIPC, vcf, tbi -> [ metaIPC + [tools: "stitch"], vcf, tbi ] }

    emit:
    vcf_tbi  = ch_vcf_tbi                        // channel:   [ [id, chr], vcf, tbi ]
    versions = ch_versions                       // channel:   [ versions.yml ]

}
