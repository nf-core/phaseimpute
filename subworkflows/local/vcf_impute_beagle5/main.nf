include { BEAGLE5_BEAGLE                           } from '../../../modules/nf-core/beagle5/beagle'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_BEAGLE  } from '../../../modules/nf-core/bcftools/index'
include { BCFTOOLS_VIEW                            } from '../../../modules/nf-core/bcftools/view'


workflow VCF_IMPUTE_BEAGLE5 {
    take:
    ch_input  // channel: [ [id, chr], vcf, tbi ]
    ch_panel  // channel: [ [id, chr], vcf, tbi ]  
    ch_map    // channel: [ [chr], map]

    main:
    ch_versions = Channel.empty()

    // Branch input files based on format 
    ch_input
        .branch { meta, vcf, tbi ->
            bcf: vcf.toString().contains('.bcf')
            vcf: vcf.toString().contains('.vcf')
        }
        .set { ch_input_branched }

    // Convert BCF to VCF if necessary
    BCFTOOLS_VIEW(
        ch_input_branched.bcf,
        [],
        [],
        []
    )
    ch_versions = ch_versions.mix(BCFTOOLS_VIEW.out.versions.first())

    // Combine VCF files 
    ch_ready_vcf = ch_input_branched.vcf
        .mix(BCFTOOLS_VIEW.out.vcf.join(BCFTOOLS_VIEW.out.csi))

    // Prepare input channels for BEAGLE5 by combining VCF, panel, and map files
    ch_beagle_input = ch_ready_vcf
    .map { meta, vcf, tbi ->
        [ meta.chr, meta, vcf, tbi ]
    }
    .combine(
        ch_panel.map { meta, vcf, idx ->
            [ meta.chr, meta, vcf ]
        },
        by: 0
    )
    .combine(
        ch_map.map { meta, map ->
            [ meta.chr, map ]
        },
        by: 0
    )
    .map { chr, target_meta, vcf, tbi, panel_meta, panel_vcf, map ->
        [ target_meta + { panel: panel_meta.id }, vcf, panel_vcf, map, [], [] ]
    }


    // Run BEAGLE5 imputation
    BEAGLE5_BEAGLE(ch_beagle_input)
    ch_versions = ch_versions.mix(BEAGLE5_BEAGLE.out.versions.first())

    // Index the imputed VCF files
    BCFTOOLS_INDEX_BEAGLE(BEAGLE5_BEAGLE.out.vcf)
    ch_versions = ch_versions.mix(BCFTOOLS_INDEX_BEAGLE.out.versions.first())


    ch_imputed_vcf_tbi = BEAGLE5_BEAGLE.out.vcf
        .join(BCFTOOLS_INDEX_BEAGLE.out.csi)
        .map{ meta, vcf, index -> [meta + [tools: "beagle5"], vcf, index] }

    

    emit:
    vcf_tbi  = ch_imputed_vcf_tbi // channel: [ [id, chr, tools], vcf, tbi ]
    versions = ch_versions // channel: [ versions.yml ]
}