include { BEAGLE5_BEAGLE } from '../../../modules/nf-core/beagle5/beagle'
include { BCFTOOLS_INDEX } from '../../../modules/nf-core/bcftools/index'
include { BCFTOOLS_VIEW  } from '../../../modules/nf-core/bcftools/view'


workflow VCF_IMPUTE_BEAGLE5 {
    take:
    ch_input  // channel: [ [id], vcf, tbi ]
    ch_panel  // channel: [ [id, chr], vcf, tbi ]
    ch_map    // channel: [ [chr], map]

    main:
    ch_versions = channel.empty()

    // Branch input files based on format
    ch_input
        .branch { _meta, vcf, _tbi ->
            bcf: vcf.toString().contains('.bcf')
            vcf: vcf.toString().contains('.vcf')
        }
        .set { ch_input_branched }

    // Convert BCF to VCF if necessary
    BCFTOOLS_VIEW(
        ch_input_branched.bcf,
        [], [], []
    )
    ch_versions = ch_versions.mix(BCFTOOLS_VIEW.out.versions.first())

    // Combine VCF files
    ch_ready_vcf = ch_input_branched.vcf
        .mix(BCFTOOLS_VIEW.out.vcf.join(BCFTOOLS_VIEW.out.csi))

    // Prepare input channels for BEAGLE5 by combining VCF, panel, and map files
    ch_beagle_input = ch_ready_vcf
        .combine(ch_panel)
        .map {
            metaI, input_vcf, input_index, metaPC, panel_vcf, panel_index -> [
            metaPC.subMap("chr"), metaI + [ panel: metaPC.id, chr: metaPC.chr ],
            input_vcf, input_index, panel_vcf, panel_index
        ]}
        .combine(
            ch_map.map { metaC, map -> [ metaC.subMap("chr"), map ]},
            by: 0
        )
        .map { _metaC, metaIPC, input_vcf, input_index, panel_vcf, panel_index, map ->
            [ metaIPC, input_vcf, input_index, panel_vcf, panel_index, map, [], [] ]
        }.view()


    // Run BEAGLE5 imputation
    BEAGLE5_BEAGLE(ch_beagle_input)
    ch_versions = ch_versions.mix(BEAGLE5_BEAGLE.out.versions.first())

    // Index the imputed VCF files
    BCFTOOLS_INDEX(BEAGLE5_BEAGLE.out.vcf)
    ch_versions = ch_versions.mix(BCFTOOLS_INDEX.out.versions.first())


    ch_vcf_index = BEAGLE5_BEAGLE.out.vcf
        .join(
            BCFTOOLS_INDEX.out.tbi
                .mix(BCFTOOLS_INDEX.out.csi)
        )
        .map{ meta, vcf, index -> [meta + [tools: "beagle5"], vcf, index] }

    emit:
    vcf_index  = ch_vcf_index // channel: [ [id, chr, tools], vcf, index ]
    versions   = ch_versions  // channel: [ versions.yml ]
}
