include { GLIMPSE2_CONCORDANCE        } from '../../../modules/nf-core/glimpse2/concordance'
include { GAWK                        } from '../../../modules/nf-core/gawk'
include { ADD_COLUMNS                 } from '../../../modules/local/addcolumns'
include { GUNZIP                      } from '../../../modules/nf-core/gunzip'

workflow VCF_CONCORDANCE_GLIMPSE2 {

    take:
        ch_vcf_emul   // VCF file with imputed genotypes [ [id, panel, tool, chr], vcf, csi]
        ch_vcf_truth  // VCF file with truth genotypes   [ [id, panel, chr], vcf, csi]
        ch_vcf_freq   // VCF file with panel frequencies [ [panel], vcf, csi]
        ch_region     // Regions to process              [ [chr, region], region]

    main:

    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()

    ch_concordance = ch_vcf_emul
        .map{metaIPTC, vcf, csi -> [metaIPTC.subMap("id", "panel"), metaIPTC, vcf, csi]}
        .combine(ch_vcf_truth
            .map{metaIPC, vcf, csi -> [ metaIPC.subMap("id", "panel"), vcf, csi ]}
            , by: 0
        )
        .combine(ch_vcf_freq)
        .combine(ch_region.map{[it[1]]}.collect().toList())
        .map{metaI, metaIPTC, emul, e_csi, truth, t_csi, metaP, freq, f_csi, regions ->
            [metaIPTC, emul, e_csi, truth, t_csi, freq, f_csi, [], regions]
        }

    GLIMPSE2_CONCORDANCE (
        ch_concordance,
        [[], [], params.bins, [], []],
        params.min_val_gl, params.min_val_dp
    )
    ch_versions = ch_versions.mix(GLIMPSE2_CONCORDANCE.out.versions.first())

    ch_multiqc_files = ch_multiqc_files.mix(GLIMPSE2_CONCORDANCE.out.errors_cal.map{meta, txt -> [txt]})
    ch_multiqc_files = ch_multiqc_files.mix(GLIMPSE2_CONCORDANCE.out.errors_grp.map{meta, txt -> [txt]})
    ch_multiqc_files = ch_multiqc_files.mix(GLIMPSE2_CONCORDANCE.out.errors_spl.map{meta, txt -> [txt]})
    ch_multiqc_files = ch_multiqc_files.mix(GLIMPSE2_CONCORDANCE.out.rsquare_grp.map{meta, txt -> [txt]})
    ch_multiqc_files = ch_multiqc_files.mix(GLIMPSE2_CONCORDANCE.out.rsquare_spl.map{meta, txt -> [txt]})
    ch_multiqc_files = ch_multiqc_files.mix(GLIMPSE2_CONCORDANCE.out.rsquare_per_site.map{meta, txt -> [txt]})

    GUNZIP(GLIMPSE2_CONCORDANCE.out.errors_grp)
    ch_versions = ch_versions.mix(GUNZIP.out.versions.first())

    ADD_COLUMNS(GUNZIP.out.gunzip)
    ch_versions = ch_versions.mix(ADD_COLUMNS.out.versions.first())

    GAWK(
        ADD_COLUMNS.out.txt
            .map{meta, txt -> [["id":"TestQuality"], txt]}
            .groupTuple(),
        []
    )
    ch_versions = ch_versions.mix(GAWK.out.versions.first())

    emit:
    stats           = GAWK.out.output             // [ [all], txt ]
    versions        = ch_versions                 // channel: [ versions.yml ]
    multiqc_files   = ch_multiqc_files
}
