include { GLIMPSE2_CONCORDANCE        } from '../../../modules/nf-core/glimpse2/concordance'
include { GAWK                        } from '../../../modules/nf-core/gawk'
include { ADDCOLUMNS                 } from '../../../modules/local/addcolumns'
include { GUNZIP                      } from '../../../modules/nf-core/gunzip'
include { GAWK as GAWK_ERROR_SPL      } from '../../../modules/nf-core/gawk'
include { GAWK as GAWK_RSQUARE_SPL    } from '../../../modules/nf-core/gawk'

workflow VCF_CONCORDANCE_GLIMPSE2 {

    take:
        ch_vcf_emul   // VCF file with imputed genotypes [ [id, panel, tool], vcf, csi]
        ch_vcf_truth  // VCF file with truth genotypes   [ [id, panel], vcf, csi]
        ch_vcf_freq   // VCF file with panel frequencies [ [panel, chr], vcf, csi]
        ch_region     // Regions to process              [ [chr, region], region]

    main:

    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()

    ch_concordance = ch_vcf_emul
        .map{metaIPTC, vcf, csi -> [metaIPTC.subMap("id"), metaIPTC, vcf, csi]}
        .combine(ch_vcf_truth
            .map{metaIPC, vcf, csi -> [ metaIPC.subMap("id"), vcf, csi ]}
            , by: 0
        )
        .combine(ch_vcf_freq)
        .combine(ch_region.map{[it[1]]}.collect().toList())
        .map{ _metaI, metaIPTC, emul, e_csi, truth, t_csi, _metaP, freq, f_csi, regions ->
            [metaIPTC, emul, e_csi, truth, t_csi, freq, f_csi, [], regions]
        }

    GLIMPSE2_CONCORDANCE (
        ch_concordance,
        [[], [], params.bins, [], [], params.min_val_gl, params.min_val_dp]
    )
    ch_versions = ch_versions.mix(GLIMPSE2_CONCORDANCE.out.versions.first())

    GAWK_ERROR_SPL(
        GLIMPSE2_CONCORDANCE.out.errors_spl,
        [],
        false
    )
    ch_versions = ch_versions.mix(GAWK_ERROR_SPL.out.versions.first())

    GAWK_RSQUARE_SPL(
        GLIMPSE2_CONCORDANCE.out.rsquare_spl,
        [],
        false
    )
    ch_versions = ch_versions.mix(GAWK_RSQUARE_SPL.out.versions.first())

    ch_multiqc_files = ch_multiqc_files.mix(GLIMPSE2_CONCORDANCE.out.errors_cal.map{ _meta, txt -> [txt]})
    ch_multiqc_files = ch_multiqc_files.mix(GLIMPSE2_CONCORDANCE.out.errors_grp.map{ _meta, txt -> [txt]})
    ch_multiqc_files = ch_multiqc_files.mix(GAWK_ERROR_SPL.out.output.map{ _meta, txt -> [txt]})
    ch_multiqc_files = ch_multiqc_files.mix(GLIMPSE2_CONCORDANCE.out.rsquare_grp.map{ _meta, txt -> [txt]})
    ch_multiqc_files = ch_multiqc_files.mix(GAWK_RSQUARE_SPL.out.output.map{ _meta, txt -> [txt]})
    ch_multiqc_files = ch_multiqc_files.mix(GLIMPSE2_CONCORDANCE.out.rsquare_per_site.map{ _meta, txt -> [txt]})

    GUNZIP(GLIMPSE2_CONCORDANCE.out.errors_grp)
    ch_versions = ch_versions.mix(GUNZIP.out.versions.first())

    ADDCOLUMNS(GUNZIP.out.gunzip)
    ch_versions = ch_versions.mix(ADDCOLUMNS.out.versions.first())

    GAWK(
        ADDCOLUMNS.out.txt
            .toSortedList { a, b -> a[0].id <=> b[0].id }
            .map { sorted_list ->
                def all_files = sorted_list.collect { it[1] }
                [["id": "AllSamples"], all_files]
            },
        [],
        false
    )
    ch_versions = ch_versions.mix(GAWK.out.versions.first())

    emit:
    stats           = GAWK.out.output             // [ [all], txt ]
    versions        = ch_versions                 // channel: [ versions.yml ]
    multiqc_files   = ch_multiqc_files
}
