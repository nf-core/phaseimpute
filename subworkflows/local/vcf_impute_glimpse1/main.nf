
include { BAM_GL_BCFTOOLS                    } from '../bam_gl_bcftools'
include { GLIMPSE_PHASE                      } from '../../../modules/nf-core/glimpse/phase'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_1 } from '../../../modules/nf-core/bcftools/index'
include { GLIMPSE_LIGATE                     } from '../../../modules/nf-core/glimpse/ligate'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_2 } from '../../../modules/nf-core/bcftools/index'

workflow VCF_IMPUTE_GLIMPSE1 {

    take:
    ch_input        // channel (mandatory): [ [id], bam, bai ]
    ch_sites_tsv    // channel (mandatory): [ [panel, chr, region], sites, tsv ]
    ch_panel        // channel (mandatory): [ [panel, chr, region], vcf, tbi ]
    ch_chunks       // channel  (optional): [ [chr], region1, region2 ]
    ch_fasta        // channel (mandatory): [ [genome], fa, fai ]

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // Glimpse1 subworkflow
    BAM_GL_BCFTOOLS( // Compute GL for input data once per panel by chromosome
        ch_input,
        ch_sites_tsv,
        ch_fasta
    )
    ch_multiqc_files = ch_multiqc_files.mix(BAM_GL_BCFTOOLS.out.multiqc_files)
    ch_versions = ch_versions.mix(BAM_GL_BCFTOOLS.out.versions)

    samples_file = Channel.of([[]]).collect()
    gmap_file    = Channel.of([[]]).collect()

    ch_phase_input = BAM_GL_BCFTOOLS.out.vcf // [metaIPC, vcf, index]
        .map {metaIPC, vcf, index -> [metaIPC.subMap("panel", "chr"), metaIPC, vcf, index] }
        .combine(ch_panel
            .map{
                metaPC, vcf, index ->
                [["panel": metaPC.id, "chr": metaPC.chr], vcf, index]
            },
            by: 0
        )
        .combine(samples_file)
        .combine(gmap_file)
        .map { metaPC, metaIPC, vcf, index, panel, p_index, sample, gmap ->
            [metaPC.subMap("chr"), metaIPC, vcf, index, panel, p_index, sample, gmap]}
        .combine(ch_chunks
            .map {metaCR, regionin, regionout -> [metaCR.subMap("chr"), metaCR, regionin, regionout]},
            by: 0
        )
        .map{
            metaC, metaIPC, vcf, index, panel, p_index, sample, gmap, metaCR, regionin, regionout
            -> [metaIPC + ["region": regionin], vcf, index, sample, regionin, regionout, panel, p_index, gmap]
        }

    GLIMPSE_PHASE ( ch_phase_input ) // [meta, vcf, index, sample, regionin, regionout, ref, ref_index, map]
    ch_versions = ch_versions.mix(GLIMPSE_PHASE.out.versions )

    BCFTOOLS_INDEX_1 ( GLIMPSE_PHASE.out.phased_variants )
    ch_versions = ch_versions.mix( BCFTOOLS_INDEX_1.out.versions )

    // Ligate all phased files in one and index it
    ligate_input = GLIMPSE_PHASE.out.phased_variants
        .groupTuple()
        .join( BCFTOOLS_INDEX_1.out.csi.groupTuple() )

    GLIMPSE_LIGATE ( ligate_input )
    ch_versions = ch_versions.mix(GLIMPSE_LIGATE.out.versions )

    BCFTOOLS_INDEX_2 ( GLIMPSE_LIGATE.out.merged_variants )
    ch_versions = ch_versions.mix( BCFTOOLS_INDEX_2.out.versions )


    ch_imputed_vcf_tbi = GLIMPSE_LIGATE.out.merged_variants
        .join(BCFTOOLS_INDEX_2.out.csi)
        .map{ metaIPCR, vcf, csi -> [metaIPCR + [tools: "Glimpse1"], vcf, csi] }

    emit:
    vcf_tbi             = ch_imputed_vcf_tbi    // channel: [ [id, chr], vcf, tbi ]
    versions            = ch_versions           // channel: [ versions.yml ]
    multiqc_files       = ch_multiqc_files      // channel: [ multiqc_files.yml ]
}
