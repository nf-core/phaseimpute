include { GLIMPSE2_CHUNK                         } from '../../../modules/nf-core/glimpse2/chunk'
include { SHAPEIT5_PHASECOMMON                   } from '../../../modules/nf-core/shapeit5/phasecommon'
include { SHAPEIT5_LIGATE                        } from '../../../modules/nf-core/shapeit5/ligate'
include { BCFTOOLS_INDEX as VCF_BCFTOOLS_INDEX_1 } from '../../../modules/nf-core/bcftools/index'
include { BCFTOOLS_INDEX as VCF_BCFTOOLS_INDEX_2 } from '../../../modules/nf-core/bcftools/index'

workflow VCF_PHASE_SHAPEIT5 {

    take:
    ch_vcf        // channel (mandatory) : [ [panel_id, chr], vcf, index, pedigree ]
    ch_region     // channel (mandatory) : [ [chr, region], region ]
    ch_ref        // channel (optional)  : [ [id, chr], vcf, index ]
    ch_scaffold   // channel (optional)  : [ [id, chr], vcf, index ]
    ch_map        // channel (optional)  : [ [chr], map]
    chunk_model   // channel (mandatory) : [ model ]

    main:

    ch_versions = Channel.empty()

    // Chunk with Glimpse2
    ch_input_glimpse2 = ch_vcf
        .map{
            metaPC, vcf, csi, _pedigree -> [metaPC.subMap("chr"), metaPC, vcf, csi]
        }
        .combine(ch_region.map{ metaCR, region -> [metaCR.subMap("chr"), region]}, by:0)
        .combine(ch_map.map{ metaCR, map -> [metaCR.subMap("chr"), map]}, by: 0)
        .map{
            _metaC, metaIC, vcf, csi, region, gmap -> [metaIC, vcf, csi, region, gmap]
        }
    GLIMPSE2_CHUNK ( ch_input_glimpse2, chunk_model )
    ch_versions = ch_versions.mix( GLIMPSE2_CHUNK.out.versions.first() )

    // Rearrange channels
    ch_chunks_glimpse2 = GLIMPSE2_CHUNK.out.chunk_chr
        .splitCsv(
            header: [
                'ID', 'Chr', 'RegionBuf', 'RegionCnk', 'WindowCm',
                'WindowMb', 'NbTotVariants', 'NbComVariants'
            ], sep: "\t", skip: 0
        )
        .map { metaIC, it -> [metaIC, it["RegionBuf"], it["RegionCnk"]]}

    ch_phase_input = ch_vcf
        .combine(ch_chunks_glimpse2, by:0)
        .map{
            metaIC, vcf, csi, pedigree, regionbuf, regioncnk -> [metaIC.subMap("chr"), metaIC, vcf, csi, pedigree, regionbuf, regioncnk]
        }
        .combine(ch_map.map{ metaCR, map -> [metaCR.subMap("chr"), map]}, by:0)
        .map { _metaC, metaIC, vcf, index, pedigree, regionbuf, regioncnk, gmap ->
            [metaIC + [chunk: regioncnk], vcf, index, pedigree, regionbuf, gmap]
        }

    SHAPEIT5_PHASECOMMON (
        ch_phase_input, ch_ref,
        ch_scaffold
    )
    ch_versions = ch_versions.mix(SHAPEIT5_PHASECOMMON.out.versions.first())

    VCF_BCFTOOLS_INDEX_1(SHAPEIT5_PHASECOMMON.out.phased_variant)
    ch_versions = ch_versions.mix(VCF_BCFTOOLS_INDEX_1.out.versions.first())

    ch_ligate_input = SHAPEIT5_PHASECOMMON.out.phased_variant
        .join(VCF_BCFTOOLS_INDEX_1.out.csi, failOnMismatch:true, failOnDuplicate:true)
        .map{ meta, vcf, csi -> [meta.subMap("panel_id", "chr"), [vcf, meta.chunk], csi]}
        .groupTuple()
        .map{ meta, vcf, csi ->
                [ meta,
                vcf
                    .sort { a, b ->
                        def aStart = a.last().split("-")[-1].toInteger()
                        def bStart = b.last().split("-")[-1].toInteger()
                        aStart <=> bStart
                    }
                    .collect{it.first()},
                csi]}

    SHAPEIT5_LIGATE(ch_ligate_input)
    ch_versions = ch_versions.mix(SHAPEIT5_LIGATE.out.versions.first())

    VCF_BCFTOOLS_INDEX_2(SHAPEIT5_LIGATE.out.merged_variants)
    ch_versions = ch_versions.mix(VCF_BCFTOOLS_INDEX_2.out.versions.first())

    ch_vcf_tbi_join = SHAPEIT5_LIGATE.out.merged_variants
        .join(VCF_BCFTOOLS_INDEX_2.out.csi)

    emit:
    vcf_tbi             = ch_vcf_tbi_join         // channel: [ [id, chr], vcf, csi ]
    versions            = ch_versions             // channel: [ versions.yml ]
}
