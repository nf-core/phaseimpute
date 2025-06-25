include { BEAGLE5_BEAGLE                           } from '../../../modules/nf-core/beagle5/beagle'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_BEAGLE } from '../../../modules/nf-core/bcftools/index'
include { BCFTOOLS_VIEW                           } from '../../../modules/nf-core/bcftools/view'
include { BCFTOOLS_CONCAT                         } from '../../../modules/nf-core/bcftools/concat'

workflow VCF_IMPUTE_BEAGLE5 {
    take:
    ch_input  // channel: [ [id, chr], vcf, tbi ]
    ch_panel  // channel: [ [id, chr], vcf, tbi ]  
    ch_map    // channel: [ [chr], map]

    main:
    ch_versions = Channel.empty()

    // Garder les métadonnées du chromosome et normaliser
    ch_input_clean = ch_input
        .map { meta, vcf, tbi -> 
            def chr = meta.chr
            
            // Si pas de chr dans meta, essayer d'extraire du contexte
            if (!chr || chr == meta.id) {
                // Extraire depuis le nom du fichier VCF en fallback
                def vcf_name = vcf.getName().toLowerCase()
                if (vcf_name.contains('chr22') || vcf_name.contains('.22.')) {
                    chr = "22"
                } else if (vcf_name.contains('chr21') || vcf_name.contains('.21.')) {
                    chr = "21"
                } else {
                    // Fallback par défaut
                    chr = "22"
                }
            }
            
            [[id: meta.id, chr: chr], vcf, tbi]
        }

    // Convertir BCF en VCF si nécessaire
    BCFTOOLS_VIEW(
        ch_input_clean,
        [],
        [],
        []
    )
    ch_versions = ch_versions.mix(BCFTOOLS_VIEW.out.versions.first())

    // Joindre VCF et TBI
    ch_converted_vcf_tbi = BCFTOOLS_VIEW.out.vcf
        .join(BCFTOOLS_VIEW.out.tbi)

    // Normaliser les chromosomes pour le matching
    def normalize_chr = { chr ->
        chr.toString().replaceAll(/^chr/, "")
    }

    // Associer chaque VCF avec son panel et map correspondant par chromosome
    ch_beagle_input = ch_converted_vcf_tbi
        .map { meta, vcf, tbi -> 
            [normalize_chr(meta.chr), meta, vcf, tbi]
        }
        .combine(
            ch_panel.map { meta, vcf, idx -> 
                [normalize_chr(meta.chr), vcf]
            },
            by: 0
        )
        .combine(
            ch_map.map { meta, map -> 
                [normalize_chr(meta.chr), map]
            },
            by: 0
        )
        .map { chr, meta, vcf, tbi, panel_vcf, map_file ->
            [meta, vcf, panel_vcf, map_file]
        }

    // Exécuter BEAGLE5 pour chaque combinaison
    ch_beagle_input
        .multiMap { meta, vcf, panel, map ->
            input: [meta, vcf]
            panel: panel
            map: map ?: []
        }
        .set { ch_beagle_channels }

    BEAGLE5_BEAGLE(
        ch_beagle_channels.input,
        ch_beagle_channels.panel,
        ch_beagle_channels.map,
        [],
        []
    )
    ch_versions = ch_versions.mix(BEAGLE5_BEAGLE.out.versions.first())

    // Indexer les résultats
    BCFTOOLS_INDEX_BEAGLE(BEAGLE5_BEAGLE.out.vcf)
    ch_versions = ch_versions.mix(BCFTOOLS_INDEX_BEAGLE.out.versions.first())

    ch_imputed_vcf_tbi = BEAGLE5_BEAGLE.out.vcf
        .join(BCFTOOLS_INDEX_BEAGLE.out.csi)
        .map{ meta, vcf, index -> 
            [meta + [tools: "beagle5"], vcf, index]
        }

    emit:
    vcf_tbi  = ch_imputed_vcf_tbi
    versions = ch_versions
}