include { IMPUTE5_CHUNK      } from '../../../modules/nf-core/impute5/chunk/main'
include { IMPUTE5_CONVERTREF } from '../../../modules/nf-core/impute5/convertref/main'
include { IMPUTE5_IMPUTE     } from '../../../modules/nf-core/impute5/impute/main'

workflow VCF_IMPUTE_IMPUTE5 {
    take:
    ch_target_data    // tuple val(meta), path(target_bcf), path(target_index), val(region)
    ch_ref_data       // tuple val(meta), path(ref_bcf), path(ref_index), val(region)
    ch_genetic_map    // tuple val(meta), path(genetic_map), val(region)

    main:
    ch_versions = Channel.empty()

    // Étape 1: Convertir la référence au format XCF
    IMPUTE5_CONVERTREF (
        ch_ref_data
    )
    ch_versions = ch_versions.mix(IMPUTE5_CONVERTREF.out.versions)

    // Étape 2: Créer les chunks pour l'imputation
    // Utiliser la région du target pour la cohérence
    ch_chunk_input = ch_target_data
        .join(IMPUTE5_CONVERTREF.out.xcf_file, by: 0)
        .map { meta, target_bcf, target_index, target_region, ref_xcf, ref_index, ref_bin, ref_fam ->
            // Réorganiser pour correspondre exactement à l'input du module
            tuple(meta, ref_xcf, ref_index, ref_bin, ref_fam, target_bcf, target_index, target_region)
        }

    IMPUTE5_CHUNK (
        ch_chunk_input
    )
    ch_versions = ch_versions.mix(IMPUTE5_CHUNK.out.versions)

    // Étape 3: Préparer les données pour l'imputation
    // CORRECTION MAJEURE: Utiliser la même région partout et simplifier la logique
    ch_impute_input = ch_target_data
        .join(IMPUTE5_CONVERTREF.out.xcf_file, by: 0)
        .join(ch_genetic_map, by: 0)
        .map { meta, target_bcf, target_index, target_region, ref_xcf, ref_index, ref_bin, ref_fam, gen_map, map_region ->
            // Utiliser target_region pour l'imputation
            tuple(meta, target_bcf, target_index, ref_xcf, ref_index, ref_bin, ref_fam, gen_map, target_region)
        }

    // Étape 4: Effectuer l'imputation
    IMPUTE5_IMPUTE (
        ch_impute_input
    )
    ch_versions = ch_versions.mix(IMPUTE5_IMPUTE.out.versions)

    emit:
    imputed_bcf = IMPUTE5_IMPUTE.out.bcf      // tuple val(meta), path("*_imputed.bcf")
    imputed_log = IMPUTE5_IMPUTE.out.log      // tuple val(meta), path("*_imputed.log")
    chunks      = IMPUTE5_CHUNK.out.chunks    // tuple val(meta), path("chunks_*.txt")
    xcf_files   = IMPUTE5_CONVERTREF.out.xcf_file // référence convertie
    versions    = ch_versions                 // versions de tous les outils
}
