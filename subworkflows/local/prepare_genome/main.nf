//
// Prepare reference genome files
//

include { SAMTOOLS_FAIDX } from '../../../modules/nf-core/samtools/faidx'

workflow PREPARE_GENOME {
    take:
    genome         // genome name
    fasta_path     // path to fasta
    fasta_fai_path // path to fasta_fai
    fasta_gzi_path // path to fasta_gzi

    main:

    def is_compressed = fasta_path.toString().endsWith('.gz') || fasta_path.toString().endsWith('.bgz')
    ch_fasta  = channel.of([
        [genome: genome],
        file(fasta_path, checkIfExists:true),
        []
    ])

    if (fasta_fai_path) {
        ch_fai = channel.of(file(fasta_fai_path, checkIfExists:true))
    } else {
        SAMTOOLS_FAIDX(ch_fasta, false)
        ch_fai = SAMTOOLS_FAIDX.out.fai.map{ _meta, fasta_fai -> fasta_fai }
    }
    if (is_compressed) {
        if (fasta_gzi_path) {
            ch_gzi = channel.of(file(fasta_gzi_path, checkIfExists:true))
        } else if (!fasta_fai_path) {
            ch_gzi = SAMTOOLS_FAIDX.out.gzi.map{ _meta, gzi -> gzi }
        } else {
            SAMTOOLS_FAIDX(ch_fasta, false)
            ch_gzi = SAMTOOLS_FAIDX.out.gzi.map{ _meta, gzi -> gzi }
        }
    } else {
        ch_gzi = channel.of([])
    }

    ch_fasta_fai_gzi = ch_fasta
        .map{ meta, fasta, _fai -> [meta, fasta] }
        .combine(ch_fai)
        .combine(ch_gzi)
        .collect()

    emit:
    ch_fasta_fai_gzi
}
