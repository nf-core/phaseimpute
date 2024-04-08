include { QUILT_QUILT     } from '../../../modules/nf-core/quilt/quilt/main'

workflow IMPUTE_QUILT {

    take:
    ch_hap_legend                            // channel: [ val(meta),hap, legend ]
    ch_input                                //  channel: [ val(meta),bam, bai ]
    ch_chunks                              //   channel: [ val(meta), start_coordinate, end_coordinate, number ]


    main:

    ch_versions = Channel.empty()

    posfile = []
    phasefile = []
    posfile_phasefile = [[id: null], posfile, phasefile]
    genetic_map_file = []
    fasta = [[id:'test'], []]
    def ngen = 100
    def buffer = 10000

    ch_bam_bamlist = ch_input

    if (genetic_map_file.isEmpty()) {
        ch_hap_chunks = ch_hap_legend.combine(ch_chunks, by:0).map { it + ngen + buffer + [[]] }
    } else {
        println "genetic_map_file is not empty"
        // Add ngen and buffer
        ch_hap_chunks = ch_hap_legend.join(ch_chunks, by:0).join(genetic_map_file)
    }

    ch_quilt = ch_bam_bamlist.combine(ch_hap_chunks)
    ch_quilt_input = ch_quilt.map { it.take(4) + it.drop(5) }

    // Add metamap with chromosome information
    ch_quilt_input = ch_quilt_input
                        .map{ meta, bam, bai, bamlist, hap, legend, chr, start, end, ngen2, buffer2, genetic ->
                        return [['id': meta.id, 'chr': chr] , bam, bai, bamlist, hap, legend, chr, start, end, ngen2, buffer2, genetic]
                        }

    // Run QUILT
    QUILT_QUILT ( ch_quilt_input, posfile_phasefile, fasta )






    emit:
    ch_imputedvcf              = QUILT_QUILT.out.vcf                    // channel:  [ meta, vcf ]
    versions                   = ch_versions                           // channel:   [ versions.yml ]
}
