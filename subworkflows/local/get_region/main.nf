include { SAMTOOLS_FAIDX              } from '../../../modules/nf-core/samtools/faidx/main'

workflow GET_REGION {
    take:
        input_region // Region string to use ["all", "chr1", "chr1:0-1000"]
        ch_fasta     // [[meta], fasta, fai]

    main:
        ch_versions      = Channel.empty()

        // Gather regions to use and create the meta map
        if (input_region ==~ '^(chr)?[0-9XYM]+$' || input_region == "all") {
            ch_regions = ch_fasta.map{it -> it[2]}
                .splitCsv(header: ["chr", "size", "offset", "lidebase", "linewidth", "qualoffset"], sep: "\t")
                .map{it -> [chr:it.chr, region:"0-"+it.size]}
            if (input_region != "all") {
                ch_regions = ch_regions.filter{it.chr == input_region}
            }
            ch_regions = ch_regions
                .map{ [[chr: it.chr, region: it.chr + ":" + it.region], it.chr + ":" + it.region]}
        } else {
            if (input_region ==~ '^chr[0-9XYM]+:[0-9]+-[0-9]+$') {
                ch_regions = Channel.from([input_region])
                    .map{ [[chr: it.split(":")[0], "region": it], it]}
            } else {
                error "Invalid input_region: ${input_region}"
            }
        }

    emit:
        regions           = ch_regions       // channel: [ meta, region ]
        versions          = ch_versions      // channel: [ versions.yml ]
}
