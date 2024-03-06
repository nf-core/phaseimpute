include { SAMTOOLS_FAIDX              } from '../../../modules/nf-core/samtools/faidx/main'

workflow GET_REGION {
    take:
        input_region // Region string to use ["all", "chr1", "chr1:0-1000"]
        ch_fasta     // [meta, fasta, fai]
    
    main:
        ch_versions      = Channel.empty()
        ch_multiqc_files = Channel.empty()
        // Gather regions to use and create the meta map
        if (input_region ==~ '^chr[0-9XYM]+$' || input_region == "all") {
            if (ch_fasta[2] == null) {
                SAMTOOLS_FAIDX(ch_fasta[0..1], Channel.of([[],[]]))
                ch_versions      = ch_versions.mix(SAMTOOLS_FAIDX.out.versions.first())
                ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_FAIDX.out.fai.collect{it[1]})
                ch_fasta[2]      = SAMTOOLS_FAIDX.out.fai
            }
            ch_regions = ch_fasta[2]
                .splitCsv(header: ["chr", "size", "offset", "lidebase", "linewidth", "qualoffset"], sep: "\t")
            if (input_region != "all") {
                ch_regions = ch_regions.filter{meta, rows -> rows.chr == input_region}
            }
            ch_regions = ch_regions
                .map{ meta, row -> [meta + ["chr": row.chr], row.chr + ":0-" + row.size]}
                .map{ metaC, region -> [metaC + ["region": region], region]}
        } else {
            if (input_region ==~ '^chr[0-9XYM]+:[0-9]+-[0-9]+$') {
                ch_regions = Channel.from([input_region])
                    .map{ region -> [["region": region], region]}
            } else {
                error "Invalid input_region: ${input_region}"
            }
        }
    emit:
        ch_regions        = ch_regions       // channel: [ meta, region ]
        versions          = ch_versions      // channel: [ versions.yml ]
        multiqc_files     = ch_multiqc_files // channel: [ multiqc_report.html ]
}
