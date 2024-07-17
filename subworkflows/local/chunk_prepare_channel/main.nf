workflow CHUNK_PREPARE_CHANNEL {

    take:
    ch_chunks // channel:   [ [id, chr], txt ]
    tool

    main:

    ch_versions      = Channel.empty()

    if(tool == "glimpse"){
        ch_chunks = ch_chunks.map { chr, txt -> [chr, file(txt)]}
                .splitCsv(header: ['ID', 'Chr', 'RegionIn', 'RegionOut', 'Size1', 'Size2'], sep: "\t", skip: 0)
                .map { meta, it -> [meta, it["RegionIn"], it["RegionOut"]]}
    }

    if(tool == "quilt") {
        ch_chunks = ch_chunks.map { chr, txt -> [chr, file(txt)]}
            .splitCsv(header: ['ID', 'Chr', 'RegionIn', 'RegionOut', 'Size1', 'Size2'], sep: "\t", skip: 0)
            .map { meta, it ->
                def startEnd = it["RegionIn"].split(':')[1].split('-')
                [metamap, metamap.chr, startEnd[0], startEnd[1]]
            }
    }



    emit:
    chunks  = ch_chunks // channel:   [ [meta], regionstart, regionend ]

}
