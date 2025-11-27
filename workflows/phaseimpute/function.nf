def chunkPrepareChannel(ch_chunks, tool) {
    if(tool == "glimpse"){
        return ch_chunks.map { chr, txt -> [chr, file(txt)]}
            .splitCsv(header: ['ID', 'Chr', 'RegionIn', 'RegionOut', 'Size1', 'Size2'], sep: "\t", skip: 0)
            .map { meta, it -> [meta, it["RegionIn"], it["RegionOut"]]}
    } else if(tool == "quilt") {
        def ch_chunks_branched = ch_chunks.branch{ meta, txt ->
            txt: txt != []
                return [meta, txt]
            empty: true
                return [meta, txt]
        }
        def ch_chunks_txt = ch_chunks_branched.txt.map { metaC, txt -> [metaC, file(txt)]}
            .splitCsv(header: ['ID', 'Chr', 'RegionIn', 'RegionOut', 'Size1', 'Size2'], sep: "\t", skip: 0)
            .map { metaC, it ->
                def startEnd = it["RegionIn"].split(':')[1].split('-')
                [ metaC, metaC.chr, startEnd[0], startEnd[1] ]
            }
        def ch_chunks_empty = ch_chunks_branched.empty.map { metaC, _empty -> [
            metaC, metaC.chr, [], []
        ]}
        return ch_chunks_txt.mix(ch_chunks_empty)
    } else {
        error "Only 'glimpse' and 'quilt' output format are supported. Got ${tool}"
    }
}
