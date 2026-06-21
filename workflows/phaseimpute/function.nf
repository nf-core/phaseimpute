def chunkPrepareChannel(ch_chunks, ch_region, tool) {
    def ch_chunks_branched = ch_chunks
        .map{ metaPC, txt -> [ metaPC.subMap("chr"), metaPC, txt] }
        .combine(
            ch_region
                .map{ metaCR, region -> [ metaCR.subMap("chr"), metaCR, region] },
            by: 0
        )
        .branch{ _metaC, metaPC, txt, _metaCR, region ->
            txt: txt != []
                return [metaPC, txt]
            empty: true
                return [metaPC, region]
        }
    if(tool == "glimpse1"){
        def ch_chunks_txt = ch_chunks_branched.txt.map { metaPC, txt -> [metaPC, file(txt)]}
            .splitCsv(header: ['ID', 'Chr', 'RegionIn', 'RegionOut', 'Size1', 'Size2'], sep: "\t", skip: 0)
            .map { meta, it -> [meta, it["RegionIn"], it["RegionOut"]]}
        def ch_chunks_region = ch_chunks_branched.empty.map{
            metaPC, region -> [ metaPC, region, region ]
        }
        return ch_chunks_txt.mix(ch_chunks_region)
    } else if(tool == "quilt") {
        def ch_chunks_txt = ch_chunks_branched.txt.map { metaC, txt -> [metaC, file(txt)]}
            .splitCsv(header: ['ID', 'Chr', 'RegionIn', 'RegionOut', 'Size1', 'Size2'], sep: "\t", skip: 0)
            .map { metaC, it ->
                def startEnd = it["RegionIn"].split(':')[1].split('-')
                [ metaC, metaC.chr, startEnd[0], startEnd[1] ]
            }
        def ch_chunks_region = ch_chunks_branched.empty.map{
            metaPC, region ->
            def startEnd = region.split(':')[1].split('-')
            [ metaPC, metaPC.chr, startEnd[0], startEnd[1] ]
        }
        return ch_chunks_txt.mix(ch_chunks_region)
    } else {
        error "ERROR: Only 'glimpse1' and 'quilt' output format are supported. Got ${tool}"
    }
}

def chRegionToBed(ch_regions) {
    def ch_bed = ch_regions
        .map{ _meta, region ->
            def chr=region.split(":")[0]
            def pos=region.split(":")[1]
            def start=pos.split("-")[0]
            def end=pos.split("-")[1]
            "${chr}\t${start}\t${end}"
        }
        .collectFile(name: "regions.bed", newLine: true)
        .map{ bed -> [[chr : "all"], bed] }

    return ch_bed
}
