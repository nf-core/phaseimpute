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

    def ch_chunks_in_out = channel.of()

    ch_chunks_in_out = ch_chunks_branched.txt.map { metaPC, txt -> [metaPC, file(txt)]}
        .splitCsv(sep:"\t", skip:0)
        .map{ meta, row ->
            if (row.size() == 6 ) {
                return [ meta, row[2], row[3] ] // header is 'ID', 'Chr', 'RegionIn', 'RegionOut', 'Size1', 'Size2'
            }
            if (row.size() == 8) {
                return [ meta, row[3], row[2] ] // header is 'ID', 'Chr', 'RegionBuff', 'RegionCnk', 'WindowCm', 'WindowMb', 'NbTotVariants', 'NbComVariants'
            }
            error "Chunks csv should have either 6 columns (glimpse V1 format) or 8 columns (glimpse V2 format)"
        }

    if(tool == "glimpse1"){
        def ch_chunks_region = ch_chunks_branched.empty.map{
            metaPC, region -> [ metaPC, region, region ]
        }
        return ch_chunks_in_out.mix(ch_chunks_region)
    } else if(tool == "quilt") {
        return ch_chunks_in_out
            .map{meta, region_in, _region_out -> [meta, region_in] }
            .mix(ch_chunks_branched.empty)
            .map { meta, it ->
                def startEnd = it.split(':')[1].split('-')
                [ meta, meta.chr, startEnd[0], startEnd[1] ]
            }
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
