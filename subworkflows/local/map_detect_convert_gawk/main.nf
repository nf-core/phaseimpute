include { MAPAUTODETECT           } from '../../../modules/local/mapautodetect'
include { MAPCONVERT              } from '../../../modules/local/mapconvert'

workflow MAP_DETECT_CONVERT_GAWK {
    take:
    ch_map          // channel: [ [id], map ]
    sep             // val: separator used in map
    header          // boolean: header present in map
    colnames        // string: columns order in map

    main:

    // Split the channel into empty and non-empty map files
    ch_map_branched = ch_map
        .branch { meta, map_file ->
            empty: !map_file || map_file.isEmpty()
                return [meta, map_file]
            valid: true
                return [meta, map_file]
        }

    // Convert null values to channels
    ch_sep = sep != null ? channel.value(sep) : channel.value(null)
    ch_header = header != null ? channel.value(header) : channel.value(null)
    ch_colnames = colnames != null ? channel.value(colnames) : channel.value(null)

    MAPAUTODETECT(ch_map_branched.valid, ch_sep, ch_header, ch_colnames)

    MAPCONVERT(ch_map_branched.valid.join(MAPAUTODETECT.out.detected))

    ch_map_glimpse = MAPCONVERT.out.glimpse_map.mix(ch_map_branched.empty)
    ch_map_stitch  = MAPCONVERT.out.stitch_map.mix(ch_map_branched.empty)
    ch_map_plink   = MAPCONVERT.out.plink_map.mix(ch_map_branched.empty)
    ch_map_minimac = MAPCONVERT.out.minimac_map.mix(ch_map_branched.empty)

    emit:
    map_glimpse = ch_map_glimpse
    map_stitch  = ch_map_stitch
    map_plink   = ch_map_plink
    map_minimac = ch_map_minimac
}
