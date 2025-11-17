process MAPCONVERT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.3.0' :
        'biocontainers/gawk:5.3.0' }"

    input:
    tuple val(meta), path(map_file)
    val(sep)
    val(header)
    val(col_names)

    output:
    tuple val(meta), path('*.glimpse.map'), emit: glimpse_map
    tuple val(meta), path('*.plink.map')  , emit: plink_map
    tuple val(meta), path('*.stitch.map') , emit: stitch_map
    tuple val(meta), path('*.minimac.map'), emit: minimac_map
    path "versions.yml"                   , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Decompress if necessary (pipe directly)
    MAP_CMD="cat ${map_file}"
    [[ ${map_file} == *.gz ]] && MAP_CMD="gunzip -c ${map_file}"

    SEP="${sep}"
    HDR=${header}
    CHR=${meta.chr}
    COLS=${col_names}

    # Single AWK to extract columns, validate, compute rate, and produce all outputs
    \$MAP_CMD | awk --bignum -v FS="\$SEP" -v chr="\$CHR" -v hdr="\$HDR" -v cols="\$COLS" '
    BEGIN {
        n = split(cols, c, ",")
        for (i=1;i<=n;i++) colname[i] = tolower(c[i])
        # Print headers for tool-specific outputs

        glimpse = "${prefix}.glimpse.map"
        minimac = "${prefix}.minimac.map"
        plink   = "${prefix}.plink.map"
        stitch  = "${prefix}.stitch.map"

        # Print headers with correct delimiters
        printf "pos\\tchr\\tcM\\n" > glimpse
        printf "#chr\\tposition\\tGenetic_Map(cM)\\n" > minimac
        printf "position COMBINED_rate.cM.Mb. Genetic_Map.cM.\\n" > stitch

        data_nr = 0
        init_cm_set = 0
    }

    NR==1 && hdr=="true" { next }

    {
        delete val
        for (i=1; i<=NF && i<=n; i++) val[colname[i]] = \$i

        chr_file = ("chr" in val  ? val["chr"] : chr)
        pos      = ("pos" in val  ? val["pos"] : "NA")
        cm       = ("cm"  in val  ? val["cm"]  : "NA")
        id       = ("id"  in val  ? val["id"]  : ".")

        if (pos=="NA" || cm=="NA" || pos=="" || cm=="") {
            print "Error: Position and cM missing" > "/dev/stderr"; exit 1
        }
        if (chr_file != chr) {
            print "Error: Chromosome mismatch" > "/dev/stderr"; exit 1
        }

        data_nr++

        # first data row: capture initial cM offset
        if (data_nr == 1) {
            prev_pos = pos
            prev_cm  = cm
            prev_id  = id
            init_cm  = cm
            init_cm_set = 1

            # output non-stitch files for the current row
            printf "%s\\t%s\\t%s\\n", pos, chr, cm >> glimpse
            printf "%s\\t%s\\t%s\\n", chr, pos, cm >> minimac
            printf "%s %s %s %s\\n", chr, id, cm, pos >> plink
            next
        }

        # compute forward rate for previous row (interval prev -> current)
        delta_bp = pos - prev_pos
        delta_cm = cm - prev_cm
        rate_prev = (delta_bp > 0 ? (delta_cm / delta_bp * 1e6) : 0)

        # adjusted cM for previous row (subtract init offset)
        adj_prev_cm = prev_cm - init_cm

        # write STITCH row for previous entry (position, rate, adjusted cM)
        printf "%.0f %.15f %.12f\\n", prev_pos, rate_prev, adj_prev_cm >> stitch

        # outputs for current row (non-STITCH)
        printf "%s\\t%s\\t%s\\n", pos, chr, cm >> glimpse
        printf "%s\\t%s\\t%s\\n", chr, pos, cm >> minimac
        printf "%s %s %s %s\\n", chr, id, cm, pos >> plink

        # shift to next
        prev_pos = pos
        prev_cm  = cm
        prev_id  = id
    }

    END {
        if (data_nr == 0) exit 0

        # last row: no forward interval -> rate = 0
        adj_prev_cm = prev_cm - init_cm
        printf "%.0f %.15f %.12f\\n", prev_pos, 0, adj_prev_cm >> stitch
    }'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.glimpse.map
    touch ${prefix}.plink.map
    touch ${prefix}.stitch.map
    touch ${prefix}.minimac.map

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//')
    END_VERSIONS
    """
}
