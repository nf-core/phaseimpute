process MAPAUTODETECT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.3.0' :
        'biocontainers/gawk:5.3.0' }"

    input:
    tuple val(meta), path(map_file), val(sep), val(header), val(colnames)

    output:
    tuple val(meta), env("DETECTED_SEP"), env("DETECTED_HEADER"), env("DETECTED_COLS"), emit: detected
    tuple val("${task.process}"), val('gawk'), eval("awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//'"), topic: versions, emit: versions_gawk

    when:
        task.ext.when == null || task.ext.when

    script:
    def auto_detect_sep = (sep == null || sep == "" || sep == [])
    def auto_detect_header = (header == null || header == "" || header == [])
    def auto_detect_colnames = (colnames == null || colnames == "" || colnames == [])

    """
    # Decompress if necessary (pipe directly)
    if [[ ${map_file} == *.gz ]]; then
        MAP_CMD="gunzip -c ${map_file}"
    else
        MAP_CMD="cat ${map_file}"
    fi

    # Get first line
    if [[ ${map_file} == *.gz ]]; then
        FIRST_LINE=\$(gunzip -c "${map_file}" 2>/dev/null | head -1 || true)
    else
        FIRST_LINE=\$(head -1 "${map_file}" || true)
    fi

    # Check if first line is empty
    if [ -z "\$FIRST_LINE" ]; then
        echo "Error: First line is empty" >&2
        exit 1
    fi

    # Auto-detect separator
    if [ "${auto_detect_sep}" = "true" ]; then
        TAB_COUNT=\$(echo "\$FIRST_LINE" | tr -cd '\t' | wc -c)
        COMMA_COUNT=\$(echo "\$FIRST_LINE" | tr -cd ',' | wc -c)
        SEMICOLON_COUNT=\$(echo "\$FIRST_LINE" | tr -cd ';' | wc -c)
        SPACE_COUNT=\$(echo "\$FIRST_LINE" | tr -cd ' ' | wc -c)

        if [ \$TAB_COUNT -gt \$COMMA_COUNT ] && [ \$TAB_COUNT -gt \$SEMICOLON_COUNT ] && [ \$TAB_COUNT -gt \$SPACE_COUNT ]; then
            DETECTED_SEP="\t"
        elif [ \$COMMA_COUNT -gt \$SEMICOLON_COUNT ] && [ \$COMMA_COUNT -gt \$SPACE_COUNT ]; then
            DETECTED_SEP=","
        elif [ \$SEMICOLON_COUNT -gt \$SPACE_COUNT ]; then
            DETECTED_SEP=";"
        else
            DETECTED_SEP=" "
        fi
    else
        DETECTED_SEP="${sep}"
    fi

    # Validate separator was detected
    if [ -z "\$DETECTED_SEP" ]; then
        echo "Error: Failed to detect separator" >&2
        exit 1
    fi

    # Auto-detect header
    if [ "${auto_detect_header}" = "true" ]; then
        HEADER_KEYWORDS="pos|cm|snp|position|COMBINED_rate|Genetic_Map|rate"
        if echo "\$FIRST_LINE" | grep -qiE "\$HEADER_KEYWORDS"; then
            DETECTED_HEADER="true"
        else
            DETECTED_HEADER="false"
        fi
    else
        DETECTED_HEADER="${header}"
    fi

    # Auto-detect column names
    if [ "${auto_detect_colnames}" = "true" ]; then
        if [ "\$DETECTED_HEADER" = "true" ]; then
            # Map header names to standard column names
            DETECTED_COLS=\$(echo "\$FIRST_LINE" | awk -F"\$DETECTED_SEP" '{
                for(i=1; i<=NF; i++) {
                    col = tolower(\$i)
                    if (col ~ /position/) {
                        name = "pos"
                    }
                    else if (col ~ /combined_rate/ || col ~ /cm\\.mb/) {
                        name = "rate"
                    }
                    else if (col ~ /genetic_map/ || col ~ /^cm\$/) {
                        name = "cm"
                    }
                    # Keep original name for others
                    else {
                        name = col
                    }
                    printf (i==1 ? "%s" : ",%s"), name
                }
            }')
        else
            NUM_COLS=\$(echo "\$FIRST_LINE" | awk -F"\$DETECTED_SEP" '{print NF}')
            if [ "\$NUM_COLS" -eq 3 ]; then
                DETECTED_COLS="chr,pos,cm"
            elif [ "\$NUM_COLS" -eq 4 ]; then
                DETECTED_COLS="chr,id,cm,pos"
            else
                echo "Error: Cannot auto-detect column names for \$NUM_COLS columns without header" >&2
                exit 1
            fi
        fi
    else
        DETECTED_COLS="${colnames}"
    fi

    echo "Detected: SEP=\$DETECTED_SEP, HEADER=\$DETECTED_HEADER, COLS=\$DETECTED_COLS"
    """

    stub:
    def auto_detect_sep = (sep == null || sep == "" || sep == [])
    def auto_detect_header = (header == null || header == "" || header == [])
    def auto_detect_colnames = (colnames == null || colnames == "" || colnames == [])
    """
    if [ "${auto_detect_sep}" = "true" ]; then
        DETECTED_SEP=${sep}
    else
        DETECTED_SEP=","
    fi

    if [ "${auto_detect_header}" = "true" ]; then
        DETECTED_HEADER=${header}
    else
        DETECTED_HEADER=true
    fi

    if [ "${auto_detect_colnames}" = "true" ]; then
        DETECTED_COLS="${colnames}"
    else
        DETECTED_COLS="pos,chr,id"
    fi
    """
}
