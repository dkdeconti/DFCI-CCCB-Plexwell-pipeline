##############################################################################
# Workflow Definition
##############################################################################

workflow PlexwellProcessing {
    File test_text = "foo.txt" 
    Array[Array[String]] test_array = read_tsv(test_text)

    scatter (line in test_array) {
        call test_map {
            input:
                colA = line[0],
                colB = line[1]
        }
    }
    call write_foo_map {
        input:
            to_map = test_map.mapping_tsv
    }
    call echo_the_map {
        input:
            the_map = write_foo_map.the_map,
            key = "A"
    }
}

##############################################################################
# Task Definitions
##############################################################################

task test_map {
    String colA
    String colB

    command {
        echo -e "${colA}\t${colB}" > ${colA}.txt
    }
    output {
        File mapping_tsv = "${colA}.txt"
    }
}

task write_foo_map {
    Array[File] to_map

    command {
        cat ${sep=' ' to_map} > output.tsv
    }
    output {
        File output_tsv = "output.tsv"
        Map[String, String] the_map = read_map("output.tsv")
    }
}

task echo_the_map {
    Map[String, String] the_map
    String key
    String val = the_map[key]

    command {
        echo "${val}" > output.txt
    }
    output {
        File output_file = "output.txt"
    }
}