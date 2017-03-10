##############################################################################
# Workflow Definition
##############################################################################

workflow PlexwellProcessing {
    File inputs_tsv
    Array[Array[String]] inputs_array = read.tsv(inputs_tsv)

    scatter (inputs for inputs_array) {
        File bam = inputs[0]
        File fastqgz = inputs[1]

        call GetRGValues {
            input:
                bam = bam,
                fastqgz = fastqgz
        }
    }
    call MapRGFileLocations {
        input:
            rg_loc_to_map = GetRGValues.rg_loc
    }
    scatter {}
}

##############################################################################
# Task Definitions
##############################################################################

task GetRGValues {
    File bam
    File fastqgz

    command {
        python get_rg_values.py ${fastqgz} > ${bam}.RG.txt
        echo -e "${bam}\t${bam}.RG.txt" > ${bam}.RG_loc.txt
    }
    output {
        File rg_values = "${bam}.RG.txt"
        File rg_loc = "${bam}.RG_loc.txt"
    }
}

task MapRGFileLocations {
    Array[File] rg_loc_to_map

    command {
        cat ${sep=' ' rg_loc_to_map} > RG_loc_map.tsv
    }
    output {
        File rg_loc_map_tsv = "RG_loc_map.tsv"
        Map[String, String] rg_loc_map = read_map("RG_loc_map.tsv")
    }
}