##############################################################################
# Workflow Definition
##############################################################################

workflow PlexwellProcessing {
    File inputs_tsv
    Array[Array[String]] inputs_array = read_tsv(inputs_tsv)

    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File ref_bwt
    File ref_sa
    File ref_amb
    File ref_ann
    File ref_pac

    scatter (inputs in inputs_array) {
        String sample_name = inputs[0]
        File fastqgz = inputs[1]

        call GetRGValues {
            input:
                sample_name = sample_name,
                fastqgz = fastqgz
        }
    }
    call MapPairedFastq {
        input:
            inputs_tsv = inputs_tsv
    }
    call MapRGFileLocations {
        input:
            rg_loc_to_map = GetRGValues.rg_loc
    }
    scatter (inputs in inputs_array) {
        String sample_name = inputs[0]
        File fastqgz = inputs[1]
        Array[File] rg_values = read_lines(MapRGFileLocations.rg_loc_map[sample_name])
        String rgid = rg_values[0]
        String rgpl = rg_values[1]
        String rglb = rg_values[2]
        String rgsm = rg_values[3]
        String rgcn = rg_values[4]
        String rgpu = rg_values[5]

        call AlignFastq {
            input:
                sample_name = sample_name,
                first_fastq = fastqgz,
                second_fastq = MapPairedFastq.mapped_pairs[sample_name],
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                ref_bwt = ref_bwt,
                ref_sa = ref_sa,
                ref_amb = ref_amb,
                ref_ann = ref_ann,
                ref_pac = ref_pac,
        }
        call AddRGFixAndSort {
            input:
                input_bam = AlignFastq.output_bam,
                sample_name = sample_name,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                ref_dict = ref_dict,
                rgid = rgid,
                rgpl = rgpl,
                rglb = rglb,
                rgsm = rgsm,
                rgcn = rgcn,
                rgpu = rgpu
        }
    }
}

##############################################################################
# Task Definitions
##############################################################################

task GetRGValues {
    String sample_name
    File fastqgz

    command {
        python /run/media/ddeconti/Data/ddeconti/farline/FireCloud/get_rg_values.py ${fastqgz} > ${sample_name}.RG.txt
        echo -e "${sample_name}\t${sample_name}.RG.txt" > ${sample_name}.RG_loc.txt
    }
    output {
        File rg_values = "${sample_name}.RG.txt"
        File rg_loc = "${sample_name}.RG_loc.txt"
    }
}

task MapRGFileLocations {
    Array[File] rg_loc_to_map

    command {
        cat ${sep=' ' rg_loc_to_map} > RG_loc_map.tsv
    }
    output {
        File rg_loc_map_tsv = "RG_loc_map.tsv"
        Map[String, String] rg_loc_map = read_map(rg_loc_map_tsv)
    }
}

task MapPairedFastq {
    File inputs_tsv

    command {
        cut -f1,3 ${inputs_tsv} > mapped_pairs.tsv
    }
    output {
        File mapped_pairs_file = "${mapped_pairs.tsv}"
        Map[String, String] mapped_pairs = read_map(mapped_pairs_file)
    }
}

task AlignFastq {
    File ref_fasta
    File ref_fasta_index
    File ref_bwt
    File ref_sa
    File ref_amb
    File ref_ann
    File ref_pac

    String sample_name
    File first_fastq
    File second_fastq

    command {
        ~/bin/bwa mem ${ref_fasta} ${first_fastq} ${second_fastq} > ${sample_name}.sam
        java -Xmx2500m -jar picard.jar \
            SamFormatConverter \
            I=${sample_name}.sam \
            O=${sample_name}.bam
    }
    output {
        File output_bam = "${sample_name}.bam"
    }
}

task AddRGFixAndSort {
    File input_bam
    String sample_name

    File ref_fasta
    File ref_fasta_index
    File ref_dict

    String rgid
    String rgpl
    String rglb
    String rgsm
    String rgcn
    String rgpu

    command {
        ~/bin/java -jar ~/bin/picard.jar \
            I=${input_bam} \
            O=/dev/stdout \
            RGID=${rgid} \
            RGPL=${rgpl} \
            RGLB=${rglb} \
            RGSM=${rgsm} \
            RGCN=${rgcn} \
            RGPU=${rgpu} |
        ~/bin/java -Xmx4000m -jar ~/bin/picard.jar \
            SortSam \
            INPUT=/dev/stdin \
            OUTPUT=/dev/stdout \
            SORT_ORDER="coordinate" \
            CREATE_INDEX=false \
            CREATE_MD5_FILE=false | \
        ~/bin/java -Xmx100m -jar ~/bin/picard.jar \
            SetNmAndUqTags \
            INPUT=/dev/stdin \
            OUTPUT=${sample_name}.sorted.bam \
            CREATE_INDEX=true \
            CREATE_MD5_FILE=false \
            REFERENCE_SEQUENCE=${ref_fasta};
    }
    output {
        File output_bam = "${sample_name}.sorted.bam"
    }
}