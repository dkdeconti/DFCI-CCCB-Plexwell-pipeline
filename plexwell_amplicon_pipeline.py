#!/usr/bin/python
'''
Main entry for plexwell amplicon analysis.
'''

from collections import defaultdict
import argparse
import gzip
import itertools
import os
import jinja2
import re
import subprocess
import sys


def align_reads(dir_map, fastq_tsv, ref_genome, bin_map, dry=False):
    '''
    Aligns reads from fastqs.
    '''
    fastq_map = map_fastq_from_tsv(fastq_tsv)
    samplename_map = dict((v, k) for k in fastq_map for v in fastq_map[k])
    fastqs = itertools.chain.from_iterable(fastq_map.values())
    mapped_fastqs = map_fastq_pairs(fastqs)
    bams = []
    for first, second in mapped_fastqs:
        sample_name = samplename_map[first]
        basename = re.sub(r'_R[1-2]_[0-9]+\.fastq.gz', '', first)
        basename = '/'.join([dir_map["indbamdir"], basename.split('/')[-1]])
        try:
            bwa = bin_map["bwa"]
            samtools = bin_map["samtools"]
        except KeyError as err:
            sys.stderr.write(err + "\n")
            sys.exit()
        bam = '.'.join([basename, "bam"])
        rg_vals = get_rg_values(sample_name, first)
        cmd = ' '.join([bwa, 'mem -R \"%s\"' % rg_vals,
                        ref_genome, first, second, '|',
                        samtools, "view -bht", ref_genome, "|",
                        samtools, "sort", ">", bam])
        if dry:
            sys.stdout.write(cmd + '\n')
        else:
            subprocess.call(cmd, shell=True)
        bams.append(bam)
    return map_bams_to_sample(samplename_map, bams, dir_map["indbamdir"])


def call_variants_from_bams(bams, ref_genome, dir_map, bin_map, dry=False):
    '''
    Pipes samtools mpileup of bam to varscan.
    '''
    java = bin_map["java"]
    samtools = bin_map["samtools"]
    varscan = bin_map["varscan"]
    snp_vcfs = []
    indel_vcfs = []
    for bam in bams:
        snps = '/'.join([dir_map["vcfdir"],
                         re.sub(r"\.bam", ".snps.vcf", bam.split('/')[-1])])
        indels = '/'.join([dir_map["vcfdir"],
                           re.sub(r"\.bam", ".indels.vcf",
                                  bam.split('/')[-1])])
        cmd1 = ' '.join([samtools, "mpileup -B -f", ref_genome, bam, "|",
                         java, "-jar", varscan,
                         "mpileup2snp --p-value 99e-02 --output-vcf 1",
                         ">", snps])
        cmd2 = ' '.join([samtools, "mpileup -B -f", ref_genome, bam, "|",
                         java, "-jar", varscan,
                         "mpileup2indel --p-value 99e-02 --output-vcf 1",
                         ">", indels])
        if dry:
            sys.stdout.write(cmd1 + '\n' + cmd2 + '\n')
        else:
            subprocess.call(cmd1, shell=True)
            subprocess.call(cmd2, shell=True)
        snp_vcfs.append(snps)
        indel_vcfs.append(indels)
    return snp_vcfs, indel_vcfs


def create_symlink(root, link_path, filename):
    '''
    Creates a symlink given root path, final path, and filename.
    '''
    os.symlink('/'.join([root, filename]),
               '/'.join([link_path, filename]))


def create_report(snps, indels, dir_map):
    '''
    Injects data into html template with jinja2.
    '''
    proj_dir = dir_map["projdir"]
    report_dir = dir_map["reportdir"]
    env = jinja2.Environment(loader=jinja2.FileSystemLoader(proj_dir))
    template = env.get_template("template.html")
    for variants in 
    pass


def get_rg_values(sample_name, fastq):
    '''
    Parses gzipped fastq for RG values and returns RG str for BWA.
    '''
    with gzip.open(fastq, 'rU') as handle:
        line = handle.readline()
        arow = line.strip('\n').split()
        info = arow[0].split(':')[1:]
        instrument_id = info[0]
        run_id = info[1]
        flowcell_id = info[2]
        flowcell_lane = info[3]
        index_seq = arow[1].split(':')[3]
        rgid = '.'.join([sample_name, flowcell_id, flowcell_lane])
    rglb = '.'.join([sample_name, run_id])
    rgpu = '.'.join([instrument_id,
                     flowcell_lane,
                     index_seq])
    rgsm = sample_name
    rgcn = "DFCI-CCCB"
    rgpl = "ILLUMINA"
    rg_vals = "@RG\\tID:" + rgid + "\\tPL:" + rgpl + "\\tLB:" + \
              rglb + "\\tSM:" + rgsm + "\\tCN:" + rgcn + "\\tPU:" + rgpu
    return rg_vals


def map_bams_to_sample(samplename_map, bams, bam_dir):
    '''
    Maps bam to sample_names.
    '''
    sample_to_bam_map = {}
    for key, value in samplename_map.items():
        bam = re.sub(r'_R[1-2]_[0-9]+\.fastq.gz', '', key).split('/')[-1] + ".bam"
        new_key = '/'.join([bam_dir, bam])
        sample_to_bam_map[new_key] = value
    bam_to_sample_map = defaultdict(list)
    for bam in bams:
        if bam in sample_to_bam_map:
            samplename = sample_to_bam_map[bam]
            bam_to_sample_map[samplename].append(bam)
    return bam_to_sample_map


def map_fastq_from_tsv(fastq_tsv):
    '''
    Parses tsv into map of sample name and lane specific fastq.
    '''
    fastq_map = defaultdict(list)
    with open(fastq_tsv, 'rU') as handle:
        for line in handle:
            arow = line.strip('\n').split('\t')
            sample = arow[0]
            fastq = arow[1]
            fastq_map[sample].append(fastq)
    return fastq_map


def map_fastq_pairs(fastqs):
    '''
    Gets fastqs from fastq dir and returns dict (k=first, v=second).
    '''
    first_read_fastqs = [f for f in fastqs if re.search("_R1_", f)]
    second_read_fastqs = [f for f in fastqs if re.search("_R2_", f)]
    basenames_map = {f.replace("_R2_", "") : f for f in second_read_fastqs}
    return map_read_by_basename(first_read_fastqs, basenames_map)


def map_read_by_basename(first_reads, basenames_map):
    '''
    Maps first reads to second reads and maps.
    '''
    mapped_fastq = []
    in_map = set([])
    for first in first_reads:
        basename = first.replace("_R1_", "")
        if basename in basenames_map:
            second = basenames_map[basename]
            mapped_fastq.append([first, second])
            in_map.add(basename)
        else:
            mapped_fastq.append([first, ""])
    seconds = [[basenames_map[b], ""] for b in basenames_map if b not in in_map]
    mapped_fastq.extend(seconds)
    return mapped_fastq


def merge_bams(bams_to_sample_map, dir_map, bin_map, dry=False):
    '''
    Merges the individual bams into a sample bam.
    '''
    merged_bams = []
    try:
        samtools = bin_map["samtools"]
    except KeyError as err:
        sys.stderr.write(err + '\n')
        sys.exit()
    for samplename, bams in bams_to_sample_map.items():
        input_bams = ' '.join(bams)
        output_bam = '/'.join([dir_map["bamdir"], samplename + ".bam"])
        cmd1 = ' '.join([samtools, "merge", output_bam, input_bams])
        cmd2 = ' '.join([samtools, "index", output_bam])
        if dry:
            sys.stdout.write(cmd1 + '\n' + cmd2 + '\n')
        else:
            subprocess.call(cmd1, shell=True)
            subprocess.call(cmd2, shell=True)
        merged_bams.append(output_bam)
    return merged_bams


def parse_vcf(vcf, is_indel=False, dry=False):
    '''
    Parses VCFs to a map of select features.
    '''
    header = ["chrom", "position", "ref allele", "alt allele",
              "alt freq", "p value", "total depth", "ref_depth",
              "alt depth"]
    variants = [header]
    with open(vcf, 'rU') as handle:
        for line in handle:
            if line[0] == "#":
                continue # skip header
            vrow = line.strip('\n').split('\t')
            chrom = vrow[0]
            pos = vrow[1]
            if is_indel:
                pos = str(int(vrow[1]) + 1)
            ref = vrow[3]
            alt = vrow[4]
            v_info = vrow[9].split(':')
            freq = v_info[6]
            pval = v_info[7]
            total_depth = v_info[3]
            ref_depth = v_info[4]
            allele_depth = v_info[5]
            variant = [chrom, pos, ref, alt, freq, pval, total_depth, 
                      ref_depth, allele_depth]
            #sys.stdout.write('\t'.join(variant) + '\n')
            variants.append(variant)
    return variants


def realign_indels(bams, ref_genome, bin_map, dry=False):
    '''
    Realigns indels with GATK.
    '''
    java = bin_map["java"]
    gatk = bin_map["gatk"]
    realn_bams = []
    for bam in bams:
        realn_intervals = re.sub(r"\.bam", "_indelrealn.intervals", bam)
        realn_bam = re.sub(r"\.bam", ".indelrealn.bam", bam)
        cmd1 = ' '.join([java, "-jar", gatk, "-T RealignerTargetCreator",
                         "-R", ref_genome, "-I", bam, "-o", realn_intervals])
        cmd2 = ' '.join([java, "-jar", gatk, "-T IndelRealigner",
                         "-R", ref_genome, "-I", bam,
                         "-targetIntervals", realn_intervals, "-o", realn_bam])
        if dry:
            sys.stdout.write(cmd1 + '\n' + cmd2 + '\n')
        else:
            subprocess.call(cmd1, shell=True)
            subprocess.call(cmd2, shell=True)
        realn_bams.append(realn_bam)
    return realn_bams


def setup_dir(cur_dir, out_dir_name):
    '''
    Sets up output directory in project directory.
    '''
    if cur_dir == '.':
        cur_dir = os.getcwd()
    if not os.path.isdir(cur_dir):
        sys.stderr.write("Error: project directory path does not exist.\n")
        sys.exit()
    try:
        os.makedirs('/'.join([cur_dir, out_dir_name]))
    except OSError as err:
        sys.stderr.write("%s\n" % err)
        sys.stderr.write("Error: output directory already exists.\n")
        sys.exit()
    out_dir = '/'.join([cur_dir, out_dir_name])
    bam_dir = '/'.join([out_dir, "bam_files"])
    indbam_dir = '/'.join([bam_dir, "individual_bam_files"])
    pileup_dir = '/'.join([out_dir, "pileups"])
    vcf_dir = '/'.join([out_dir, "vcfs"])
    report_dir = '/'.join([out_dir, "report_html"])
    for folder in [out_dir, bam_dir, indbam_dir, pileup_dir, vcf_dir,
                   report_dir]:
        os.makedirs(folder)
    return {"bamdir": bam_dir,
            "outdir": out_dir,
            "projdir": cur_dir,
            "indbamdir": indbam_dir,
            "pileupdir": pileup_dir,
            "vcfdir": vcf_dir,
            "reportdir": report_dir}


def main():
    '''
    Parses CLI args and central dispatch of functions.
    '''
    parser = argparse.ArgumentParser(description="Plexwell amplicon analysis")
    parser.add_argument("fastqtsv", metavar="FastqTSV",
                        help="lane specific fastq TSV")
    parser.add_argument("-d", "--projectdir",
                        metavar="PROJECT_DIR",
                        help="project directory; default: .")
    parser.add_argument("-o", "--outdir",
                        metavar="OUTPUT_DIR",
                        help="output dir; relative to -d; default: plexout")
    parser.add_argument("-g", "--genome",
                        help="genome [hg19,]; default: hg19")
    parser.add_argument("--dry",
                        action="store_true",
                        help="Creates dirs, but no file creation")
    parser.set_defaults(projectdir=".", outdir="plexout", genome="hg19")
    args = parser.parse_args()
    # Setup process
    dir_map = setup_dir(args.projectdir, args.outdir)
    #bwa = "/ifs/labs/cccb/projects/cccb/apps/bwa-0.7.12/bwa"
    #samtools = "/ifs/labs/cccb/projects/cccb/apps/samtools-1.2-5/samtools"
    #bwa = "~/bin/bwa"
    #samtools = "~/bin/samtools"
    bin_map = {"samtools": "~/bin/samtools",
               "bwa": "~/bin/bwa",
               "java": "~/bin/java",
               "gatk": "~/bin/gatk.jar",
               "varscan": "~/bin/varscan.jar",}
    # TODO change ref at later point
    #ref_map = {"hg19": ""}
    indvbams_to_samplename_map = align_reads(dir_map, args.fastqtsv,
                                             "hg19.fa", bin_map,
                                             dry=args.dry)
    merged_bams = merge_bams(indvbams_to_samplename_map, dir_map, bin_map,
                             dry=args.dry)
    indelrealn_bams = realign_indels(merged_bams, "hg19.fa", bin_map,
                                     dry=args.dry)
    snps, indels = call_variants_from_bams(indelrealn_bams, "hg19.fa", dir_map,
                                           bin_map, dry=args.dry)
    for snp in snps:
        parse_vcf(snp)


if __name__ == "__main__":
    main()
