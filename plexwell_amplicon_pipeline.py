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
import shutil
import subprocess
import sys


def align_reads(fastq_tsv, ref_genome, dir_map, bin_map, suffix_map, dry=False):
    '''
    Aligns reads from fastqs.
    '''
    bwa = bin_map["bwa"]
    samtools = bin_map["samtools"]
    fastq_map = map_fastq_from_tsv(fastq_tsv)
    inverted_fastq_map = dict((v, k) for k in fastq_map for v in fastq_map[k])
    fastqs = itertools.chain.from_iterable(fastq_map.values())
    fastq_pairs = map_fastq_pairs(fastqs)
    bams = defaultdict(list)
    for first, second in fastq_pairs.items():
        sample_name = inverted_fastq_map[first]
        bam_basename = re.sub(r'_R[1-2]_[0-9]+\.fastq.gz',
                              '', first).split('/')[-1]
        bam = dir_map["indbamdir"] + '/' + bam_basename + suffix_map["bam"]
        rg_vals = get_rg_values(sample_name, first)
        cmd = ' '.join([bwa, 'mem -R \"%s\"' % rg_vals,
                        ref_genome, first, second, '|',
                        samtools, "view -bht", ref_genome, "|",
                        samtools, "sort", ">", bam])
        if dry:
            sys.stdout.write(cmd + '\n')
        else:
            subprocess.call(cmd, shell=True)
        bams[sample_name].append(bam)
    return bams


def call_variants(pileups_map, ref_genome, dir_map, bin_map, suffix_map, dry=False):
    '''
    Pipes samtools mpileup of bam to varscan.
    '''
    java = bin_map["java"]
    varscan = bin_map["varscan"]
    snp_map = {}
    indel_map = {}
    for samplename, pileup in pileups_map.items():
        snps = '/'.join([dir_map["vcfdir"],
                         samplename + suffix_map["snps"]])
        indels = '/'.join([dir_map["vcfdir"],
                           samplename + suffix_map["indels"]])
        cmd1 = ' '.join([java, "-jar", varscan, "mpileup2snp", pileup,
                         "--p-value 99e-02 --output-vcf 1", ">", snps])
        cmd2 = ' '.join([java, "-jar", varscan, "mpileup2indel", pileup,
                         "--p-value 99e-02 --output-vcf 1", ">", indels])
        if dry:
            sys.stdout.write(cmd1 + '\n' + cmd2 + '\n')
        else:
            subprocess.call(cmd1, shell=True)
            subprocess.call(cmd2, shell=True)
        snp_map[samplename] = snps
        indel_map[samplename] = indels
    return snp_map, indel_map


def cluster_filter(variant_map, cluster_map):
    '''
    Filters variants by intersection with clusters.
    '''
    filtered_variant_map = {}
    for samplename, variants in variant_map.items():
        clusters = cluster_map[samplename]
        filtered_variants = [v for i, v in enumerate(variants)
                             if i == 0 or has_intersect(v, clusters)]
        filtered_variant_map[samplename] = filtered_variants
        print "DEBUG VARIANTS", filtered_variant_map
        print "DEBUG CLUSTERS", clusters
    return filtered_variant_map


def cluster_regions(bams_map, min_mean_depth, dir_map, bin_map, suffix_map,
                    dry=False):
    '''
    Determines specific amplicon by min_mean_coverage.
    '''
    bedtools = bin_map["bedtools"]
    bed_map = {}
    cluster_map = {}
    for samplename, bam in bams_map.items():
        bed = '/'.join([dir_map["coveragedir"],
                        samplename + suffix_map["coveragebed"]])
        cmd = ' '.join([bedtools, "merge -i", bam, "| head -n -1 |",
                        bedtools, "coverage -a - -b", bam, ">", bed])
        if dry:
            sys.stdout.write(cmd + '\n')
        else:
            subprocess.call(cmd, shell=True)
        bed_map[samplename] = bed
        cluster_map[samplename] = filter_clusters(bed, min_mean_depth)
    return bed_map, cluster_map


def create_report(snps, indels, dir_map, dry=False):
    '''
    Injects data into html template with jinja2.
    '''
    this_dir = os.path.dirname(os.path.realpath(__file__))
    lib_dir = os.path.join(this_dir, 'lib')
    report_dir = dir_map["reportdir"]
    lib_destination = os.path.join(report_dir, 'lib')
    report = '/'.join([report_dir, "html_report.html"])
    env = jinja2.Environment(loader=jinja2.FileSystemLoader(this_dir))
    # can't seem to find template.html
    template = env.get_template("template.html")
    print snps
    samples = {samplename : {"header" : snps[samplename][0],
                             "snps" : snps[samplename][1:],
                             "indels" : indels[samplename][1:],
                             "samplename": samplename.split('/')[-1]}
               for samplename in snps.keys()}
    context = {"samples": samples}
    if not dry:
        with open(report, 'w') as outfile:
            outfile.write(template.render(context))
        shutil.copytree(lib_dir, lib_destination)


def filter_clusters(bed, min_mean_depth):
    '''
    Parses BED file to filter for high coverage regions.
    '''
    #clusters = defaultdict(list)
    clusters = {} # for some reason defaultdict is buggy
    with open(bed, 'rU') as bedfile:
        for line in bedfile:
            arow = line.strip('\n').split('\t')
            chrom = arow[0]
            try:
                begin = int(arow[1])
                end = int(arow[2])
                total_depth = int(arow[3])
                bases_covered = int(arow[5])
            except ValueError as err:
                sys.stderr.write(str(err) + "\nError in coverage bed file.\n")
                sys.exit()
            if total_depth/float(bases_covered) >= min_mean_depth:
                #clusters[chrom].append((begin, end))
                if chrom in clusters:
                    clusters[chrom].append((begin, end))
                else:
                    clusters[chrom] = [(begin, end)]
    return clusters


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


def has_intersect(variant, clusters):
    '''
    Checks for any intersects between variant and any cluster.
    '''
    intersected = False
    vk = variant["chrom"]
    try:
        vpos = int(variant["pos"])
    except ValueError as err:
        sys.stderr.write(str(err) + "\n")
        sys.exit()
    return_bool = False
    if vk in clusters:
        for cluster in clusters[vk]:
            begin = cluster[0]
            end = cluster[1]
            if begin <= vpos <= end:
                return_bool = True
                break
    print "DEBUG BOOL", return_bool, vk, vpos, clusters
    return return_bool
    # elegance not working...
    #return vk in clusters or (len(clusters[vk]) > 0 
    # and any(b[0] <= vpos <= b[1] for b in clusters[vk]))


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
    basename_second_map = {f.replace("_R2_", "") : f 
                           for f in second_read_fastqs}
    mapped_pairs = {f : basename_second_map[f.replace("_R1_", "")]
                    for f in first_read_fastqs
                    if f.replace("_R1_", "") in basename_second_map}
    first_only = {f : "" for f in first_read_fastqs
                  if f.replace("_R1_", "") not in basename_second_map}
    second_only = {f : "" for f in second_read_fastqs
                   if f.replace("_R2_", "" not in mapped_pairs)}
    mapped_pairs.update(first_only)
    mapped_pairs.update(second_only)
    return mapped_pairs


def merge_bams(indv_bams_map, dir_map, bin_map, suffix_map, dry=False):
    '''
    Merges the individual bams into a sample bam.
    '''
    merged_bams = {}
    try:
        samtools = bin_map["samtools"]
    except KeyError as err:
        sys.stderr.write(err + '\n')
        sys.exit()
    for samplename, bams in indv_bams_map.items():
        input_bams = ' '.join(bams)
        output_bam = '/'.join([dir_map["bamdir"], 
                               samplename + suffix_map["bam"]])
        cmd1 = ' '.join([samtools, "merge", output_bam, input_bams])
        cmd2 = ' '.join([samtools, "index", output_bam])
        if dry:
            sys.stdout.write(cmd1 + '\n' + cmd2 + '\n')
        else:
            subprocess.call(cmd1, shell=True)
            subprocess.call(cmd2, shell=True)
        merged_bams[samplename] = output_bam
    return merged_bams


def parse_vcf(vcf, is_indel=False, dry=False):
    '''
    Parses VCFs to a map of select features.
    '''
    header = {"chrom": "chrom",
              "pos": "position",
              "ref": "ref allele",
              "alt": "alt allele",
              "freq": "alt freq",
              "pval": "p value",
              "total_depth": "total depth",
              "ref_depth": "ref depth",
              "allele_depth": "alt depth"}
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
            variant = {"chrom": chrom,
                       "pos": pos,
                       "ref": ref,
                       "alt": alt,
                       "freq": freq,
                       "pval": pval,
                       "total_depth": total_depth,
                       "ref_depth": ref_depth,
                       "allele_depth": allele_depth}
            variants.append(variant)
    return variants


def pileup(bam_map, ref_genome, dir_map, bin_map, suffix_map, dry=False):
    '''
    Creates pileups from bam files with samtools mpileup.
    '''
    samtools = bin_map["samtools"]
    pileup_map = {}
    for samplename, bam in bam_map.items():
        pileup = '/'.join([dir_map["pileupdir"],
                           samplename + suffix_map["pileup"]])
        cmd = ' '.join([samtools, "mpileup -B -f", ref_genome, bam,
                        ">", pileup])
        if dry:
            sys.stdout.write(cmd + "\n")
        else:
            subprocess.call(cmd, shell=True)
        pileup_map[samplename] = pileup
    return pileup_map


def realign_indels(bam_map, ref_genome, dir_map, bin_map, suffix_map,
                   dry=False):
    '''
    Realigns indels with GATK.
    '''
    java = bin_map["java"]
    gatk = bin_map["gatk"]
    realn_bams = {}
    for samplename, bam in bam_map.items():
        realn_intervals = dir_map["bamdir"] + '/' + samplename + \
                          suffix_map["indelrealnintervals"]
        realn_bam = dir_map["bamdir"] + '/' + samplename + \
                    suffix_map["indelrealn"]
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
        realn_bams[samplename] = realn_bam
    return realn_bams


def setup_dir(cur_dir, out_dir_name, dry=False):
    '''
    Sets up output directory in project directory.
    '''
    if cur_dir == '.':
        cur_dir = os.getcwd()
    if not os.path.isdir(cur_dir):
        sys.stderr.write("Error: project directory path does not exist.\n")
        sys.exit()
    out_dir = '/'.join([cur_dir, out_dir_name])
    bam_dir = '/'.join([out_dir, "bam_files"])
    indbam_dir = '/'.join([bam_dir, "individual_bam_files"])
    pileup_dir = '/'.join([out_dir, "pileups"])
    vcf_dir = '/'.join([out_dir, "vcfs"])
    report_dir = '/'.join([out_dir, "report_html"])
    coverage_dir = '/'.join([out_dir, "coverage"])
    if not dry:
        try:
            for folder in [out_dir, bam_dir, indbam_dir, pileup_dir,
                           coverage_dir, vcf_dir, report_dir]:
                os.makedirs(folder)
        except OSError as err:
            sys.stderr.write("%s\n" % err)
            sys.stderr.write("Error: %s directory already exists.\n" % folder)
            sys.exit()
    return {"bamdir": bam_dir,
            "outdir": out_dir,
            "projdir": cur_dir,
            "indbamdir": indbam_dir,
            "pileupdir": pileup_dir,
            "vcfdir": vcf_dir,
            "coveragedir": coverage_dir,
            "reportdir": report_dir}


def main():
    '''
    Parses CLI args and central dispatch of functions.
    '''
    # Arg parsing
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
    parser.add_argument("-m", "--mindepth",
                        type=int,
                        metavar="MIN_MEAN_DEPTH",
                        help="min mean coverage for clusters; default 200")
    parser.add_argument("--dry",
                        action="store_true",
                        help="Creates dirs, but no file creation")
    parser.set_defaults(projectdir=".", outdir="plexout", genome="hg19",
                        mindepth=200)
    args = parser.parse_args()
    # Set up globally used maps.
    bin_map = {"samtools": "~/bin/samtools",
            "bwa": "~/bin/bwa",
            "java": "~/bin/java",
            "gatk": "~/bin/gatk.jar",
            "varscan": "~/bin/varscan.jar",
            "bedtools": "~/bin/bedtools"}
    suffix_map = {"indelrealn": ".indelrealn.bam",
                  "indelrealnintervals": ".indelrealn.intervals",
                  "coveragebed": ".bed",
                  "bam": ".bam",
                  "snps": ".snps.vcf",
                  "indels": ".indels.vcf",
                  "pileup": ".pileup"}
    dir_map = setup_dir(args.projectdir, args.outdir, dry=args.dry)
    # Processing
    indv_bams_map = align_reads(args.fastqtsv, "hg19.fa", dir_map, bin_map,
                                suffix_map, dry=args.dry)
    merged_bams_map = merge_bams(indv_bams_map, dir_map, bin_map,
                                 suffix_map, dry=args.dry)
    realn_bams_map = realign_indels(merged_bams_map, "hg19.fa", dir_map,
                                    bin_map, suffix_map, dry=args.dry)
    _, cluster_map = cluster_regions(realn_bams_map, int(args.mindepth),
                                     dir_map, bin_map, suffix_map, dry=args.dry)
    pileups_map = pileup(realn_bams_map, "hg19.fa", dir_map, bin_map,
                         suffix_map, dry=args.dry)
    snp_vcf_map, indel_vcf_map = call_variants(pileups_map, "hg19.fa",
                                               dir_map, bin_map, suffix_map,
                                               dry=args.dry)
    snp_map = {samplename : parse_vcf(vcf, dry=args.dry)
               for samplename, vcf in snp_vcf_map.items()}
    indel_map = {samplename : parse_vcf(vcf, is_indel=True, dry=args.dry)
                 for samplename, vcf in indel_vcf_map.items()}
    f_snp_map, f_indel_map = [cluster_filter(m, cluster_map)
                              for m in [snp_map, indel_map]]
    create_report(f_snp_map, f_indel_map, dir_map)


if __name__ == "__main__":
    main()
