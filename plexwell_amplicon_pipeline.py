#!/usr/bin/python
'''
Main entry for plexwell amplicon analysis.
'''

from collections import defaultdict
import argparse
import os
import re
import subprocess
import sys


def align_reads(project_dir, bam_dir, ref_genome):
    '''
    Aligns reads from fastqs.
    '''
    fastq_dir = '/'.join([project_dir, "fastq_symlinks"])
    # Checks for existence of fastq_dir, otherwise creates it in output_dir.
    if not os.path.isdir(fastq_dir):
        sys.stderr.write("Error: fastq symlink dir does not exist.\n")
        sys.stderr.write("       Creating symlink dir.\n")
        # TODO create symlinks for fastq
    mapped_fastqs = get_fastqs(fastq_dir)
    for first, second in mapped_fastqs:
        basename = re.sub('_R[1-2]_.final.fastq.gz', '', first)
        bwa = "/ifs/labs/cccb/projects/cccb/apps/bwa-0.7.12/bwa"
        samtools = "/ifs/labs/cccb/projects/cccb/apps/samtools-1.2-5/samtools"
        sam = '.'.join([basename, "sam"])
        first, second, sam = ['/'.join([bam_dir, b])
                              for b in [first, second, sam]]
        cmd = ' '.join([bwa, 'mem -t 10', ref_genome, first, second, '|',
                        samtools, "view -bht", "|",
                        samtools, "sort -@ 10", ">", sam])
        print(cmd)


def create_symlink(root, link_path, file):
    '''
    Creates a symlink given root path, final path, and filename.
    '''
    os.symlink('/'.join([root, file]),
               '/'.join([link_path, file]))


def get_fastqs(fastq_dir):
    '''
    Gets fastqs from fastq dir and returns dict (k=first, v=second).
    '''
    fastqs = [f for f in os.listdir(fastq_dir)
              if os.path.isfile(os.path.join(fastq_dir, f))]
    first_read_fastqs = [f for f in fastqs if re.search("_R1_", f)]
    second_read_fastqs = [f for f in fastqs if re.search("_R2_", f)]
    basenames_map = {f.replace("_R2_", "") : f for f in second_read_fastqs}
    return map_fastq_pairs(first_read_fastqs, basenames_map)


def map_fastq_pairs(first_reads, basenames_map):
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
    os.makedirs(bam_dir)
    return {"bamdir": bam_dir}


def main():
    '''
    Parses CLI args and central dispatch of functions.
    '''
    parser = argparse.ArgumentParser(description="Plexwell amplicon analysis")
    parser.add_argument("-d", "--projectdir",
                        metavar="PROJECT_DIR",
                        help="project directory; default: .")
    parser.add_argument("-o", "--outdir",
                        metavar="OUTPUT_DIR",
                        help="output dir; relative to -d; default: plexout")
    parser.add_argument("-g", "--genome",
                        help="genome [hg19,]; default: hg19")
    parser.set_defaults(projectdir=".", outdir="plexout", genome="hg19")
    args = parser.parse_args()

    # Setup process
    setup_dir(args.projectdir, args.outdir)
    align_reads(args.projectdir, args.outdir, "foo.fa")


if __name__ == "__main__":
    main()
