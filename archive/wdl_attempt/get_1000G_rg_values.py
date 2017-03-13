#! /usr/bin/python
'''
Parses gzipped FASTQ file for RG values from read headers.
'''

import gzip
import sys


def print_rg(filename, sample_name):
    '''
    Parses first line plus sample_name to print RG values.
    '''
    #@ERR031964.6 B81C73ABXX:2:1101:1270:2549#TTAGGCAN/1
    with gzip.open(filename, "rU") as handle:
        line = handle.readline()
        arow = line.strip('\n').split()
        info = arow[1][:-2].split(':')
        instrument_id = info[0]
        #run_id = info[1]
        #flowcell_id = info[2]
        flowcell_lane = info[1]
        index_seq = arow[1][:-2].split('#')[1]
    rgid = '.'.join([sample_name, flowcell_lane])
    rglb = '.'.join([sample_name, instrument_id, flowcell_lane])
    rgpu = '.'.join([instrument_id,
                     flowcell_lane,
                     index_seq])
    rgsm = sample_name
    rgcn = "DFCI-CCCB"
    rgpl = "ILLUMINA"
    #out = "@RG\\tID:" + rgid + "\\tPL:" + rgpl + "\\tLB:" + \
    #      rglb + "\\tSM:" + rgsm + "\\tCN:" + rgcn + "\\tPU:" + \
    #      rgpu
    #sys.stdout.write(out)
    out = [rgid, rgpl, rglb, rgsm, rgcn, rgpu]
    sys.stdout.write('\n'.join(out) + '\n')


def main(sa):
    '''
    Parses CLI input
    '''
    sample_name = sa[0]
    filename = sa[1]
    print_rg(filename, sample_name)


if __name__ == "__main__":
    main(sys.argv[1:])
