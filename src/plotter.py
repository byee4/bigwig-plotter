#!/usr/bin/env python
# encoding: utf-8
'''
@author:     brian
@copyright:  2023 Yeo lab. All rights reserved.
@license:    license
@contact:    bay001@health.ucsd.edu
@deffield    updated: 2023-04-07
'''
import matplotlib
matplotlib.use('Agg')

from matplotlib import rc
rc('text', usetex=False)
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import logging
import os
import argparse
from collections import OrderedDict
from . import ReadDensity
from .make_bigwig_files import genome_coverage_bed, sort_bedgraph, bed_graph_to_big_wig, get_norm_constant, check_for_index

from subprocess import call

logger = logging.getLogger('plot_features')

def plot_features(features, bam_file, output_file, width=7, height=2):
    """
    params: 
        features: dict
    """
    coords = features.keys()
    w = width
    h = height*len(coords)  # scale with number of features to plot
    
    fig = plt.figure(figsize=(w, h))
    fig.suptitle(f"Overlapping {bam_file} with {len(coords)} features")

    # Set up the Grid (one row per feature)
    axs = []
    full_grid = gridspec.GridSpec(
        len(coords), 1, height_ratios=[1] * len(coords),
    )
    full_grid.update(wspace=0.025, hspace=1) # set the spacing between axes. 

    
    i = 0  # counter for number of rows/features we want to plot
    for key, values in features.items():
        axs.append(fig.add_subplot(full_grid[i]))
        axs[i].plot(features[key])
        axs[i].set_title(key)
        i += 1

    fig.savefig(output_file)
    
    
def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--bam",
        required=False
    )
    parser.add_argument(
        "--output",
        required=False,
        default='image.svg'
    )
    parser.add_argument(
        "--five_prime",
        required=False,
        default=False,
        action='store_true',
        help='Compute coverage using the 5 prime end [1nt] of each read (default: False, use entire read)'
    )
    parser.add_argument(
        "--chrom_sizes",
        required=True,
    )
    parser.add_argument(
        "--rpm",
        required=False,
        action='store_true',
        default=False,
        help='RPM-normalize values (default: False, leave unscaled)'
    )
    parser.add_argument(
        "--regions-bed", "--regions_bed", "-r",
        required=True,
        help='BED6 file containing regions to overlap BAM signal with. Each region will be its own plot.'
    )
    parser.add_argument(
        "--width",
        required=False,
        default=12,
        type=float,
        help='Width of the total figure (default: 12)'
    )
    parser.add_argument(
        "--height",
        required=False,
        default=2,
        type=float,
        help='Height PER REGION (default: 2). So if there are 5 regions to plot, the total height will be 10.'
    )
    
    # Process arguments
    args = parser.parse_args()
    output_file = args.output
    bam_file = args.bam
    five_prime = args.five_prime
    chrom_sizes = args.chrom_sizes
    bed_file = args.regions_bed
    rpm = args.rpm
    width = args.width
    height = args.height
    
    check_for_index(bam_file)
    
    # Scale values if rpm else leave unscaled
    if rpm:
        norm_constant = get_norm_constant(bam_file)
        neg_norm_constant = norm_constant * -1
        bedGraphFilePos = os.path.splitext(bam_file)[0] + ".norm.pos.bg"
        bedGraphFileNeg = os.path.splitext(bam_file)[0] + ".norm.neg.bg"

        bedGraphFilePosSorted = os.path.splitext(bam_file)[0] + ".sorted.norm.pos.bg"
        bedGraphFileNegSorted = os.path.splitext(bam_file)[0] + ".sorted.norm.neg.bg"

        bw_pos = os.path.splitext(bam_file)[0] + ".norm.pos.bw"
        bw_neg = os.path.splitext(bam_file)[0] + ".norm.neg.bw"
    
    else:
        norm_constant = 1
        neg_norm_constant = -1
        bedGraphFilePos = os.path.splitext(bam_file)[0] + ".pos.bg"
        bedGraphFileNeg = os.path.splitext(bam_file)[0] + ".neg.bg"

        bedGraphFilePosSorted = os.path.splitext(bam_file)[0] + ".sorted.pos.bg"
        bedGraphFileNegSorted = os.path.splitext(bam_file)[0] + ".sorted.neg.bg"

        bw_pos = os.path.splitext(bam_file)[0] + ".pos.bw"
        bw_neg = os.path.splitext(bam_file)[0] + ".neg.bw"
        
    # Generate bigwig files

    if not os.path.exists(bw_pos) or not os.path.exists(bw_neg):
        
        genome_coverage_bed(
            in_bam=bam_file,
            out_bed_graph=bedGraphFilePos,
            strand="+", 
            scale=norm_constant,
            five_prime=args.five_prime
        )
        sort_bedgraph(bedGraphFilePos, bedGraphFilePosSorted)
        bed_graph_to_big_wig(bedGraphFilePosSorted, chrom_sizes, bw_pos)

        genome_coverage_bed(
            in_bam=bam_file,
            out_bed_graph=bedGraphFileNeg,
            strand="-", 
            scale=neg_norm_constant,
            five_prime=five_prime
        )
        sort_bedgraph(bedGraphFileNeg, bedGraphFileNegSorted)
        bed_graph_to_big_wig(bedGraphFileNegSorted, chrom_sizes, bw_neg)
    
    rdd = ReadDensity.ReadDensity(pos=bw_pos, neg=bw_neg, bam=bam_file)
    
    features = {}
    with open(bed_file, 'r') as f:
        for line in f:
            chrom, start, end, name, score, strand = line.rstrip().split('\t')
            start = int(start)
            end = int(end)
            features[f"{chrom}:{start}-{end}:{strand}"] = rdd.values(chrom, start, end, strand, zeroes=True)

    plot_features(features, os.path.basename(bam_file), output_file, width, height)
if __name__ == "__main__":
    main()
