#!/bin/env python

"""
Created on May 3, 2016

Module that helps containerize the CLIP density information.

@author: Gabe, Brian
"""

import numpy as np
import pyBigWig
import pysam
import os
import math

class Density:
    def __init__(self, bam, name=None):
        if bam is not None:
            self.bam = pysam.AlignmentFile(bam)
        else:
            print("warning no bam file!")
            
        self.name = name if name is not None else os.path.basename(bam)
        
    def values(self, chrom, start, end, strand):
        return 0
    
    def pseudocount(self):
        """
        Returns the minimum normalized pseudocount of 1 read.

        Returns
        -------
        rpm : float
        """
        return 1000000.0 / self.bam.mapped

    def total_mapped(self):
        """
        Returns the number of mapped reads

        Returns
        -------
        mapped :
        """
        return self.bam.mapped

    def rpm_to_r(self, rpm):
        """
        Returns the raw read representation given a pseudocount

        Parameters
        ----------
        rpm : float
            rpm-normalized read density
        Returns
        -------
        r : float
        """
        return (rpm * 1000000.0) / self.bam.count()
    

class ReadDensity(Density):
    """
    ReadDensity class
    Attributes:
        self.pos(positive *.bw file)
        self.neg(negative *.bw file)
    """

    def __init__(self, pos, neg, name=None, bam=None):
        try:
            super().__init__(bam, name)
            self.pos = pyBigWig.open(pos)
            self.neg = pyBigWig.open(neg)

        except Exception as e:
            print("couldn't open the bigwig files!")
            print(e)


    def values(self, chrom, start, end, strand, zeroes=False):
        """

        Parameters
        ----------
        chrom : basestring
            (eg. chr1)
        start : int
            0-based start (first position in chromosome is 0)
        end : int
            1-based end (last position is not included)
        strand : str
            either '+' or '-'
        zeroes : boolean
            if True, return nan values (zero signal) as 0
        Returns
        -------
        densites : list
            values corresponding to density over specified positions.
        """

        try:
            if strand == "+":
                values = self.pos.values(chrom, start, end)
            elif strand == "-":
                values = list(reversed(self.neg.values(chrom, start, end)))
            else:
                print("Strand neither + or -")
                return 1
            
            if zeroes:
                return [0 if math.isnan(x) else x for x in values]
            else:
                return values
        except RuntimeError:
            # usually occurs when no chromosome exists in the bigwig file
            return [np.NaN] * abs(start - end)


class ReadDensityUnstranded(Density):

    def __init__(self, bigwig, name=None):
        try:
            super().__init__(bam, name)
            self.bigwig = pyBigWig.open(bigwig)
        except Exception as e:
            print("couldn't open the bigwig files!")
            print(e)

            
    def values(self, chrom, start, end):
        """

        Parameters
        ----------
        chrom : basestring
            (eg. chr1)
        start : int
            0-based start (first position in chromosome is 0)
        end : int
            1-based end (last position is not included)

        Returns
        -------
        densites : list
            values corresponding to density over specified positions. Assume all values positive.
        """

        try:
            return self.bigwig.values(chrom, start, end)
        except RuntimeError:
            # usually occurs when no chromosome exists in the bigwig file
            return [np.NaN] * abs(start - end)
