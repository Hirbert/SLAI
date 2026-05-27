#!coding:utf-8
import pysam
import re
import argparse
#----------------------------------------------
#Project: Commercial Inspection - Basic Process - Basic Analysis
#Description: Core module of enzyme digestion filtering program - pysam, cannot be used alone, secondary development of pysam to facilitate use by main program
#Usage: Cannot be used alone
#Author: Luo Lei
#------------------------------------------------

#--------------pysam method
#alignment_file=pysam.AlignmentFile(bam)#Returns an iterator, each iteration is one read
#allreads=bamfile_handle.fetch('chr1',start=2036839,end=2036840) 0-based coordinate, even if it's a D, it is considered mapped, this should calculate the originally mappable region
#sequence=read.query_sequence Returns the read sequence in the bam file, same orientation as reference sequence, includes softclip but not hardclip
#align_pair=read.get_aligned_pairs(with_seq) Mismatched sites are lowercase, but in rare cases mismatches are also uppercase, inexplicable
#
#
#
class bamhunter(object):
    '''
    v0: 20220729
    This class is written to capture mismatched reads
    Functions that should be included:
    1. List of reads covering certain point(s)
    2. Input readid and position, return alignment status at that position
    '''
    def __init__(self,bamfile):
        self.bamfile=bamfile
        self.pysam_bam=self.__readbam(self.bamfile)


    def __readbam(self,bam):
        '''
        Read bam using pysam method, saves time and memory
        '''
        read_info_dict={}
        read_map_dict={}
        bamfile_handle=pysam.AlignmentFile(bam)
        return bamfile_handle
    def get_pos_reads(self,chrname,start,end):
        '''
        Specify a position, get all reads covering this position, regardless of whether the base is aligned
        '''
#        allreads=bamfile_handle.fetch('chr1',start=2036839,end=2036840)
#        When called externally, pos is 1-based coordinate
        start-=1
        end-=1
        allreads=self.pysam_bam.fetch(chrname,start=start,end=end)
        return allreads#Returns all reads covering this position
                    
    def get_read_pos(self,refpos,read_pysam):
        '''
        Specify a read, get at this position: which base of the read, which base of the reference, and what the reference base is
        '''
        align_status=read_pysam.get_aligned_pairs(with_seq=True)
        refpos-=1#External call is 1-based
        query_seq=read_pysam.query_sequence
        for pileup in align_status:
            if refpos==pileup[1]:
                readpos=pileup[0]
                try:
                    readbase=query_seq[readpos]
                except(TypeError):
                    return False
                return [readpos+1,readbase,refpos+1,pileup[2].upper()]
        return False
        
    def new_get_read_pos(self,refpos,read_pysam,align_status):
        '''
        Specify a read, get at this position: which base of the read, which base of the reference, and what the reference base is
        '''
        query_seq=read_pysam.query_sequence
        for pileup in align_status:
            if refpos==pileup[2] and pileup[3]!="-":
                readpos=pileup[0]
                try:
                    readbase=query_seq[readpos]
                except(TypeError):
                    return False
                return [readpos+1,readbase,refpos+1,pileup[3].upper()]
    def get_read_base(self,readpos,read_pysam):
        '''
        Specify a base position of a read, output the reference genome base corresponding to this base
        '''
        query_seq=read_pysam.query_sequence
        refbase=query_seq[readpos-1]
        return refbase
    def get_read_seq(self,read_pysam,start,end):
        '''
        Extract the read sequence within the specified interval
        '''
        start-=1#Convert to 0-based
        query_seq=read_pysam.query_sequence
        if start==end:
            return '-'
        if end=="end":
            return query_seq[start:]
        else:
            return query_seq[start:end]
