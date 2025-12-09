#!coding:utf-8
import pysam
import re
import argparse
#----------------------------------------------
#Project: 商检-基础流程-基础分析
#Description: 酶切过滤程序的核心模块-pysam, 不可单独使用，对pysam进行二次开发，方便主程序使用
#Usage: 不可单独使用
#Author: 骆磊
#------------------------------------------------

#--------------pysam method
#alignment_file=pysam.AlignmentFile(bam)#返回迭代器，每一次迭代是一条read
#allreads=bamfile_handle.fetch('chr1',start=2036839,end=2036840) 以0为坐标起始, 即使是个D，也算mapping上了，这个应该算的是原本能map上的区域
#sequence=read.query_sequence 返回的是bam文件中的read序列，与ref序列同向，包含softclip但不含hardclip
#align_pair=read.get_aligned_pairs(with_seq) 不一致的位点是小写，但有少数情况不一致也是大写，莫名其妙
#
#
#
class bamhunter(object):
    '''
    v0: 20220729
    这个class是为了捕捉错配的reads而写的
    应该包含的功能有:
    1. cover到某个/些点的reads列表
    2. 输入readid和位置就能返回这个位置的比对情况
    '''
    def __init__(self,bamfile):
        self.bamfile=bamfile
        self.pysam_bam=self.__readbam(self.bamfile)


    def __readbam(self,bam):
        '''
        pysam的方式读入bam，节省时间和内存
        '''
        read_info_dict={}
        read_map_dict={}
        bamfile_handle=pysam.AlignmentFile(bam)
        return bamfile_handle
    def get_pos_reads(self,chrname,start,end):
        '''
        指定一个位点，获取cover这个位点的所有reads，不管这个碱基有没有比对上
        '''
#        allreads=bamfile_handle.fetch('chr1',start=2036839,end=2036840)
#        在外面调用时pos是以1为坐标起始的
        start-=1
        end-=1
        allreads=self.pysam_bam.fetch(chrname,start=start,end=end)
        return allreads#返回覆盖到这个位点的所有reads
                    
    def get_read_pos(self,refpos,read_pysam):
        '''
        指定一条reads，获取在这个位点，read的第几个碱基，ref的第几个碱基，ref碱基是什么
        '''
        align_status=read_pysam.get_aligned_pairs(with_seq=True)
        refpos-=1#外面调用时是1-based
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
        指定一条reads，获取在这个位点，read的第几个碱基，ref的第几个碱基，ref碱基是什么
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
        指定一条read的第几个碱基，输出这个碱基对应的参考基因组的碱基
        '''
        query_seq=read_pysam.query_sequence
        refbase=query_seq[readpos-1]
        return refbase
    def get_read_seq(self,read_pysam,start,end):
        '''
        取出指定区间的read序列
        '''
        start-=1#换成0-based
        query_seq=read_pysam.query_sequence
        if start==end:
            return '-'
        if end=="end":
            return query_seq[start:]
        else:
            return query_seq[start:end]
