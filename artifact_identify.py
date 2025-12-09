#!coding:utf-8
#----------------------------------------------
#Project: 商检-基础流程-基础分析
#Description: 酶切过滤程序的核心模块-比对部分, 也可单独使用，输出snv附近的负链比对结果
#Usage: ./artifact_identify.py -h
#Author: 骆磊
#
#
#
#------------------------------------------------
import sys
import os
import pysam
import re
import SmithWaterman as SmithWaterman
import argparse
import support_reads as SR



def argument_parser():
    '''
    外部传参
    '''
    parser = argparse.ArgumentParser(description="check whether snv is potentially from artifact reads")
    parser.add_argument('--snv',help='file or snv(chr1:1000:A:T mean A->T occured at chr1 1000)',required=True)
    #/beegfs/work/commercial_test/cupcake/databases/gatk_bundle/2.8/hg19/ucsc.hg19.noconfig.fasta
    parser.add_argument('--fasta',help='fasta file used',required=True)
    parser.add_argument('--outdir',help='outdir, default current directory',default='./')
    parser.add_argument('--outfile',help='prefix of outfie',required=True)
    parser.add_argument('--flank',type=int,help='search region',default=300)
    parser.add_argument('--present',action='store_true',help='show artifact')
    argv=vars(parser.parse_args())
    return argv


def complement_reverse(seq):
    '''
    输入序列，返回对应的反向互补序列
    '''
    tmp=seq.upper()
    tmp=tmp.replace('A', 't')
    tmp=tmp.replace('T', 'a')
    tmp=tmp.replace('C', 'g')
    tmp=tmp.replace('G', 'c')
    tmp=tmp[::-1]
    return tmp.upper()

def read_fasta(fasta):
    '''
    返回字典{chr1:'AGCATGCTAGCTACGATC...',....}
    '''
    ref_dict={}#{chr:sequence,...}
    fasta_handle=pysam.FastxFile(fasta)
    for each_chr in fasta_handle:
        ref_dict[each_chr.name]=each_chr.sequence.upper()
    return ref_dict

def get_fasta_sequence(ref_dict,chrname,begin,end):
    '''
    输入基因组上的坐标，根据坐标返回序列
    '''
    return ref_dict[chrname][begin-1:end]

def structure_identify(site,ref_dict):
    '''
    提取snv附近的参考基因组序列，记为A，获取A的反向互补序列，记为B，将A和B用动态规划的方法进行比对
    '''
    sites=site.split(':')
    chrname=sites[0]
    pos=int(sites[1])
    ref_base=sites[2]
    alt_base=sites[3]
    refseq=get_fasta_sequence(ref_dict,chrname,pos-flank,pos+flank-1+len(ref_base))#变异以及上下游的参考序列
    query_1=get_fasta_sequence(ref_dict,chrname,pos-flank,pos-1)#变异上游的参考序列
    query_2=get_fasta_sequence(ref_dict,chrname,pos+len(ref_base),pos+flank-1+len(ref_base))#变异下游的参考序列
    query=complement_reverse(query_1+alt_base+query_2)#变异序列的反向互补序列，用于和参考序列做局部比对


    #局部比对
    (m_score,m_cigar,start1,end1,mseq1,start2,end2,mseq2,x,y,z,similarity)=SmithWaterman.local_alignment(refseq,query,0,flank+1)#cover snv的序列的反向互补比对到参考基因组序列
    #-------------------------------------------------------------------------------------------
    #(m_score,m_cigar,start1,end1,mseq1,start2,end2,mseq2,x,y,z)=SmithWaterman.local_alignment(refseq,query,0,0)
    #m_score 比对分数
    #m_cigar 比对详情
    #start1 取出的参考基因组正向序列比对起始位点
    #end1 取出的参考基因组正向序列比对终止位点
    #mseq1 比对到的参考基因组正向序列
    #start2 变异替换后的序列反向互补序列起始位点(正向的终止位点)
    #end2 变异替换后的序列反向互补序列终止位点(正向的起始位点)
    #mseq2 变异替换后的序列反向互补序列
    #-----------------------------------------------------------------------------


    #比对后分析
    if not int(m_score):#无回文假阳基础的位点
        return None

    #将比对的结果坐标映射到参考基因组上
    new_start1= pos - flank - 1 +int(start1)#参考基因组序列上实际的比对起始位点
    new_end1=pos - flank - 1 + int(end1)#参考基因组序列上实际的比对终止位点
    start2_r=2*flank+len(alt_base)-int(end2)+1#变异替换后的序列反向互补序列在正向序列的起始位点
    end2_r=2*flank+len(ref_base)-int(start2)+1#变异替换后的序列反向互补序列在正向序列的终止位点
    new_start2=pos-flank-1 +start2_r#变异替换后的序列反向互补序列在基因组上正向序列的起始位点
    new_end2=pos-flank-1+end2_r#变异替换后的序列反向互补序列在基因组上正向序列的终止位点
    reverse_query_seq=complement_reverse(mseq2)#变异替换后的序列的反向互补序列比对上的区域的反向互补序列（正向）
    return [chrname,pos,ref_base,alt_base,m_score,m_cigar,new_start1,new_end1,new_start2,mseq1,new_end2,reverse_query_seq,mseq2,similarity]

def read_ref_alignment(read,chrname,pos,ref_base,alt_base,ref_dict,flank,keysite):
    refseq=get_fasta_sequence(ref_dict,chrname,pos-flank,pos+flank-1+len(ref_base))#变异以及上下游的参考序列
    query=complement_reverse(read)#变异序列的反向互补序列，用于和参考序列做局部比对
    keysite=len(read)-keysite+1


    #局部比对
    (m_score,m_cigar,start1,end1,mseq1,start2,end2,mseq2,x,y,z,similarity)=SmithWaterman.local_alignment(refseq,query,0,keysite)#cover snv的序列的反向互补比对到参考基因组序列
    #-------------------------------------------------------------------------------------------
    #(m_score,m_cigar,start1,end1,mseq1,start2,end2,mseq2,x,y,z)=SmithWaterman.local_alignment(refseq,query,0,0)
    #m_score 比对分数
    #m_cigar 比对详情
    #start1 取出的参考基因组正向序列比对起始位点
    #end1 取出的参考基因组正向序列比对终止位点
    #mseq1 比对到的参考基因组正向序列
    #start2 变异替换后的序列反向互补序列起始位点(正向的终止位点)
    #end2 变异替换后的序列反向互补序列终止位点(正向的起始位点)
    #mseq2 变异替换后的序列反向互补序列
    #-----------------------------------------------------------------------------


    #比对后分析
    if not int(m_score):#无回文假阳基础的位点
        return None

    #将比对的结果坐标映射到参考基因组上
    new_start1= pos - flank - 1 +int(start1)#参考基因组序列上实际的比对起始位点
    new_end1=pos - flank - 1 + int(end1)#参考基因组序列上实际的比对终止位点
    new_start2=len(read)-int(end2)+1 # read反向互补之前的起始位点
    new_end2=len(read)-int(start2)+1 # read反向互补之前的终止位点
    reverse_mseq1=complement_reverse(mseq1)
    reverse_mseq2=complement_reverse(mseq2)
    return [chrname,pos,ref_base,alt_base,m_score,m_cigar,new_start1,new_end1,mseq1,reverse_mseq1,new_start2,new_end2,reverse_mseq2,mseq2,similarity]
def reverse_cigar(cigar):
    '''
    输出反向的cigar
    '''
    new_cigar=''
    pattern=re.compile('((\d+)([SDIM]))')#这里的S指的是subsititution，bam中的cigar没有碱基的替换
    cigar_parse=re.findall(pattern,cigar)
    cigar_parse.reverse()
    for each_cigar in cigar_parse:
        new_cigar+=each_cigar[0]
    return new_cigar

def present_artifact(match_list):
    '''
    输出比对结果，作为模块时不执行此函数，单独使用此脚本时使用，将比对结果输出到标准输出，并用下划线标记snv
    '''
    (chrname,pos,ref_base,alt_base,m_score,m_cigar,new_start1,new_end1,new_start2,mseq1,new_end2,reverse_query_seq,mseq2,similarity)=match_list
    #ref高亮
    part1=mseq1[0:pos-new_start1]
    part2=mseq1[pos-new_start1:pos-new_start1+len(ref_base)]
    part3=mseq1[pos-new_start1+len(ref_base):]
    new_mseq1=part1+'\033[4m'+part2+'\033[0m'+part3
    #ref reverse
    mseq1_reverse=complement_reverse(mseq1)
    part1=mseq1_reverse[0:pos-new_start1]
    part2=mseq1_reverse[pos-new_start1:pos-new_start1+len(alt_base)]
    part3=mseq1_reverse[pos-new_start1+len(alt_base):]
    new_mseq1_reverse=part1+'\033[4m'+part2+'\033[0m'+part3
    #query高亮
    part1=mseq2[0:pos-new_start2]
    part2=mseq2[pos-new_start2:pos-new_start2+len(alt_base)]
    part3=mseq2[pos-new_start2+len(alt_base):]
    new_mseq2=part1+'\033[4m'+part2+'\033[0m'+part3

    #  ref cover snv
    ref_cover_snv_tmp2=get_fasta_sequence(ref_dict,chrname,new_start2,new_end2)
    part1=ref_cover_snv_tmp2[0:pos-new_start2]
    part2=ref_cover_snv_tmp2[pos-new_start2:pos-new_start2+len(ref_base)]
    part3=ref_cover_snv_tmp2[pos-new_start2+len(ref_base):]
    ref_cover_snv=part1+'\033[4m'+part2+'\033[0m'+part3

    #用alt替代ref型，alt cover snv
    part2=alt_base
    alt_cover_snv=part1+'\033[4m'+part2+'\033[0m'+part3
    true_seq=part1+part2+part3
    new_cigar=reverse_cigar(m_cigar)


    print('\t'.join([chrname,str(pos),ref_base,alt_base,str(m_score),m_cigar,new_cigar,str(similarity),str(new_start2),str(new_end2),ref_cover_snv,alt_cover_snv,str(new_start1),str(new_end1),mseq1,mseq1_reverse]))

if __name__=='__main__':
    argv=argument_parser()
    snv=argv['snv']
    fasta=argv['fasta']
    outdir=argv['outdir']
    outflie=argv['outfile']
    flank=argv['flank']
    present=argv['present']
#   读fasta
    ref_dict=read_fasta(fasta)

#   处理snv，是一个位点还是一个文件
    match_list=[]
    if ':' in snv:
        match_result=structure_identify(snv,ref_dict)
        if match_result:
            match_list.append(match_result)
    elif os.path.exists(snv):
        with open(snv) as r:
            header=r.readline().strip("\n")
            header_info=header.split('\t')
            while True:
                line=r.readline().strip("\n")
                if not line:break
                lineInfo=line.split('\t')
                site=':'.join([lineInfo[SR.snv_header_index(header_info,'chrom')],lineInfo[SR.snv_header_index(header_info,'pos_raw')],lineInfo[SR.snv_header_index(header_info,'ref')],lineInfo[SR.snv_header_index(header_info,'alt')]])
                match_result=structure_identify(site,ref_dict)
                if match_result:
                    match_list.append(match_result)
    else:
        print("Wrong parameter of --snv")
        exit()
    if present:
        for eachsnv in match_list:
            present_artifact(eachsnv)

