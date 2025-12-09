#!coding:utf-8

#----------------------------------------------
#Project: 商检-基础流程-基础分析
#Description: 酶切过滤程序的核心模块-比对部分, 也可单独使用，输出snv附近的负链比对结果
#Usage: ./artifact_identify.py -h
#Author: 骆磊
#------------------------------------------------

#-----------------
#v0.4, support_reads从v2修改为v3
#v0.5, support_reads_v5, reads分级
import pysam
import argparse
import re
import os
import bamhunter
# from multi_thread import *  # 调用的脚本函数有问题，不会用到
from artifact_identify import *
import Needleman_Wunsch as NW
import support_reads as sr
import time



def argument_parser():
    parser = argparse.ArgumentParser(description="input bam file and snv file to find out which read support this snv")
    parser.add_argument('--bamfile',help='input bam file',required=True)
    #/beegfs/work/commercial_test/cupcake/databases/gatk_bundle/2.8/hg19/ucsc.hg19.noconfig.fasta
    parser.add_argument('--fasta',help='fasta file used',required=True)
    parser.add_argument('--outdir',help='outdir, default current directory',default='./')
    parser.add_argument('--outfile',help='prefix of outfie',required=True)
    parser.add_argument('--snv',help='snv file',required=True)
    parser.add_argument('--flank',help='flank of snv to match',type=int,default=300)
    argv=vars(parser.parse_args())
    return argv

def find_passengers(seq1,seq2,seq3):
    '''
    寻找是否匹配会引起passengers snv，并且read是否支持passengers snv
    '''
    match_cigar_1_2=NW.local_alignment(seq2,seq1)
    match_cigar_1_3=NW.local_alignment(seq3,seq1)

    pattern=re.compile('((\d+)([SHMXNDI]))')
    cigar_parse_1_2=re.findall(pattern,match_cigar_1_2)
    cigar_parse_1_3=re.findall(pattern,match_cigar_1_3)

    query_pos=0
    ref_pos=0
    cigar_status_1_2=[]
    for each_cigar in cigar_parse_1_2:
        if each_cigar[2]=='I':#消耗reads，ref不变
            cigar_status_1_2.append([query_pos,seq2[query_pos:query_pos+int(each_cigar[1])],ref_pos,'-'])#ref是1based
            query_pos+=int(each_cigar[1])
            #query_pos和query_pos+int(each[1])是ins在reads上的起始和终止位点
            #ref_pos+1是插入位置, 应该是在ref_pos和ref_pos+1中间插入
            #"-",表示ref型是'-'
            #-----------------------
        #D会消耗ref的碱基，但reads的碱基位置不变，所以只有ref的指针需要变化
        if each_cigar[2]=='D':#消耗ref，reads不变
            cigar_status_1_2.append([query_pos,'-',ref_pos,seq1[ref_pos:ref_pos+int(each_cigar[1])]])
            ref_pos+=int(each_cigar[1])
        #mismatch
        if each_cigar[2]=='S':
            for i in range(int(each_cigar[1])):
                cigar_status_1_2.append([query_pos,seq2[query_pos],ref_pos,seq1[ref_pos]])
                ref_pos+=1
                query_pos+=1
        if each_cigar[2]=='M':
            for i in range(int(each_cigar[1])):
                ref_pos+=1
                query_pos+=1

    query_pos=0
    ref_pos=0
    cigar_status_1_3=[]
    for each_cigar in cigar_parse_1_3:
        if each_cigar[2]=='I':#消耗reads，ref不变
            cigar_status_1_3.append([query_pos,seq3[query_pos:query_pos+int(each_cigar[1])],ref_pos,'-'])#ref是1based
            query_pos+=int(each_cigar[1])
            #query_pos和query_pos+int(each[1])是ins在reads上的起始和终止位点
            #ref_pos+1是插入位置, 应该是在ref_pos和ref_pos+1中间插入
            #"-",表示ref型是'-'
            #-----------------------
        #D会消耗ref的碱基，但reads的碱基位置不变，所以只有ref的指针需要变化
        if each_cigar[2]=='D':#消耗ref，reads不变
            cigar_status_1_3.append([query_pos,'-',ref_pos,seq1[ref_pos:ref_pos+int(each_cigar[1])]])
            ref_pos+=int(each_cigar[1])
        #match or mismatch
        if each_cigar[2]=='S':
            for i in range(int(each_cigar[1])):
                cigar_status_1_3.append([query_pos,seq3[query_pos],ref_pos,seq1[ref_pos]])
                ref_pos+=1
                query_pos+=1
        if each_cigar[2]=='M':
            for i in range(int(each_cigar[1])):
                ref_pos+=1
                query_pos+=1

    passengers=[]
    for eachvariation in cigar_status_1_3:
        if eachvariation in cigar_status_1_2:
            passengers.append(eachvariation)

    return passengers


def get_complement_region(ref_sequence,ref_pos,sequence,cigar,match_results_f):
    '''
    将read匹配上的区域对应到基因组上
    '''
    (chrname,pos,ref_base,alt_base,m_score,m_cigar,new_start1,new_end1,mseq1,reverse_mseq1,start2,end2,mseq2,reverse_mseq2,similarity)=match_results_f
    m_cigar=reverse_cigar(m_cigar)
    start2-=1
    end2-=1
    ref_pos-=1
    new_start1-=1
    new_end1-=1
    query_pos=0#作为指针，用于标记当前处理的碱基相对于reads起始碱基的位置
    pattern=re.compile('((\d+)([SHMXNDI]))')#提取cigar列信息的pattern
    cigar_parse=re.findall(pattern,cigar)
    read_map_status=[]
    read_complement_pos=[]#存储read的complement的位点

    for each_cigar in cigar_parse:#cigar中每个元素的处理, 每个子列表的第一个元素是不需要的[['',10,'S'],['',33,'M'],['',2,'I'],['',15,'M'],['',6,'D'],['',29,'M'],['',1,'D'],['',11,'M']]

        #S部分, reads序列会在bam中展示,起始位点其实是在reads中间的某个碱基
        if each_cigar[2]=='S':#跳过，不修改任何东西，但reads的坐标需要改动,ref坐标不变
            query_pos+=int(each_cigar[1])
        #H对解析bam无影响，bam中不会显示H部分的序列
        if each_cigar[2]=='H': #直接跳过
            pass
        #I会消耗reads的碱基，但ref的位置不变，所以只有reads的指针需要变化
        if each_cigar[2]=='I':#消耗reads，ref不变
            read_map_status.append([query_pos,sequence[query_pos:query_pos+int(each_cigar[1])],ref_pos+1,"-"])
            query_pos+=int(each_cigar[1])
            #query_pos和query_pos+int(each[1])是ins在reads上的起始和终止位点
            #ref_pos+1是插入位置, 应该是在ref_pos和ref_pos+1中间插入
            #"-",表示ref型是'-'
            #-----------------------
        #D会消耗ref的碱基，但reads的碱基位置不变，所以只有ref的指针需要变化
        if each_cigar[2]=='D':#消耗ref，reads不变
            read_map_status.append([query_pos,"-",ref_pos+1,ref_sequence[ref_pos:ref_pos+int(each_cigar[1])]])
            ref_pos+=int(each_cigar[1])

        #match时，ref和alt的指针都需要变化，这里先不管mismatch的情况
        if each_cigar[2]=='M':
            for i in range(int(each_cigar[1])):
                read_map_status.append([query_pos,sequence[query_pos],ref_pos+1,ref_sequence[ref_pos]])
                ref_pos+=1
                query_pos+=1
    k=0
    for each_base_map in read_map_status:
        #I会消耗reads的碱基，但ref的位置不变，所以只有reads的指针需要变化
        if each_base_map[0] >=start2 and each_base_map[0]<=end2:
            read_complement_pos.append(each_base_map[0])

            if k==0:
                complement_start=each_base_map[2]
                complement_end=each_base_map[2]
                k=1
            else:
                complement_end=each_base_map[2]


    try:

        read_complement_start=min(read_complement_pos)# 0-based
        read_complement_end=max(read_complement_pos)# 0-based
    except:
        return 0,0,0,0,0,0,0,0,0
    #local alignment cigar
    read_complement_cigar_status=[]
    ref_seq_c=reverse_mseq1# 匹配部分的基因组序列的反向互补序列
    ref_start_c=new_end1# 由于反向互补，起始位置其实是由终止位置来的
    ref_end_c=new_start1
    query_pos_c=start2

    m_cigar_parse=re.findall(pattern,m_cigar)
    for each_m_cigar in m_cigar_parse:
        if each_m_cigar[2]=='I':#消耗reads，ref不变
            read_complement_cigar_status.append([query_pos_c,ref_start_c+1])#ref是1based
            query_pos_c+=int(each_m_cigar[1])
            #query_pos和query_pos+int(each[1])是ins在reads上的起始和终止位点
            #ref_pos+1是插入位置, 应该是在ref_pos和ref_pos+1中间插入
            #"-",表示ref型是'-'
            #-----------------------
        #D会消耗ref的碱基，但reads的碱基位置不变，所以只有ref的指针需要变化
        if each_m_cigar[2]=='D':#消耗ref，reads不变
            read_complement_cigar_status.append([query_pos_c,ref_start_c+1])
            ref_start_c-=int(each_m_cigar[1])
        #match or mismatch
        if each_m_cigar[2]=='M' or each_m_cigar[2]=='S':
            for i in range(int(each_m_cigar[1])):
                read_complement_cigar_status.append([query_pos_c,ref_start_c+1])
                ref_start_c-=1
                query_pos_c+=1
    tmp_read_complement_cigar_status=[]
    for each_base in read_complement_cigar_status:
        if each_base[0] < read_complement_start or each_base[0] >read_complement_end:
            continue
        tmp_read_complement_cigar_status.append(each_base)
    ref_max=tmp_read_complement_cigar_status[0][1]
    ref_min=tmp_read_complement_cigar_status[-1][1]
    ref_alt_complement_seq=ref_sequence[ref_min-1:ref_max]

    try:
        return complement_start,complement_end,ref_sequence[complement_start-1:complement_end],ref_min,ref_max,ref_alt_complement_seq,read_complement_start+1,read_complement_end+1,sequence[read_complement_start:read_complement_end+1]
    except:
        print(start2)
        print(end2)
        print(ref_pos)
        print(read_map_status)
        return complement_start,complement_end





def filter_snv(snv,snv_information_dict):
    '''
    这个函数用于去掉认为不是回文阳性的snv, 频率>30或者是黑名单的位点不做处理
    '''
    freq=snv_information_dict[snv][0]
    black_flag=snv_information_dict[snv][1]
    if not freq:
        freq=0
    if freq > 30:# freq over 30 will not be analysed
        return False
#    if 'homologouse_black' in black_flag or 'repeat_black' in black_flag or 'site_black' in black_flag:
    if "black_site" == black_flag:
        return False

    return True

def read_vcf_with_support_reads(snv):
    '''
    读入带有support reads信息的snv文件，并存入字典
    '''
    all_support_reads_dict={}#{readid:None,...}
    snv_support_reads_dict={}
    with open(snv) as snv_handle:
        header=snv_handle.readline().strip('\n')
        #get black flag
        header_info=header.split('\t')
        while True:
            line=snv_handle.readline()
            if not line:break
            lineInfo=line.split('\t')
            key='\t'.join([lineInfo[sr.snv_header_index(header_info,'chrom')],lineInfo[sr.snv_header_index(header_info,'pos_raw')],lineInfo[sr.snv_header_index(header_info,'ref')],lineInfo[sr.snv_header_index(header_info,'alt')]])
            value=lineInfo[-3].split(',')
#            for eachread in value:
#                all_support_reads_dict[eachread]=None
            snv_support_reads_dict[key]=value
#    return header,snv_support_reads_dict,all_support_reads_dict
    return header,snv_support_reads_dict

def change_sequence_by_snv(start,end,region,snv):#1-based
    '''
    比对前，用alt带替ref序列
    '''
    (chrname,snv_pos,ref_base,alt_base)=snv.split('\t')
    snv_pos=int(snv_pos)
    snv_len=len(ref_base)
    end=end-start
    snv_pos=snv_pos-start
    start=0
    left_seq=region[0:snv_pos]
    wait_seq=region[snv_pos:snv_pos+snv_len]
    right_seq=region[snv_pos+snv_len:]
    return left_seq+alt_base+right_seq

#def read_bam_support_reads(bamfile_handle,supportreads_dict):
#    readid_sequence_dict={}
#    read_ref_pair_pos_dict={}
#    for read in bamfile_handle:
#        readid=read.query_name+'_'+str(read.flag)
#        if not  readid in supportreads_dict:continue
#        cigar=read.cigarstring


def get_position_read_ref(read_ref_pair_pos_dict,readid,pos):
    '''
    从reads的碱基位点推断基因组上的位点
    '''
    for eachpos in read_ref_pair_pos_dict[readid]:
        if eachpos[0]==pos-1:
            return eachpos[1]+1
    raise Exception()

def get_position_ref_read(read_ref_pair_pos_dict,readid,pos):
    '''
    从ref位点推断reads位点
    '''
    aligned_pairs=read_ref_pair_pos_dict[readid][3].get_aligned_pairs()
    for eachpos in aligned_pairs:
        if eachpos[1]==pos-1:
        #chr15    90630344    C    G    15    E100044124L1C030R01000347883_163    90630333    4S94M2S    40    20M    1.0    90630334    90630353    CCCAGCGTACCCTGGGCCAG    CTGGCCCAGGGTACGCTGGG    1    20    CTGGCCCAGGGTACGCTGGG    CCCAGCGTACCCTGGGCCAG    CTGGCCCAGGGTACGCTGGGCCAGGATGTCTGACTGCACATCTCCGTCATAGTTCTTGCAGGCCCACACAAAGCCACCCGAAGACTTGAGGACCTGAGAT
#[(0, 90631630), (1, 90631631), (2, 90631632), (3, 90631633), (4, 90631634), (5, 90631635), ******(None, 90631636)******, (6, 90631637), (7, 90631638), (8, 90631639), (9, 90631640), (10, 90631641), (11, None), (12, 90631642), (13, 90631643), (14, 90631644), (15, 90631645), (16, 90631646), (17, 90631647), (18, 90631648), (19, 90631649), (20, 90631650), (21, 90631651), (22, 90631652), (23, 90631653), (24, 90631654), (25, 90631655), (26, 90631656), (27, 90631657),(28, 90631658), (29, 90631659), (30, 90631660), (31, 90631661), (32, 90631662), (33, 90631663), (34, 90631664), (35, 90631665), (36, 90631666), (37, 90631667), (38, 90631668), (39, 90631669), (40, 90631670), (41, 90631671), (42, 90631672), (43, 90631673), (44, 90631674), (45, 90631675), (46, 90631676), (47, 90631677), (48, 90631678), (49, 90631679), (50, 90631680), (51, 90631681), (52, 90631682), (53, 90631683), (54, 90631684), (55, 90631685), (56, 90631686), (57, 90631687), (58, 90631688), (59, 90631689), (60, 90631690), (61, 90631691), (62, 90631692), (63, 90631693), (64, 90631694), (65, 90631695), (66, 90631696), (67, 90631697), (68, 90631698), (69, 90631699), (70, 90631700), (71, 90631701), (72, 90631702), (73, 90631703), (74, 90631704), (75, 90631705), (76, 90631706), (77, 90631707), (78, 90631708), (79, 90631709), (80, 90631710), (81, 90631711), (82, 90631712), (83, 90631713), (84, 90631714), (85, 90631715), (86, 90631716), (87, 90631717), (88, 90631718)] 90631636发生的Deletion
            if eachpos[0] or eachpos[0] == 0:#如果是del，发生del的read pos是none
                return eachpos[0]+1
            else:
                return last_pos+1
        last_pos=eachpos[0]
        #insert在read边缘的情况，返回1
    if aligned_pairs[0][1]==None and aligned_pairs[0][0]==0:
        return 1
    print(aligned_pairs)
    print(readid)
    print(pos)
    raise Exception()
def get_refseq(chrname,pos,ref,alt,ref_dict,flank):
    '''
    获取snv附近的参考基因组序列
    '''
    upstream=ref_dict[chrname][pos-flank:pos]
    downstream=ref_dict[chrname][pos+len(ref):pos+len(ref)+flank]
    refseq=upstream+ref+downstream
    return refseq

def reverse_cigar(cigarstring):
    pattern=re.compile('((\d+)([SHMXNDI]))')
    cigar_parse=re.findall(pattern,cigarstring)
    cigar_parse.reverse()
    new_cigar=""
    for each_cigar in cigar_parse:
        new_cigar+=each_cigar[0]
    return new_cigar

def main(bam,fasta,outdir,outfile,snvfile,flank):
    out_prefix=os.path.join(outdir,outfile)
    final_outfile=out_prefix+'.snv_supportReads.txt'
    #-------------------初始准备
    #读取fasta
#    print("Read fasta")
    ref_dict=read_fasta(fasta)
#    print("Read fasta over")
    out_prefix=os.path.join(outdir,outfile)


    #处理bam获取每个snv的support reads，这步可以跳过以用于兼容不同方式获取的support reads
    bamfile=bamhunter.bamhunter(bam)
    (snv_dict,snv_indel_dict,snv_complex_dict,raw_snv_dict,raw_snv_list,snv_header,snv_info_dict)=sr.parser_snv_file(snvfile)
    snv_support_dict=sr.single_base_snv(bamfile,snv_dict,ref_dict)
    snv_indel_support_dict=sr.indel_snv(bamfile,snv_indel_dict,ref_dict)
    snv_complex_support_dict=sr.complex_snv(bamfile,snv_complex_dict,ref_dict)
    w=open(final_outfile,'w')
    w.write(snv_header+'\t'+'\t'.join(['support_reads_original','support_reads_num_nondedup','support_reads_num_dedup'])+'\n')
    for eachsnv in raw_snv_list:
        if eachsnv in snv_support_dict:
            supportreads=','.join(snv_support_dict[eachsnv][0])
            count_dup=str(len(snv_support_dict[eachsnv][0]))
            count_nondup=str(len(snv_support_dict[eachsnv][1]))
        elif eachsnv in snv_indel_support_dict:
            supportreads=','.join(snv_indel_support_dict[eachsnv][0])
            count_dup=str(len(snv_indel_support_dict[eachsnv][0]))
            count_nondup=str(len(snv_indel_support_dict[eachsnv][1]))
        elif eachsnv in snv_complex_support_dict:
            supportreads=','.join(snv_complex_support_dict[eachsnv][0])
            count_dup=str(len(snv_complex_support_dict[eachsnv][0]))
            count_nondup=str(len(snv_complex_support_dict[eachsnv][1]))
        else:
            supportreads='No_support_reads'
            count_dup='0'
            count_nondup='0'

        w.write(raw_snv_dict[eachsnv]+'\t'+'\t'.join([supportreads,count_dup,count_nondup])+'\n')
    w.close()
    del snv_support_dict
    del snv_indel_support_dict
    del snv_complex_support_dict

    print("Reading vcf...")
    #(header,snv_support_reads_dict,supportreads_dict)=read_vcf_with_support_reads(final_outfile)
    (header,snv_support_reads_dict)=read_vcf_with_support_reads(final_outfile)
    time.sleep(10)
    print("Reading vcf over")


    alignment_results=out_prefix+'_artifact_match.txt'
    w=open(alignment_results,'w')
    for eachsnv in raw_snv_list:
        if filter_snv(eachsnv,snv_info_dict):#根据snv的信息预判断是否为回文假阳，如果不是，那就不参与后面的分析，直接判为阴性
            pass
        else:
            continue
        reaidlist=snv_support_reads_dict[eachsnv]
        if not reaidlist:continue
        (chrname,pos,ref,alt)=eachsnv.split('\t')
        for read in reaidlist:
            if not read:continue
            cigar=sr.READID_SEQUENCE_DICT[read][0]
            read_start=int(sr.READID_SEQUENCE_DICT[read][2])
            snv_pos_in_read=get_position_ref_read(sr.READID_SEQUENCE_DICT,read,int(pos))
            match_results=read_ref_alignment(sr.READID_SEQUENCE_DICT[read][3].query_sequence,chrname,int(pos),ref,alt,ref_dict,flank,snv_pos_in_read)
            if match_results:
                chrname,pos,ref_base,alt_base,m_score,m_cigar,new_start1,new_end1,mseq1,reverse_mseq1,start2,end2,mseq2,reverse_mseq2,similarity=match_results
                new_m_cigar=reverse_cigar(m_cigar)
                (complement_raw_start,complement_raw_end,complement_raw_region,complement_mate_start,complement_mate_end,complement_mate_region,complement_read_start,complement_read_end,complement_read_region)=get_complement_region(ref_dict[chrname],read_start,sr.READID_SEQUENCE_DICT[read][3].query_sequence,cigar,match_results)
                if not complement_raw_region:
                    continue
                complement_raw_region_modify=change_sequence_by_snv(complement_raw_start,complement_raw_end,complement_raw_region,eachsnv)#变异应该都必定存在complement_raw_region范围内，如果不是，这里会报错，应该不存在这种情
                passengers=find_passengers(complement_raw_region_modify,complement_reverse(complement_mate_region),complement_read_region)#[[8, 'C', 8, 'T'], [11, 'A', 11, 'G'], [20, 'T', 20, 'G']]
                really_passengers=[]
                for eachpassenger in passengers:
                    fix_pos=eachpassenger[0]+complement_raw_start
                    really_passengers.append(str(fix_pos)+':'+eachpassenger[3]+':'+eachpassenger[1])
                w.write('\t'.join([chrname,str(pos),ref_base,alt_base,str(snv_pos_in_read),read,str(read_start),cigar,str(m_score),str(new_m_cigar),str(similarity),str(new_start1),str(new_end1),mseq1,reverse_mseq1,str(start2),str(end2),mseq2,reverse_mseq2,str(complement_raw_start),str(complement_raw_end),complement_raw_region,str(complement_mate_start),str(complement_mate_end),complement_mate_region,complement_reverse(complement_mate_region),str(complement_read_start),str(complement_read_end),complement_read_region,sr.READID_SEQUENCE_DICT[read][-1],','.join(really_passengers),sr.READID_SEQUENCE_DICT[read][3].query_sequence])+'\n')
#                w.write('\t'.join([chrname,str(pos),ref_base,alt_base,str(snv_pos_in_read),read,str(read_start),cigar,str(m_score),str(m_cigar),str(similarity),str(new_start1),str(new_end1),mseq1,reverse_mseq1,str(start2),str(end2),mseq2,reverse_mseq2,READID_SEQUENCE_DICT[read][3]])+'\n')

if __name__=='__main__':
    main()
