#!coding:utf-8
#------------------------------------------------
#v3: 新增过滤规则，向商检靠拢
#v4.3: rm single snv read, lowq, confulict
#------------------------------------
#Project: 商检-基础流程-基础分析
#Description: 根据snv文件和bam文件，对每个snv的support reads进行回溯,输出每个snv的support reads的ID。 support_reads_original:support reads ID; support_reads_num_nondedup: 回溯得到的support reads数，read1和read2同时支持时算两条; support_reads_num_dedup:回溯得到的support reads数，read1和read2同时支持时算一条
#Usage: python3 ./support_reads.py -h
#Author: 骆磊
import pysam
import argparse
import re
import os
import bamhunter
import time
import random


#全局变量
global READID_SEQUENCE_DICT
READID_SEQUENCE_DICT={}

def argument_parser():
    parser = argparse.ArgumentParser(description="input bam file and snv file to find out which read support this snv")
    parser.add_argument('--bamfile',help='input bam file',required=True)
    #/beegfs/work/commercial_test/cupcake/databases/gatk_bundle/2.8/hg19/ucsc.hg19.noconfig.fasta
    parser.add_argument('--fasta',help='fasta file used',required=True)
    parser.add_argument('--outdir',help='outdir, default current directory',default='./')
    parser.add_argument('--outfile',help='prefix of outfie',required=True)
    parser.add_argument('--snv',help='snv file',required=True)
    argv=vars(parser.parse_args())
    return argv

#读取fasta文件，返回以染色体为键，序列为值的字典中,同时读取的任务太多，可能会引起内存错误,为了避免，这里加入了循环,如果失败就停一会儿再重读，没有设置循环终止，因为正常情况下都能正常读取

def parser_alignment_direction(readhandle):
    if readhandle.is_secondary or readhandle.is_supplementary or readhandle.is_unmapped or readhandle.is_duplicate:
        return False
    chr1=readhandle.reference_id    
    chr2=readhandle.next_reference_id
    if chr1!=chr2:
        samechr=False
    else:
        samechr=True
    #read 是read1的情况:
    if readhandle.is_read1:
        if not samechr:#map到不同染色体，这种时候只能关注单链的方向了
            if readhandle.is_reverse:#R1
                return "R1_diffchr"
            else:
                return "F1_diffchr"
        if readhandle.mate_is_unmapped:#mate unmapped的情况
            if readhandle.is_reverse:#R1
                return "R1_unmap"
            else:
                return "F1_unmap"
        if True:#properly + unusual insert size
            if readhandle.is_reverse: # R1的情况
                if readhandle.mate_is_reverse:#R2
                    return "R1R2"
                else:
                    return "F2R1"
            else:#F1
                if readhandle.mate_is_reverse:#R2
                    return "F1R2"
                else:#F2
                    return "F1F2"
        raise Exception
    else:#当read是read2
        if not samechr:#map到不同染色体，这种时候只能关注单链的方向了
            if readhandle.is_reverse:#R2
                return "R2_diffchr"
            else:
                return "F2_diffchr"
        if readhandle.mate_is_unmapped:#mate unmapped的情况
            if readhandle.is_reverse:#R2
                return "R2_unmap"
            else:
                return "F2_unmap"
#        if readhandle.is_proper_pair or readhandle.template_length <1000:# 正常比对情况
        if True:#properly + unusual insert size
            if readhandle.is_reverse: # R2的情况
                if readhandle.mate_is_reverse:#R1
                    return "R1R2"
                else:#F1
                    return "F1R2"
            else:#F2
                if readhandle.mate_is_reverse:#R1
                    return "F2R1"
                else:#F1
                    return "F1F2"
        raise Exception
        
                
def read_fasta(fasta):
    k_count=0
    while True:
        try:
            ref_dict={}#{chr:sequence,...}
            fasta_handle=pysam.FastxFile(fasta)
            for each_chr in fasta_handle:
                ref_dict[each_chr.name]=each_chr.sequence.upper()
            return ref_dict
        except:
            k_count+=1
            print("error while reading fasta, retry "+k_count+' time')
            time.sleep(random.randint(20,60))

def read_cover_single_snv(pos,read_map_ref):#for single snv
    for eachitem in read_map_ref[1:-1]:
        if eachitem[2]==pos:
            try:
                if eachitem[-1]>=15:
                    return True
                else:
                    return False
            except TypeError:
                pass
    return False
def read_cover_LR_single_snv(pos,read_map_ref):#for single snv flank
    for eachitem in read_map_ref:
        if eachitem[2]==pos:
            if eachitem[1]==eachitem[3]:
                return True
            else:
                try:
                    if eachitem[-1]<15:
                        return True
                except:
                    pass
    return False
                
def read_cover_flank_snv(pos,read_map_ref):#for indel and complex snv
    for eachitem in read_map_ref:
        if eachitem[2]==pos:
            if eachitem[1]==eachitem[3]:
                return True
    return False

def check_reads(bamfile,readid):
    for eachread in bamfile.pysam_bam:
        if eachread.query_name==readid:
            t=[eachread.query_name,str(eachread.flag),str(eachread.reference_start),eachread.cigarstring]
            print(t)
            print(eachread.get_aligned_pairs(with_seq=True))



def single_base_snv(bamfile,snv_dict,ref_dict):
    snv_support_dict={}
    
    pattern=re.compile('((\d+)([SHMXNDI]))')#提取cigar列信息的pattern
    for eachsnv in list(snv_dict.keys()):
        (chrname,pos,refbase,altbase)=snv_dict[eachsnv]
        might_support_reads=bamfile.get_pos_reads(chrname,pos-1,pos+1)#获取覆盖这个位点的所有reads，但受限于pysam，最小只能找2个bp的覆盖reads，所以这些reads不一定真正比对上了，另外，如果这个位点是缺失，那也会被判断为没有比对上，对结果无影响，暂时不修改;如果查询的区域，是mismatch，也会被判断为没有比对上，所以要找一段范围内cover的reads，暂定为50，毕竟连续50个mismatch可能性非常低。
        #---------------------------------------------
        #为什么要找覆盖pos-1和pos+1的reads
        #第一个原因是pysam只能找cover一个“区域”的reads，最少2bp
        #第二个原因是不匹配的位点，不算cover，比如一个reads在100这个位点是mismatch，那这个reads不会被提取出来，所以如果只找cover snv单个位点，mismatch的reads都会被丢掉。
        #所以要找snv附近一个bp被覆盖到的，或的关系，但是如果前后都不是match，insertion，或者softclip或者read结尾，也找不到，所以这里认为snv前后的一个bp，必须是match。理应是这样，不然这个snv会被检测为其它变异，而不是一个单核苷酸突变,下面indel和complex突变的处理方式相同。
        #---------------------------------------------
        #要找的实际是read和ref在这个位置都有碱基，但碱基不一样
        #如果是个多碱基突变，支持多碱基突变的reads还是会支持单碱基突变
        readid_list={}
        simple_readid_list={}#support readid 不带flag
        cover_snv_readlist={}#cover这个位点的fragment列表
        simple_readid_readid={}#{simple_readid:readid}
        for eachread in might_support_reads:
            readid=eachread.query_name+'_'+str(eachread.flag)
            simple_readid=eachread.query_name
            cigar=eachread.cigarstring
            if not cigar:continue
            #解析cigar
            ref_pos=eachread.reference_start # 和ref_start初始值是一样的，但ref_pos后续用于ref的指针
            ref_start=eachread.reference_start
            sequence=eachread.query_sequence# 含softclip的序列
            m_chrname=bamfile.pysam_bam.get_reference_name(eachread.reference_id)
            try:
                assert chrname == m_chrname
            except:
                print("some trouble with read "+simple_readid+" ,maybe there is secondary aligned for this read?")
            ref_sequence=ref_dict[chrname]
            cigar_parse=re.findall(pattern,cigar)
            query_quality=eachread.query_qualities
            query_pos=0#作为指针，用于标记当前处理的碱基相对于reads起始碱基的位置
            read_map_status=[]
            for each_cigar in cigar_parse:#cigar中每个元素的处理, 每个子列表的第一个元素是不需要的[['',10,'S'],['',33,'M'],['',2,'I'],['',15,'M'],['',6,'D'],['',29,'M'],['',1,'D'],['',11,'M']]
            #S部分, reads序列会在bam中展示,起始位点其实是在reads中间的某个碱基
                if each_cigar[2]=='S':#跳过，不修改任何东西，但reads的坐标需要改动,ref坐标不变
                    query_pos+=int(each_cigar[1])
                    #----------
                    #5S10M
                    #指针位置在第一个碱基:[A]AAAATTTTTTTTTT
                    #修改指针后: AAAAA[T]TTTTTTTTT
                    #-----------
                #H对解析bam无影响，bam中不会显示H部分的序列
                if each_cigar[2]=='H': #直接跳过
                    pass
                #I会消耗reads的碱基，但ref的位置不变，所以只有reads的指针需要变化
                if each_cigar[2]=='I':#消耗reads，ref不变
                    read_map_status.append([query_pos,sequence[query_pos:query_pos+int(each_cigar[1])],ref_pos+1,"-",query_quality[query_pos-1:query_pos+int(each_cigar[1])]])
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
                        read_map_status.append([query_pos,sequence[query_pos],ref_pos+1,ref_sequence[ref_pos],query_quality[query_pos]])
                        ref_pos+=1
                        query_pos+=1
            #解析cigar结束
            #从cigar判断read有没有cover snv
            readCover=read_cover_single_snv(pos,read_map_status)
            if not readCover:continue#read没有cover到snv位点或者snv位点质量太低就再见
#            if readid=="E100047980L1C004R00203310736_163":print(read_map_status)
#            readCover_left=read_cover_flank_snv(pos-1,read_map_status)
#            readCover_right=read_cover_flank_snv(pos+1,read_map_status)
#            if not readCover_left and readCover_right:
#                continue
            fragment_direction=parser_alignment_direction(eachread)
            if not fragment_direction:
                continue
            #从解析完cigar的结果中提取所有
            key=0# 不support
            for each_position in read_map_status:#[35, '-', 151856125, 'A'], 0是0-based，2是1-based
#                if each_position[3]=="-": #insertion
#                elif each_position[1]=='-':#deletion
#                elif each_position[1]!=each_position[3]:#sn 
                
                if pos == each_position[2] and refbase==each_position[3] and altbase==each_position[1]:# single base snv support reads
                    readCover_left=read_cover_LR_single_snv(pos-1,read_map_status)
                    readCover_right=read_cover_LR_single_snv(pos+1,read_map_status)
                    if not (readCover_left and readCover_right):
                        continue
                    add_item_dict(eachread,chrname,fragment_direction)
                    if simple_readid in cover_snv_readlist and simple_readid not in simple_readid_list:#mate cover了，但是不支持变异
                        break
                    assert readid not in readid_list
                    readid_list[readid]=None
                    if simple_readid not in simple_readid_list:
                        simple_readid_list[simple_readid]=None
                    key=1#support
                    simple_readid_readid[simple_readid]=readid
            if key==0:
                if  simple_readid in simple_readid_list:#在这个read之前已经把materead加进去了,当前read不支持变异
                    del simple_readid_list[simple_readid]
                    del readid_list[simple_readid_readid[simple_readid]]

            cover_snv_readlist[simple_readid]=None
        snv_support_dict[eachsnv]=[list(readid_list.keys()),list(simple_readid_list.keys())]
    return snv_support_dict


def indel_snv(bamfile,snv_indel_dict,ref_dict):

    #处理indel
    snv_indel_support={}
    pattern=re.compile('((\d+)([SHMXNDI]))')#提取cigar列信息的pattern
    for eachsnv in list(snv_indel_dict.keys()):
        readid_list={}#support readid 带flag 
        simple_readid_list={}#support readid 不带flag
        (chrname,pos,refbase,altbase)=snv_indel_dict[eachsnv]
        might_support_reads=bamfile.get_pos_reads(chrname,pos-1,pos+len(refbase))#ins和del还不一样，inscover了3个bp，不过区别不大，忽略
        for eachread in might_support_reads:
            readid=eachread.query_name+'_'+str(eachread.flag)
            simple_readid=eachread.query_name
            cigar=eachread.cigarstring
            if not cigar:continue
            if 'H' in cigar:continue
            #解析cigar
            ref_pos=eachread.reference_start # 和ref_start初始值是一样的，但ref_pos后续用于ref的指针
            ref_start=eachread.reference_start
            sequence=eachread.query_sequence# 含softclip的序列
            query_quality=eachread.query_qualities
            assert len(query_quality)==len(sequence)
            m_chrname=bamfile.pysam_bam.get_reference_name(eachread.reference_id)
            try:
                assert chrname == m_chrname
            except:
                print("some trouble with read "+simple_readid+" ,maybe there is secondary aligned for this read?")
            ref_sequence=ref_dict[chrname]
            cigar_parse=re.findall(pattern,cigar)
            query_pos=0#作为指针，用于标记当前处理的碱基相对于reads起始碱基的位置
            read_map_status=[]
            for each_cigar in cigar_parse:#cigar中每个元素的处理, 每个子列表的第一个元素是不需要的[['',10,'S'],['',33,'M'],['',2,'I'],['',15,'M'],['',6,'D'],['',29,'M'],['',1,'D'],['',11,'M']]
            #S部分, reads序列会在bam中展示,起始位点其实是在reads中间的某个碱基
                if each_cigar[2]=='S':#跳过，不修改任何东西，但reads的坐标需要改动,ref坐标不变
                    query_pos+=int(each_cigar[1])
                    #----------
                    #5S10M
                    #指针位置在第一个碱基:[A]AAAATTTTTTTTTT
                    #修改指针后: AAAAA[T]TTTTTTTTT
                    #-----------
                #H对解析bam无影响，bam中不会显示H部分的序列
                if each_cigar[2]=='H': #直接跳过
                    pass
                #I会消耗reads的碱基，但ref的位置不变，所以只有reads的指针需要变化
                if each_cigar[2]=='I':#消耗reads，ref不变
                    read_map_status.append([query_pos,sequence[query_pos:query_pos+int(each_cigar[1])],ref_pos+1,"-",query_quality[query_pos-1:query_pos+int(each_cigar[1])]])
                    query_pos+=int(each_cigar[1])
                    #query_pos和query_pos+int(each[1])是ins在reads上的起始和终止位点
                    #ref_pos+1是插入位置, 应该是在ref_pos和ref_pos+1中间插入
                    #"-",表示ref型是'-'
                    #-----------------------
                #D会消耗ref的碱基，但reads的碱基位置不变，所以只有ref的指针需要变化
                if each_cigar[2]=='D':#消耗ref，reads不变
                    read_map_status.append([query_pos,"-",ref_pos+1,ref_sequence[ref_pos:ref_pos+int(each_cigar[1])],[40]])
                    ref_pos+=int(each_cigar[1])
                    
                #match时，ref和alt的指针都需要变化，这里先不管mismatch的情况
                if each_cigar[2]=='M':
                    for i in range(int(each_cigar[1])):
                        read_map_status.append([query_pos,sequence[query_pos],ref_pos+1,ref_sequence[ref_pos],query_quality[query_pos-1]])
                        ref_pos+=1
                        query_pos+=1
            #解析cigar结束
            #从解析完cigar的结果中提取所有indel
            
            for each_position in read_map_status:#[35, '-', 151856125, 'A']
#                if each_position[3]=="-": #insertion
#                elif each_position[1]=='-':#deletion
#                elif each_position[1]!=each_position[3]:#snv
#                    pass
                if each_position[3]=="-" or each_position[1]=='-':#ins or del
                    if pos == each_position[2] and refbase==each_position[3] and altbase==each_position[1]:
                        fragment_direction=parser_alignment_direction(eachread)
                        if not fragment_direction:
                            continue
                        tmp_base_len=len(each_position[-1])
                        tmp_low_quality_num=0
                        for each_quality in each_position[-1]:
                            if each_quality <= 10:
                                tmp_low_quality_num+=1
                        if each_position[3]=="-":#如果是indel就要考虑碱基质量值
                            if tmp_low_quality_num > float(tmp_base_len)/2:
                                continue
                        if simple_readid not in simple_readid_list:
                            simple_readid_list[simple_readid]=None
                        add_item_dict(eachread,chrname,fragment_direction)
                        assert readid not in readid_list
                        readid_list[readid]=None
                    
                else:
                    pass    
        snv_indel_support[eachsnv]=[list(readid_list.keys()),list(simple_readid_list.keys())]
    return snv_indel_support
                
def complex_snv(bamfile,snv_complex_dict,ref_dict):
    #处理复杂突变
    snv_complex_support_dict={}
    pattern=re.compile('((\d+)([SHMXNDI]))')#提取cigar列信息的pattern
    for eachsnv in list(snv_complex_dict.keys()):
        snv_complex_support_dict[eachsnv]=[[],[]]
        (chrname,pos,refbase,altbase)=snv_complex_dict[eachsnv]
        leftpos=pos-1#snv左边的位置
        rightpos=pos+len(refbase)#snv右边的位置
        might_support_reads=bamfile.get_pos_reads(chrname,leftpos,rightpos)
        cover_snv_reads={}
        simple_readid_readid={}
    

        for eachread in might_support_reads:            
            readid=eachread.query_name+'_'+str(eachread.flag)
            readid_simple=eachread.query_name
            query_quality=eachread.query_qualities
            fragment_direction=parser_alignment_direction(eachread)
            if not fragment_direction:
                continue
            cigar=eachread.cigarstring
            if not cigar:continue
            #解析cigar
            ref_pos=eachread.reference_start # 和ref_start初始值是一样的，但ref_pos后续用于ref的指针
            ref_start=eachread.reference_start
            sequence=eachread.query_sequence# 含softclip的序列
            m_chrname=bamfile.pysam_bam.get_reference_name(eachread.reference_id)
            try:
                assert chrname == m_chrname
            except:
                 raise Exception("some trouble with read "+readid_simple+" ,maybe there is secondary aligned for this read?")
            ref_sequence=ref_dict[chrname]            
            cigar_parse=re.findall(pattern,cigar)
            query_quality=eachread.query_qualities
            query_pos=0#作为指针，用于标记当前处理的碱基相对于reads起始碱基的位置
            read_map_status=[]
            for each_cigar in cigar_parse:#cigar中每个元素的处理, 每个子列表的第一个元素是不需要的[['',10,'S'],['',33,'M'],['',2,'I'],['',15,'M'],['',6,'D'],['',29,'M'],['',1,'D'],['',11,'M']]
                #S部分, reads序列会在bam中展示,起始位点其实是在reads中间的某个碱基
                if each_cigar[2]=='S':#跳过，不修改任何东西，但reads的坐标需要改动,ref坐标不变
                    query_pos+=int(each_cigar[1])
                    #----------
                    #5S10M
                    #指针位置在第一个碱基:[A]AAAATTTTTTTTTT
                    #修改指针后: AAAAA[T]TTTTTTTTT
                    #-----------
                #H对解析bam无影响，bam中不会显示H部分的序列
                if each_cigar[2]=='H': #直接跳过
                    pass
                #I会消耗reads的碱基，但ref的位置不变，所以只有reads的指针需要变化
                if each_cigar[2]=='I':#消耗reads，ref不变
                    read_map_status.append([query_pos,sequence[query_pos:query_pos+int(each_cigar[1])],ref_pos+1,"-",query_quality[query_pos:query_pos+int(each_cigar[1])]])
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
                        read_map_status.append([query_pos,sequence[query_pos],ref_pos+1,ref_sequence[ref_pos],query_quality[query_pos]])
                        ref_pos+=1
                        query_pos+=1
                 #解析cigar结束
            #从cigar判断read有没有cover snv
            #read_cover_flank_snv take place read_cover_single_snv, because of difference between single snv and complex snv
            readCover_left=read_cover_flank_snv(leftpos,read_map_status)
            readCover_right=read_cover_flank_snv(rightpos,read_map_status)
            if not readCover_left or not readCover_right:#必须同时cover两边的碱基
                continue
            #cover required refbase==readbase
            query_ref_pileup_left=bamfile.new_get_read_pos(leftpos,eachread,read_map_status)#[readpos,readbase,refpos,refbase],不止是没cover到，这个位置如果是个deletion，也会返回空值
            query_ref_pileup_right=bamfile.new_get_read_pos(rightpos,eachread,read_map_status)#[readpos,readbase,refpos,refbase])
            k=0#不cover
            if query_ref_pileup_left and query_ref_pileup_right:
                k=1#cover
            
            #需要refbase和altbase相同
            if query_ref_pileup_right and query_ref_pileup_left:
                query_seq=bamfile.get_read_seq(eachread,query_ref_pileup_left[0]+1,query_ref_pileup_right[0]-1)
                query_quality_snv_part=query_quality[query_ref_pileup_left[0]:query_ref_pileup_right[0]-1]
            elif query_ref_pileup_right and not query_ref_pileup_left:
                query_seq=bamfile.get_read_seq(eachread,1,query_ref_pileup_right[0]-1)
                query_quality_snv_part=query_quality[0:query_ref_pileup_right[0]-1]
            elif query_ref_pileup_left and not query_ref_pileup_right:
                query_seq=bamfile.get_read_seq(eachread,query_ref_pileup_left[0]+1,"end")
                query_quality_snv_part=query_quality[query_ref_pileup_left[0]:]
            else:
                continue
            if not query_seq:#只cover到了complex两边的位点，没有cover到complex位点
                continue
                
#            if eachread.query_name=='E100044124L1C014R00303299624':
#                bamfile.get_alignment(eachread)
            if query_seq==altbase:
                try:
                    low_freq=caculate_lowq_base_frequency(query_quality_snv_part)
                except:
                    print(eachsnv)
                    print(query_seq)
                    print(query_quality_snv_part)
                    print(query_ref_pileup_right)
                    print(query_ref_pileup_left)
                if low_freq>0.5:#低质量当没有cover处理
                    continue    
                if readid_simple in cover_snv_reads and readid_simple not in snv_complex_support_dict[eachsnv][1]:
                    continue
                snv_complex_support_dict[eachsnv][0].append(readid)
                if readid_simple not in snv_complex_support_dict[eachsnv][1]:
                    snv_complex_support_dict[eachsnv][1].append(readid_simple)
                add_item_dict(eachread,chrname,fragment_direction)
                simple_readid_readid[readid_simple]=readid
            else:#此reads不支持complex变异
                if k==1:#在此reads cover的前提下
                    if readid_simple in cover_snv_reads and readid_simple in snv_complex_support_dict[eachsnv][1]:#如果在此之前mate reads支持变异的话就删除这个fragment
                        del_index=snv_complex_support_dict[eachsnv][1].index(readid_simple)#删除的index
                        del snv_complex_support_dict[eachsnv][1][del_index]#删除这个fragment
                        del_index2=snv_complex_support_dict[eachsnv][0].index(simple_readid_readid[readid_simple])
                        del snv_complex_support_dict[eachsnv][0][del_index2]
            if k==1:#不受任何事情影响，只要算cover，那就加到字典里面去
                cover_snv_reads[readid_simple]=None
    return snv_complex_support_dict



def caculate_lowq_base_frequency(query_quality_f):
    tmp_base_len=len(query_quality_f)
    tmp_low_quality_num=0
    for each_quality in query_quality_f:
        if each_quality < 15:
            tmp_low_quality_num+=1
    return float(tmp_low_quality_num)/tmp_base_len
     

def snv_header_index(header_list,header):
    try:
        return header_list.index(header)
    except:
        raise Exception('Wrong header')
    

def parser_snv_file(infile):#读取snv文件，返回3个字典，2个列表
    snv_dict={}#{"chr1\t2491442\tC\tG":["chr1","2491442","C","G"],...}键和值实际是一样的，好像没什么必要用字典，但为了之后可能的修改，还是先和复杂变异的储存方式保持一致
    snv_indel_dict={}# T:TA 或者 -:A的形式
    snv_complex_dict={}#complex可能会带多个突变 {"chr1\t2491442\tCA\tGT":[["chr1","2491442","C","G"],["chr1","2491443","A","T"],...}
    raw_snv_dict={}# {"chr1\t2491442\tCA\tGT":line,...}, 将这四列作为一行的唯一标识，方便后续输出
    raw_snv_list=[]#["chr1\t2491442\tCA\tGT",...] 固定输出的顺序，保持和原snv结果一致
    snv_info_dict={}#{"chr1\t2491442\tCA\tGT":[info],...} 存储snv的信息，freq之类的
    with open(infile) as infile_handle:
        header=infile_handle.readline().strip('\n')
        header_info=header.split('\t')
        chrname_idx=snv_header_index(header_info,'chrom')
        pos_idx=snv_header_index(header_info,'pos_raw')
        refbase_idx=snv_header_index(header_info,'ref')
        altbase_idx=snv_header_index(header_info,'alt')
        type_idx=snv_header_index(header_info,'type')
        depth_idx=snv_header_index(header_info,'depth')
        support_reads_idx=snv_header_index(header_info,'support_reads')
        freq_idx=snv_header_index(header_info,'freq')
        try:
            black_flag_index=header_info.index('cutoff')
        except:
            black_flag_index=False
        while True:
            line=infile_handle.readline().strip('\n')
            if not line:break
            lineInfo=line.split('\t')
            chrname=lineInfo[chrname_idx]
            pos=lineInfo[pos_idx]
            refbase=lineInfo[refbase_idx]
            altbase=lineInfo[altbase_idx]
            freq=float(lineInfo[freq_idx])
            key='\t'.join([chrname,pos,refbase,altbase])
            raw_snv_list.append(key)
            raw_snv_dict[key]=line
            if black_flag_index:
                black_flag=lineInfo[black_flag_index]
            else:
                black_flag=""
            snv_info_dict[key]=[freq,black_flag]
            if freq>30:continue# will not to find support reads of snv freq >30
            #下面把单核苷酸突变,indel和complex种类突变分开，因为处理的方式不同
            if len(refbase)==1 and len(altbase)==1 and refbase!="-" and altbase!="-":#单核苷酸突变
                assert key not in snv_dict
                snv_dict[key]=[chrname,int(pos),refbase,altbase]
            elif (len(refbase)>1 or len(altbase)>1) and refbase[0]==altbase[0]:#indel的一种表现形式,比如 chr1 100 AT A, chr1 101 TTA T.
                pos=str(int(pos)+1)
                assert key not in snv_indel_dict
                refbase=refbase[1:]
                if not refbase:refbase='-'
                altbase=altbase[1:]
                if not altbase:altbase='-'
                snv_indel_dict[key]=[chrname,int(pos),refbase,altbase]
            elif refbase=="-" or altbase=="-":#indel用"-"表示的形式:
                assert key not in snv_indel_dict
                snv_indel_dict[key]=[chrname,int(pos),refbase,altbase]
            else:#复杂突变
                assert key not in snv_complex_dict
                snv_complex_dict[key]=[chrname,int(pos),refbase,altbase]
            snv_info_dict[key]=[freq,black_flag]
    return snv_dict,snv_indel_dict,snv_complex_dict,raw_snv_dict,raw_snv_list,header,snv_info_dict

def add_item_dict(read,chrname,fragment_direction):
    readid=read.query_name+'_'+str(read.flag)
    if readid not in READID_SEQUENCE_DICT:
        cigar=read.cigarstring
        ref_start=read.reference_start+1
        READID_SEQUENCE_DICT[readid]=[cigar,chrname,ref_start,read,fragment_direction]


if __name__=="__main__":
    argv=argument_parser()
    bamfile=argv['bamfile']
    fasta=argv['fasta']
    outdir=argv['outdir']
    outfile=argv['outfile']
    snvfile=argv['snv']


    out_prefix=os.path.join(outdir,outfile)

    #读取fasta
    ref_dict=read_fasta(fasta)

    #读bam文件
    bamfile=bamhunter.bamhunter(bamfile)
#    check_reads(bamfile,readid)
#    for eachread in bamfile.pysam_bam:
#        if eachread.query_name=='E100044124L1C014R00303299624':
#            bamfile.get_alignment(eachread)

    #解析snv文件
    (snv_dict,snv_indel_dict,snv_complex_dict,raw_snv_dict,raw_snv_list,snv_header,snv_info_dict)=parser_snv_file(snvfile)

    #单核苷酸突变支持reads检索
    snv_support_dict=single_base_snv(bamfile,snv_dict,ref_dict)
    READID_SEQUENCE_DICT={}

    snv_indel_support_dict=indel_snv(bamfile,snv_indel_dict,ref_dict)

    #complex snv
    snv_complex_support_dict=complex_snv(bamfile,snv_complex_dict,ref_dict)



    #最终输出
    final_outfile=out_prefix+'.snv_supportReads.txt'
    w=open(final_outfile,'w')
    w.write(snv_header+'\t'+'\t'.join(['support_reads_original','support_reads_num_nondedup','support_reads_num_dedup'])+'\n')
    for eachsnv in raw_snv_list:
        if eachsnv in list(snv_support_dict.keys()):
            supportreads=','.join(snv_support_dict[eachsnv][0])
            count_dup=str(len(snv_support_dict[eachsnv][0]))
            count_nondup=str(len(snv_support_dict[eachsnv][1]))
        elif eachsnv in list(snv_complex_support_dict.keys()):
            supportreads=','.join(snv_complex_support_dict[eachsnv][0])
            count_dup=str(len(snv_complex_support_dict[eachsnv][0]))
            count_nondup=str(len(snv_complex_support_dict[eachsnv][1]))
        elif eachsnv in list(snv_indel_support_dict.keys()):
            supportreads=','.join(snv_indel_support_dict[eachsnv][0])
            count_dup=str(len(snv_indel_support_dict[eachsnv][0]))
            count_nondup=str(len(snv_indel_support_dict[eachsnv][1]))
        else:
            supportreads="-"
            count_dup="-"
            count_nondup="-"
        w.write(raw_snv_dict[eachsnv]+'\t'+'\t'.join([supportreads,count_dup,count_nondup])+'\n')

    w.close()
    del snv_support_dict
    del snv_indel_support_dict
    del snv_complex_support_dict

