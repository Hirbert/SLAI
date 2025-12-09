#!coding:utf-8

#----------------------------------------------
#Project: 商检-基础流程-基础分析
#Description: 酶切过滤程序的核心模块-判断逻辑部分
#Usage: ./judge_artifact.py -h
#Author: 骆磊
#------------------------------------------------

import os
import argparse
import re
import support_reads as SR
#------------------
#read计数改为fragment计数 
#v0.4添加位置和softclip的限制
#v0.6允许softclip reads不到结尾允许有5bp的误差，位置关系放松至只需要有softclip就可以，不再要求方向。最小为9。hardclip直接认为是
#v0.7，修正频率直接过滤，过滤掉的位点不出现在最终文件中
def complement_reverse(seq):
    '''
    获取反向互补序列
    '''
    tmp=seq.upper()
    tmp=tmp.replace('A', 't')
    tmp=tmp.replace('T', 'a')
    tmp=tmp.replace('C', 'g')
    tmp=tmp.replace('G', 'c')
    tmp=tmp[::-1]
    return tmp.upper()
def argument_parser():
    parser = argparse.ArgumentParser(description="Judge one reads is artifact or not")
    parser.add_argument('--input',help='input file',required=True)
    parser.add_argument('--snv',help='input file',required=True)
    parser.add_argument('--outdir',help='input file',required=True)
    parser.add_argument('--outfile',help='input file',required=True)
    #/beegfs/work/commercial_test/cupcake/databases/gatk_bundle/2.8/hg19/ucsc.hg19.noconfig.fasta
    argv=vars(parser.parse_args())
    return argv

def artifact_pos(cigar,start,end,read,panelname):
    pattern=re.compile('((\d+)([SHMXNDI]))')#提取cigar列信息的pattern
    cigar_parse=re.findall(pattern,cigar)#3S94M [('6S', '3', 'S'), ('94M', '94', 'M')]
    pos_vars_nosoftclip=1
    pos_vars_softclip=5
    readlen=len(read)
    if cigar_parse[0][2]=='H' and not panelname=="panel43":
        return True,'5prime',1
    if cigar_parse[-1][2]=='H' and not panelname=="panel43":
        return True,'3prime',1
    if start <= 1+pos_vars_nosoftclip:#允许1个bp的误差
        pos='5prime'
        if cigar_parse[0][2]=='S' and int(cigar_parse[0][1])-(start-1)>=2:#softclip at least 2bp match 
            readtype=1# 有clip部分
        else:
            readtype=2# 无clip部分
        return True,pos,readtype
    elif end>=readlen-pos_vars_nosoftclip:#允许1个bp的误差
        pos='3prime'
        if cigar_parse[-1][2]=="S"  and int(cigar_parse[-1][1])-(readlen-end) >=2:#softclip at least 2bp match
            readtype=1# 有clip部分
        else:
            readtype=2# 无clip部分
        return True,pos,readtype
    elif 'S' in cigar:
        if cigar_parse[0][2]=='S':
            if start<=1+pos_vars_softclip and start<=1+int(cigar_parse[0][1])-2:#cover at least 2 base at softclip region
                pos='5prime'
                readtype=1# 有clip部分
                return True,pos,readtype
            elif start<=1+int(cigar_parse[0][1])-1:
                pos='5prime'
                readtype=2
                return True,pos,readtype
            else:
                return False,False,False
        elif cigar_parse[-1][2]=="S":
            if end>=readlen-pos_vars_softclip and end>=readlen-int(cigar_parse[-1][1]) +2:#cover at least 2 base at softclip region
                pos='3prime'
                readtype=1# 有clip部分
                return True,pos,readtype
            elif end>=readlen-int(cigar_parse[-1][1]) +1:
                pos='5prime'
                readtype=2
                return True,pos,readtype
            else:
               return False,False,False
    else:
        return False,False,False


    


def judge_read_artifact(readInfo,panelname):
    readid=readInfo[5]
    fragmentid='_'.join(readid.split('_')[:-1])#E100047969L1C012R00300322929_163, 这样处理应该没问题
    cigar=readInfo[7]
    matchscore=int(readInfo[8])
    ref_seq=readInfo[13]
    query_seq=readInfo[18]
    readstart=int(readInfo[15])
    readend=int(readInfo[16])
    read_seq=readInfo[-1]
    similar=float(readInfo[10])
    #根据比对位置过滤,鱼钩结构时，snv应该集中在read的一端，并且有softclip    
    snv_pos=int(readInfo[1])
    ref_base=readInfo[2]
    alt_base=readInfo[3]
    ref_start=int(readInfo[11])
    ref_end=int(readInfo[12])

    #互补配对区域
    ref_complement_region_start=int(readInfo[19])
    ref_complement_region_end=int(readInfo[20])
    passengers_field=readInfo[-2].split(',')
#    passengers_field=readInfo[-1].split(',')
    passengers=[]
    

    #indel longer than 3bp will bot be considered
    if len(ref_base)==1 and len(alt_base)>6 and ref_base[0]==alt_base[0]:#inset >5bp
        return False,False,False,fragmentid,"Non-AR1"
    if len(ref_base)>6 and len(alt_base)==1 and ref_base[0]==alt_base[0]:#deletion >5bp
        return False,False,False,fragmentid,"Non-AR1"



    for eachpassengers in passengers_field:
        passengers_pos=eachpassengers.split(':')[0]
        if passengers_pos==ref_complement_region_start or passengers_pos==ref_complement_region_end:
            pass
        else:
            if eachpassengers:
                passengers.append(eachpassengers)
    #根据比对分过滤
    if matchscore<=10:#比对分太低，这个read不该被认为是互补结合产生的
        return False,False,False,fragmentid,"Non-AR1"

    if snv_pos>ref_complement_region_end or snv_pos<ref_complement_region_start:
        return False,False,False,fragmentid,"Non-AR1"
    
    #根据比对长度过滤,应大于5
    if len(ref_seq) <=5 or len(query_seq)<=5:
        return False,False,False,fragmentid,"Non-AR1"

#    if len(ref_seq) <=9 or len(query_seq)<=9:
#        return False,False,False,fragmentid,"Non-AR2"
    #根据similar过滤，大于0.8
    if similar < 0.9:
        return False,False,False,fragmentid,"Non-AR1"
    (judge_f,pos_f,readtype_f)=artifact_pos(cigar,readstart,readend,read_seq,panelname)
    if not judge_f:
        if passengers:
            readlevel="confused"
        else:
            readlevel="Non-AR1"
    else:
        if readtype_f==1:
            readlevel="AR2"
        elif passengers and matchscore>27 :
            readlevel="AR1"
        elif len(query_seq)>=18:
            readlevel="Ambiguous-AR2"
        elif len(query_seq)>=10 and similar==1:
            readlevel="Ambiguous-AR1"
        elif len(query_seq)>=6:
            readlevel="Non-AR2"
        else:
            readlevel="unknown"
    return    judge_f,pos_f,readtype_f,fragmentid,readlevel
    
def snv_init_judge_dict(matchfile,snvfile):
    init_judge_dict={}#[total,true,false,type1,type2,5prime,3prime]
    match_dict={}
    with open(matchfile) as r:#没有header
         while True:
            line=r.readline().strip("\n")
            if not line:break
            lineInfo=line.split('\t')
            key='\t'.join(lineInfo[:4])
            try:
                match_dict[key].append(lineInfo)
            except KeyError:
                match_dict[key]=[lineInfo]
    with open(snvfile) as r:
        header=r.readline().strip("\n")#有header
        header_info=header.split('\t')
        while True:
            line=r.readline().strip("\n")
            if not line:break
            lineInfo=line.split('\t')
            key='\t'.join([lineInfo[SR.snv_header_index(header_info,'chrom')],lineInfo[SR.snv_header_index(header_info,'pos_raw')],lineInfo[SR.snv_header_index(header_info,'ref')],lineInfo[SR.snv_header_index(header_info,'alt')]])
            init_judge_dict[key]=[int(lineInfo[-1]),0,0,0,0,0,0,0,0,0,0,0,0,0]
    return init_judge_dict,match_dict

def snv_judge_dict(init_judge_dict_f,match_dict_f,match_reads_judge_file_f,panelname):
    level=['Non-AR1','unknown','confused','Non-AR2','Ambiguous-AR1','Ambiguous-AR2','AR1','AR2']
    w_matchfile=open(match_reads_judge_file_f,'w')
    for eachsnv in list(match_dict_f.keys()):
        tmp_fragment_id_dict={}#{fragmentid:""}用于存储已被计数的fragment
        for eachread in match_dict_f[eachsnv]:#这是matchresult的整行信息
            (judge,pos,readtype,fragmentid,readlevel)=judge_read_artifact(eachread,panelname)
            if fragmentid in tmp_fragment_id_dict:
                if level.index(readlevel)<level.index(tmp_fragment_id_dict[fragmentid][3]):
                    continue
            tmp_fragment_id_dict[fragmentid]=[judge,pos,readtype,readlevel,eachread[5]]
        used_fragment_id={}
        for eachread in match_dict_f[eachsnv]:#one more time
            fragmentid='_'.join(eachread[5].split('_')[:-1])#E100047969L1C012R00300322929_163, 这样处理应该没问题
            (judge,pos,readtype,readlevel,used_readid)=tmp_fragment_id_dict[fragmentid]
            if eachread[5] != used_readid:continue 
            if fragmentid in used_fragment_id:continue
            used_fragment_id[fragmentid]=""
            if judge and readlevel!="Non-AR2":
                yes_no="yes"
                init_judge_dict_f[eachsnv][1]+=1
                if readtype==1:
                    init_judge_dict_f[eachsnv][3]+=1
                if readtype==2:
                    init_judge_dict_f[eachsnv][4]+=1
                if pos=='5prime':
                    init_judge_dict_f[eachsnv][5]+=1
                if pos=='3prime':
                    init_judge_dict_f[eachsnv][6]+=1
    #read分级统计
                if readlevel=="Non-AR1":
                    init_judge_dict_f[eachsnv][8]+=1
                elif readlevel=="Non-AR2":
                    init_judge_dict_f[eachsnv][9]+=1
                elif readlevel=="Ambiguous-AR1":
                    init_judge_dict_f[eachsnv][10]+=1
                elif readlevel=="Ambiguous-AR2":
                    init_judge_dict_f[eachsnv][11]+=1
                elif readlevel=="AR1":
                    init_judge_dict_f[eachsnv][12]+=1
                elif readlevel=="AR2":
                    init_judge_dict_f[eachsnv][13]+=1
                else:
                    init_judge_dict_f[eachsnv][7]+=1
                w_matchfile.write('\t'.join(eachread)+'\t'+readlevel+'\tyes\n')
            else:
                init_judge_dict_f[eachsnv][2]+=1
                w_matchfile.write('\t'.join(eachread)+'\t'+readlevel+'\tno\n')
                if readlevel=="Non-AR1":
                    init_judge_dict_f[eachsnv][8]+=1
                elif readlevel=="Non-AR2":
                    init_judge_dict_f[eachsnv][9]+=1
        try:
            positive_reads_ratio=float(init_judge_dict_f[eachsnv][1])/init_judge_dict_f[eachsnv][0]
        except:#没有找support reads
            positive_reads_ratio=float(0)
            
        init_judge_dict_f[eachsnv].append(positive_reads_ratio)
        if positive_reads_ratio >=0.7:
            init_judge_dict_f[eachsnv].append('Positive')
        elif positive_reads_ratio <0.7 and positive_reads_ratio >=0.2:
            init_judge_dict_f[eachsnv].append('Ambiguous')
        else:
            init_judge_dict_f[eachsnv].append('Negative')

        #fix positive_reads_ratio and PN; Ambiguous-AR1  count 0.2
        ref=eachsnv.split('\t')[2]
        alt=eachsnv.split('\t')[3]
        if init_judge_dict_f[eachsnv][13]==0 and init_judge_dict_f[eachsnv][12]==0:
            coe=0
        else:
            coe=1
        if len(ref)==len(alt) and len(ref)>1 and complement_reverse(ref)==alt:
            coe=1


        fix_ar_reads=init_judge_dict_f[eachsnv][13] +init_judge_dict_f[eachsnv][12] +init_judge_dict_f[eachsnv][11] +init_judge_dict_f[eachsnv][10]*coe
        try:
            fix_positive_reads_ratio=float(fix_ar_reads)/init_judge_dict_f[eachsnv][0]
        except:
            fix_positive_reads_ratio=float(0)
        init_judge_dict_f[eachsnv].append(fix_positive_reads_ratio)
        if fix_positive_reads_ratio >=0.7:
            init_judge_dict_f[eachsnv].append('Positive')
        elif fix_positive_reads_ratio <0.7 and fix_positive_reads_ratio >=0.2:
            init_judge_dict_f[eachsnv].append('Ambiguous')
        else:
            init_judge_dict_f[eachsnv].append('Negative')


    w_matchfile.close()
    for eachsnv in list(init_judge_dict_f.keys()):
        if eachsnv not in match_dict_f:#没有reads有match
            init_judge_dict_f[eachsnv].extend(['0','Negative','0','Negative'])
        for i in range(len(init_judge_dict_f[eachsnv])):#all transfer to string format
            init_judge_dict_f[eachsnv][i]=str(init_judge_dict_f[eachsnv][i])

    return init_judge_dict_f

def write_output(output,judge_dict_f,snvfile_f):
    w=open(output,'w')
    with open(snvfile_f) as r:
        header=r.readline().strip("\n")
        header_info=header.split('\t')
        w.write(header+'\tReview_covered_reads\tArtifact_read_num\tNon_artifact_read_num\tArtifact_read_num_with_softclip\tArtifact_read_num_without_softclip\t5prime_read_num\t3prime_read_num\tOther_level\tNon-AR1\tNon-AR2\tAmbiguous-AR1\tAmbiguous-AR2\tAR1\tAR2\tArtifact_reads_ratio\tArtifact_snv_PN\tfix_Artifact_reads_ratio\tfix_Artifact_snv_PN\tfix_snv_freq\n')
        while True:
            line=r.readline().strip("\n")
            if not line:break
            lineInfo=line.split('\t')
            key='\t'.join([lineInfo[SR.snv_header_index(header_info,'chrom')],lineInfo[SR.snv_header_index(header_info,'pos_raw')],lineInfo[SR.snv_header_index(header_info,'ref')],lineInfo[SR.snv_header_index(header_info,'alt')]])
            assert key in judge_dict_f
            depth=int(lineInfo[SR.snv_header_index(header_info,'depth')])
            support_count=int(judge_dict_f[key][0])
            fix_Artifact_reads_ratio=float(judge_dict_f[key][-2])
            fix_snv_freq=(support_count*100*(1-fix_Artifact_reads_ratio))/depth
            w.write(line+'\t'+'\t'.join(judge_dict_f[key])+'\t'+str(fix_snv_freq)+'\n')
    w.close()

def write_final_output(output,judge_dict_f,snvfile_f):
    w=open(output,'w')
    with open(snvfile_f) as r:
        tmp_header=r.readline().strip("\n").split('\t')
        header='\t'.join(tmp_header[:-3])
        w.write(header+'\n')
        while True:
            line=r.readline().strip("\n")
            if not line:break
            lineInfo=line.split('\t')
            depth=int(lineInfo[SR.snv_header_index(tmp_header,'depth')])
            key='\t'.join([lineInfo[SR.snv_header_index(tmp_header,'chrom')],lineInfo[SR.snv_header_index(tmp_header,'pos_raw')],lineInfo[SR.snv_header_index(tmp_header,'ref')],lineInfo[SR.snv_header_index(tmp_header,'alt')]])
            assert key in judge_dict_f
            support_count=int(judge_dict_f[key][0])
            fix_Artifact_reads_ratio=float(judge_dict_f[key][-2])
            fix_snv_freq=(support_count*100*(1-fix_Artifact_reads_ratio))/depth
            if judge_dict_f[key][17]=="Positive":
                pn_tag="Y"
            elif judge_dict_f[key][17]=="Negative":
                pn_tag="N"
            elif judge_dict_f[key][17]=="Ambiguous":
                pn_tag="A"
            else:
                pn_tag="U"

            #---------下面的部分用于判断此snv是否应该被过滤------------------------

            #白名单header, 以防万一，还是将一二级位点和三级分开处理，设置不同的阈值
            try:
                exact_match=lineInfo[tmp_header.index('exact_match')]
            except:
                exact_match='no'
            try:
                range_match=lineInfo[tmp_header.index('range_match')]
            except:
                range_match='no'

            if exact_match!="no" or range_match!="no":
                freq_cutoff=1
            else:
                freq_cutoff=2

            #过滤规则，fliter_tag为True则过滤，False则不过滤
            if pn_tag=="N":
                fliter_tag=False
            elif pn_tag=="Y":
                fliter_tag=True
            elif pn_tag=="A":#这里目前还是先用矫正频率判断，不直接过滤掉
                if fix_snv_freq>=freq_cutoff:
                    fliter_tag=False
                else:
                    fliter_tag=True
            else:
                pass
            #-------------------------------------------------------------------

            #最终输出
            if not fliter_tag:
                w.write('\t'.join(lineInfo[:-3])+'\n')
    w.close()

if __name__=='__main__':
    my_argv=argument_parser()
    inputfile=my_argv['input']
    (init_judge_dict,match_dict)=snv_init_judge_dict(inputfile,my_argv['snv'])#[total,true,false,type1,type2,5prime,3prime]
    match_reads_judge_file=my_argv['outdir']+'/'+my_argv['outfile']+'.fix_artifact_match.txt'
    judge_dict=snv_judge_dict(init_judge_dict,match_dict,match_reads_judge_file)
    write_output(my_argv['outdir']+'/'+my_argv['outfile']+'.snv_artifact_judge.txt',judge_dict,my_argv['snv'])
    write_final_output(my_argv['outdir']+'/'+my_argv['outfile']+'.final.txt',judge_dict,my_argv['snv'])



