#!coding:utf-8

#----------------------------------------------
#Project: Commercial Inspection - Basic Process - Basic Analysis
#Description: Core module of enzyme digestion filtering program - alignment part, can also be used separately to output negative strand alignment results near SNV
#Usage: python3 artifact_identify.py -h
#Author: Luo Lei
#------------------------------------------------

#-----------------
#v0.4, support_reads modified from v2 to v3
#v0.5, support_reads_v5, read classification
import pysam
import argparse
import re
import os
import bamhunter
# from multi_thread import *  # The called script function has issues, will not be used
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
    Check whether matching causes passenger SNVs and whether the read supports passenger SNVs
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
        if each_cigar[2]=='I':#Consumes reads, ref unchanged
            cigar_status_1_2.append([query_pos,seq2[query_pos:query_pos+int(each_cigar[1])],ref_pos,'-'])#ref is 1-based
            query_pos+=int(each_cigar[1])
            #query_pos and query_pos+int(each[1]) are the start and end positions of insertion on reads
            #ref_pos+1 is the insertion position, should be inserted between ref_pos and ref_pos+1
            #"-", indicates ref type is '-'
            #-----------------------
        #D consumes reference bases, read base position unchanged, so only ref pointer needs to change
        if each_cigar[2]=='D':#Consumes ref, reads unchanged
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
        if each_cigar[2]=='I':#Consumes reads, ref unchanged
            cigar_status_1_3.append([query_pos,seq3[query_pos:query_pos+int(each_cigar[1])],ref_pos,'-'])#ref is 1-based
            query_pos+=int(each_cigar[1])
            #query_pos and query_pos+int(each[1]) are the start and end positions of insertion on reads
            #ref_pos+1 is the insertion position, should be inserted between ref_pos and ref_pos+1
            #"-", indicates ref type is '-'
            #-----------------------
        #D consumes reference bases, read base position unchanged, so only ref pointer needs to change
        if each_cigar[2]=='D':#Consumes ref, reads unchanged
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
    Map the matched region of the read to the genome
    '''
    (chrname,pos,ref_base,alt_base,m_score,m_cigar,new_start1,new_end1,mseq1,reverse_mseq1,start2,end2,mseq2,reverse_mseq2,similarity)=match_results_f
    m_cigar=reverse_cigar(m_cigar)
    start2-=1
    end2-=1
    ref_pos-=1
    new_start1-=1
    new_end1-=1
    query_pos=0#Used as a pointer to mark the current base position relative to the read start base
    pattern=re.compile('((\d+)([SHMXNDI]))')#Pattern for extracting CIGAR information
    cigar_parse=re.findall(pattern,cigar)
    read_map_status=[]
    read_complement_pos=[]#Store complement positions of the read

    for each_cigar in cigar_parse:#Process each element in CIGAR, the first element of each sublist is not needed [['',10,'S'],['',33,'M'],['',2,'I'],['',15,'M'],['',6,'D'],['',29,'M'],['',1,'D'],['',11,'M']]

        #S part, read sequence will be displayed in bam, the start position is actually somewhere in the middle of the read
        if each_cigar[2]=='S':#Skip, don't modify anything, but read coordinates need to change, ref coordinates unchanged
            query_pos+=int(each_cigar[1])
        #H has no effect on bam parsing, bam does not display the H part sequence
        if each_cigar[2]=='H': #Skip directly
            pass
        #I consumes read bases, ref position unchanged, so only read pointer needs to change
        if each_cigar[2]=='I':#Consumes reads, ref unchanged
            read_map_status.append([query_pos,sequence[query_pos:query_pos+int(each_cigar[1])],ref_pos+1,"-"])
            query_pos+=int(each_cigar[1])
            #query_pos and query_pos+int(each[1]) are the start and end positions of insertion on reads
            #ref_pos+1 is the insertion position, should be inserted between ref_pos and ref_pos+1
            #"-", indicates ref type is '-'
            #-----------------------
        #D consumes reference bases, read base position unchanged, so only ref pointer needs to change
        if each_cigar[2]=='D':#Consumes ref, reads unchanged
            read_map_status.append([query_pos,"-",ref_pos+1,ref_sequence[ref_pos:ref_pos+int(each_cigar[1])]])
            ref_pos+=int(each_cigar[1])

        #For match, both ref and alt pointers need to change, ignore mismatch case for now
        if each_cigar[2]=='M':
            for i in range(int(each_cigar[1])):
                read_map_status.append([query_pos,sequence[query_pos],ref_pos+1,ref_sequence[ref_pos]])
                ref_pos+=1
                query_pos+=1
    k=0
    for each_base_map in read_map_status:
        #I consumes read bases, ref position unchanged, so only read pointer needs to change
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
    ref_seq_c=reverse_mseq1# Reverse complement of the genomic sequence of the matched portion
    ref_start_c=new_end1# Due to reverse complement, the start position actually comes from the end position
    ref_end_c=new_start1
    query_pos_c=start2

    m_cigar_parse=re.findall(pattern,m_cigar)
    for each_m_cigar in m_cigar_parse:
        if each_m_cigar[2]=='I':#Consumes reads, ref unchanged
            read_complement_cigar_status.append([query_pos_c,ref_start_c+1])#ref is 1-based
            query_pos_c+=int(each_m_cigar[1])
            #query_pos and query_pos+int(each[1]) are the start and end positions of insertion on reads
            #ref_pos+1 is the insertion position, should be inserted between ref_pos and ref_pos+1
            #"-", indicates ref type is '-'
            #-----------------------
        #D consumes reference bases, read base position unchanged, so only ref pointer needs to change
        if each_m_cigar[2]=='D':#Consumes ref, reads unchanged
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
    This function is used to remove SNVs that are not considered palindrome positive. Sites with frequency >30 or blacklisted are not processed
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
    Read SNV file with support reads information and store in dictionary
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
    Before alignment, replace ref sequence with alt
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
    Infer genomic position from read base position
    '''
    for eachpos in read_ref_pair_pos_dict[readid]:
        if eachpos[0]==pos-1:
            return eachpos[1]+1
    raise Exception()

def get_position_ref_read(read_ref_pair_pos_dict,readid,pos):
    '''
    Infer read position from reference position
    '''
    aligned_pairs=read_ref_pair_pos_dict[readid][3].get_aligned_pairs()
    for eachpos in aligned_pairs:
        if eachpos[1]==pos-1:
        #chr15    90630344    C    G    15    E100044124L1C030R01000347883_163    90630333    4S94M2S    40    20M    1.0    90630334    90630353    CCCAGCGTACCCTGGGCCAG    CTGGCCCAGGGTACGCTGGG    1    20    CTGGCCCAGGGTACGCTGGG    CCCAGCGTACCCTGGGCCAG    CTGGCCCAGGGTACGCTGGGCCAGGATGTCTGACTGCACATCTCCGTCATAGTTCTTGCAGGCCCACACAAAGCCACCCGAAGACTTGAGGACCTGAGAT
#[(0, 90631630), (1, 90631631), (2, 90631632), (3, 90631633), (4, 90631634), (5, 90631635), ******(None, 90631636)******, (6, 90631637), (7, 90631638), (8, 90631639), (9, 90631640), (10, 90631641), (11, None), (12, 90631642), (13, 90631643), (14, 90631644), (15, 90631645), (16, 90631646), (17, 90631647), (18, 90631648), (19, 90631649), (20, 90631650), (21, 90631651), (22, 90631652), (23, 90631653), (24, 90631654), (25, 90631655), (26, 90631656), (27, 90631657),(28, 90631658), (29, 90631659), (30, 90631660), (31, 90631661), (32, 90631662), (33, 90631663), (34, 90631664), (35, 90631665), (36, 90631666), (37, 90631667), (38, 90631668), (39, 90631669), (40, 90631670), (41, 90631671), (42, 90631672), (43, 90631673), (44, 90631674), (45, 90631675), (46, 90631676), (47, 90631677), (48, 90631678), (49, 90631679), (50, 90631680), (51, 90631681), (52, 90631682), (53, 90631683), (54, 90631684), (55, 90631685), (56, 90631686), (57, 90631687), (58, 90631688), (59, 90631689), (60, 90631690), (61, 90631691), (62, 90631692), (63, 90631693), (64, 90631694), (65, 90631695), (66, 90631696), (67, 90631697), (68, 90631698), (69, 90631699), (70, 90631700), (71, 90631701), (72, 90631702), (73, 90631703), (74, 90631704), (75, 90631705), (76, 90631706), (77, 90631707), (78, 90631708), (79, 90631709), (80, 90631710), (81, 90631711), (82, 90631712), (83, 90631713), (84, 90631714), (85, 90631715), (86, 90631716), (87, 90631717), (88, 90631718)] Deletion occurring at 90631636
            if eachpos[0] or eachpos[0] == 0:#If deletion, read pos with deletion is none
                return eachpos[0]+1
            else:
                return last_pos+1
        last_pos=eachpos[0]
        #Insertion at read edge, return 1
    if aligned_pairs[0][1]==None and aligned_pairs[0][0]==0:
        return 1
    print(aligned_pairs)
    print(readid)
    print(pos)
    raise Exception()
def get_refseq(chrname,pos,ref,alt,ref_dict,flank):
    '''
    Get reference genome sequence near the SNV
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
    #-------------------Initial preparation
    #Read fasta
#    print("Read fasta")
    ref_dict=read_fasta(fasta)
#    print("Read fasta over")
    out_prefix=os.path.join(outdir,outfile)


    #Process bam to get support reads for each snv, this step can be skipped to be compatible with support reads obtained by different methods
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
        if filter_snv(eachsnv,snv_info_dict):#Pre-determine whether it is a palindrome false positive based on SNV information, if not, skip the subsequent analysis and directly judge as negative
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
                complement_raw_region_modify=change_sequence_by_snv(complement_raw_start,complement_raw_end,complement_raw_region,eachsnv)#The variant must exist within complement_raw_region range, if not, an error will occur here, this situation should not exist
                passengers=find_passengers(complement_raw_region_modify,complement_reverse(complement_mate_region),complement_read_region)#[[8, 'C', 8, 'T'], [11, 'A', 11, 'G'], [20, 'T', 20, 'G']]
                really_passengers=[]
                for eachpassenger in passengers:
                    fix_pos=eachpassenger[0]+complement_raw_start
                    really_passengers.append(str(fix_pos)+':'+eachpassenger[3]+':'+eachpassenger[1])
                w.write('\t'.join([chrname,str(pos),ref_base,alt_base,str(snv_pos_in_read),read,str(read_start),cigar,str(m_score),str(new_m_cigar),str(similarity),str(new_start1),str(new_end1),mseq1,reverse_mseq1,str(start2),str(end2),mseq2,reverse_mseq2,str(complement_raw_start),str(complement_raw_end),complement_raw_region,str(complement_mate_start),str(complement_mate_end),complement_mate_region,complement_reverse(complement_mate_region),str(complement_read_start),str(complement_read_end),complement_read_region,sr.READID_SEQUENCE_DICT[read][-1],','.join(really_passengers),sr.READID_SEQUENCE_DICT[read][3].query_sequence])+'\n')
#                w.write('\t'.join([chrname,str(pos),ref_base,alt_base,str(snv_pos_in_read),read,str(read_start),cigar,str(m_score),str(m_cigar),str(similarity),str(new_start1),str(new_end1),mseq1,reverse_mseq1,str(start2),str(end2),mseq2,reverse_mseq2,READID_SEQUENCE_DICT[read][3]])+'\n')

if __name__=='__main__':
    main()
