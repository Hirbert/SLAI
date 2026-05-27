#!coding:utf-8
#------------------------------------------------
#v3: Added filtering rules, align with commercial inspection standards
#v4.3: rm single snv read, lowq, confulict
#------------------------------------
#Project: Commercial Inspection - Basic Process - Basic Analysis
#Description: Based on SNV file and BAM file, backtrack support reads for each SNV, output the IDs of support reads for each SNV. support_reads_original: support read IDs; support_reads_num_nondedup: number of support reads obtained from backtracking, read1 and read2 counted separately when both support; support_reads_num_dedup: number of support reads obtained from backtracking, read1 and read2 counted as one when both support
#Usage: python3 ./support_reads.py -h
#Author: Luo Lei
import pysam
import argparse
import re
import os
import bamhunter
import time
import random


#Global variable
global READID_SEQUENCE_DICT,READ_CIGAR_DICT
READID_SEQUENCE_DICT={}
READ_CIGAR_DICT={}

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

#Read fasta file, return dictionary with chromosome as key and sequence as value. Reading too many tasks simultaneously may cause memory errors. To avoid this, a loop is added here. If fails, pause and retry. No loop termination set because normal reading should succeed.

def pos_softclip_ratio(ref_dict,bamfile,chrname,pos1_original,pos2_original,exclude_list=[]):
    pattern=re.compile('((\d+)([SHMXNDI]))')#Pattern for extracting CIGAR information
    depth=0
    softnum=0
    pos1=pos1_original-50
    pos2=pos2_original+50
    for eachread in bamfile.get_pos_reads(chrname,pos1,pos2):
        depth+=1
        cigar=eachread.cigarstring
        if not cigar:continue
        if not "S" in cigar:continue
        cigar_parse=re.findall(pattern,cigar)
        sequence=eachread.query_sequence# Sequence containing softclip
        readid=eachread.query_name+'_'+str(eachread.flag)
        ref_sequence=ref_dict[chrname]
        read_map_status=[]
        ref_pos=eachread.reference_start # Same as initial ref_start, but ref_pos is subsequently used as ref pointer
        ref_start=eachread.reference_start
        query_quality=eachread.query_qualities
        query_pos=0#Used as pointer to mark current base position relative to read start
        for each_cigar in cigar_parse:#Process each element in CIGAR, first element of each sublist is not needed [['',10,'S'],['',33,'M'],['',2,'I'],['',15,'M'],['',6,'D'],['',29,'M'],['',1,'D'],['',11,'M']]
        #S part, read sequence will be displayed in bam, start position is actually somewhere in the middle of the read
            if each_cigar[2]=='S':#Skip, don't modify anything, but read coordinates need to change, ref coordinates unchanged
                query_pos+=int(each_cigar[1])
                #----------
                #5S10M
                #Pointer at first base: [A]AAAATTTTTTTTTT
                #After modifying pointer: AAAAA[T]TTTTTTTTT
                #-----------
            #H has no effect on bam parsing, bam does not display H part sequence
            if each_cigar[2]=='H': #Skip directly
                pass
            #I consumes read bases, ref position unchanged, so only read pointer needs to change
            if each_cigar[2]=='I':#Consumes reads, ref unchanged
                read_map_status.append([query_pos,sequence[query_pos:query_pos+int(each_cigar[1])],ref_pos+1,"-",query_quality[query_pos-1:query_pos+int(each_cigar[1])]])
                query_pos+=int(each_cigar[1])
                #query_pos and query_pos+int(each[1]) are start and end positions of insertion on reads
                #ref_pos+1 is insertion position, should be inserted between ref_pos and ref_pos+1
                #"-", indicates ref type is '-'
                #-----------------------
            #D consumes reference bases, read base position unchanged, so only ref pointer needs to change
            if each_cigar[2]=='D':#Consumes ref, reads unchanged
                read_map_status.append([query_pos,"-",ref_pos+1,ref_sequence[ref_pos:ref_pos+int(each_cigar[1])]])
                ref_pos+=int(each_cigar[1])
                
            #For match, both ref and alt pointers need to change, ignore mismatch case for now
            if each_cigar[2]=='M':
                for i in range(int(each_cigar[1])):
                    read_map_status.append([query_pos,sequence[query_pos],ref_pos+1,ref_sequence[ref_pos],query_quality[query_pos]])
                    ref_pos+=1
                    query_pos+=1
        #CIGAR parsing end
        if cigar_parse[0][2]=="S" and float(cigar_parse[0][1])>=1:#softclip upstream
            s_len=int(cigar_parse[0][1])
            r_start=read_map_status[0][2]
#            if r_start-s_len<=pos2 and r_start>=pos1 and not readid in exclude_list:
#            print(readid)
#            print(pos1_original)
#            print(pos2_original)
#            print(r_start)
#            print(s_len)
            if r_start-s_len<=pos2_original and r_start>=pos1_original and not readid in exclude_list:
                softnum+=1
#                print(readid+'\t5\t'+str(pos1_original)+"\t"+str(pos2_original)+"\t"+str(r_start)+'\t'+str(s_len))
        if cigar_parse[-1][2]=="S" and float(cigar_parse[-1][1])>=1:
            s_len=int(cigar_parse[-1][1])
            r_end=read_map_status[-1][2]
            if r_end<=pos2_original and r_end+s_len>=pos1_original and not readid in exclude_list:
                softnum+=1
#                print(readid+'\t3'+str(pos1_original)+"\t"+str(pos2_original)+"\t"+str(r_end)+'\t'+str(s_len))
#    print(softnum)
    
    ratio=round(float(softnum)/depth,6)
    return softnum
def parser_alignment_direction(readhandle):
    if readhandle.is_secondary or readhandle.is_supplementary or readhandle.is_unmapped or readhandle.is_duplicate:
        return False
    chr1=readhandle.reference_id    
    chr2=readhandle.next_reference_id
    if chr1!=chr2:
        samechr=False
    else:
        samechr=True
    #When read is read1:
    if readhandle.is_read1:
        if not samechr:#Mapped to different chromosomes, in this case only single-strand direction can be considered
            if readhandle.is_reverse:#R1
                return "R1_diffchr"
            else:
                return "F1_diffchr"
        if readhandle.mate_is_unmapped:#mate unmapped case
            if readhandle.is_reverse:#R1
                return "R1_unmap"
            else:
                return "F1_unmap"
        if True:#properly + unusual insert size
            if readhandle.is_reverse: # R1 case
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
    else:#When read is read2
        if not samechr:#Mapped to different chromosomes, in this case only single-strand direction can be considered
            if readhandle.is_reverse:#R2
                return "R2_diffchr"
            else:
                return "F2_diffchr"
        if readhandle.mate_is_unmapped:#mate unmapped case
            if readhandle.is_reverse:#R2
                return "R2_unmap"
            else:
                return "F2_unmap"
#        if readhandle.is_proper_pair or readhandle.template_length <1000:# Normal alignment case
        if True:#properly + unusual insert size
            if readhandle.is_reverse: # R2 case
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

def read_cover_single_pos(pos,read_map_ref):
    for eachitem in read_map_ref:
        if eachitem[2]==pos:
            try:
                return 1
            except TypeError:
                pass
    return 0

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
    
    pattern=re.compile('((\d+)([SHMXNDI]))')#Pattern for extracting CIGAR information
    for eachsnv in list(snv_dict.keys()):
        (chrname,pos,refbase,altbase)=snv_dict[eachsnv]
        ratio=pos_softclip_ratio(ref_dict,bamfile,chrname,pos,pos)
        might_support_reads=bamfile.get_pos_reads(chrname,pos-1,pos+1)#Get all reads covering this position, but limited by pysam, minimum can only find reads covering 2bp, so these reads may not actually be aligned. Also, if this position is a deletion, it will also be judged as not aligned. No impact on results, temporarily not modified. If the queried region is a mismatch, it will also be judged as not aligned, so need to find reads covering a range, tentatively set to 50, since the possibility of 50 consecutive mismatches is very low.
        #---------------------------------------------
        #Why find reads covering pos-1 and pos+1?
        #First reason: pysam can only find reads covering a "region", minimum 2bp
        #Second reason: mismatched positions are not considered covered. For example, if a read has a mismatch at position 100, this read will not be extracted. So if only looking for reads covering a single SNV position, all reads with mismatches will be lost.
        #Therefore, need to find reads covering positions near the SNV, with OR relationship. But if neither the preceding nor following positions are match, insertion, softclip, or read end, they cannot be found either. So we assume that the bp before and after the SNV must be match. This should be the case, otherwise the SNV would be detected as another variant rather than a single nucleotide mutation. The same processing method applies to indel and complex mutations below.
        #---------------------------------------------
        #What we actually need to find are reads where both read and reference have bases at this position, but the bases are different
        #If it's a multi-base mutation, reads supporting the multi-base mutation will still support the single-base mutation
        readid_list={}
        simple_readid_list={}#support readid without flag
        cover_snv_readlist={}#List of fragments covering this position
        simple_readid_readid={}#{simple_readid:readid}
        cigar_readid_list={}
        for eachread in might_support_reads:
            readid=eachread.query_name+'_'+str(eachread.flag)
            simple_readid=eachread.query_name
            cigar=eachread.cigarstring
            SAtag=eachread.has_tag("SA")
            READ_CIGAR_DICT[readid]=[cigar,SAtag]
            if not cigar:continue
            #Parse CIGAR
            ref_pos=eachread.reference_start # Same as initial ref_start, but ref_pos is subsequently used as ref pointer
            ref_start=eachread.reference_start
            sequence=eachread.query_sequence# Sequence containing softclip
            m_chrname=bamfile.pysam_bam.get_reference_name(eachread.reference_id)
            try:
                assert chrname == m_chrname
            except:
                print("some trouble with read "+simple_readid+" ,maybe there is secondary aligned for this read?")
            ref_sequence=ref_dict[chrname]
            cigar_parse=re.findall(pattern,cigar)
            query_quality=eachread.query_qualities
            query_pos=0#Used as pointer to mark current base position relative to read start
            read_map_status=[]
            for each_cigar in cigar_parse:#Process each element in CIGAR, first element of each sublist is not needed [['',10,'S'],['',33,'M'],['',2,'I'],['',15,'M'],['',6,'D'],['',29,'M'],['',1,'D'],['',11,'M']]
            #S part, read sequence will be displayed in bam, start position is actually somewhere in the middle of the read
                if each_cigar[2]=='S':#Skip, don't modify anything, but read coordinates need to change, ref coordinates unchanged
                    query_pos+=int(each_cigar[1])
                    #----------
                    #5S10M
                    #Pointer at first base: [A]AAAATTTTTTTTTT
                    #After modifying pointer: AAAAA[T]TTTTTTTTT
                    #-----------
                #H has no effect on bam parsing, bam does not display H part sequence
                if each_cigar[2]=='H': #Skip directly
                    pass
                #I consumes read bases, ref position unchanged, so only read pointer needs to change
                if each_cigar[2]=='I':#Consumes reads, ref unchanged
                    read_map_status.append([query_pos,sequence[query_pos:query_pos+int(each_cigar[1])],ref_pos+1,"-",query_quality[query_pos-1:query_pos+int(each_cigar[1])]])
                    query_pos+=int(each_cigar[1])
                    #query_pos and query_pos+int(each[1]) are start and end positions of insertion on reads
                    #ref_pos+1 is insertion position, should be inserted between ref_pos and ref_pos+1
                    #"-", indicates ref type is '-'
                    #-----------------------
                #D consumes reference bases, read base position unchanged, so only ref pointer needs to change
                if each_cigar[2]=='D':#Consumes ref, reads unchanged
                    read_map_status.append([query_pos,"-",ref_pos+1,ref_sequence[ref_pos:ref_pos+int(each_cigar[1])]])
                    ref_pos+=int(each_cigar[1])
                    
                #For match, both ref and alt pointers need to change, ignore mismatch case for now
                if each_cigar[2]=='M':
                    for i in range(int(each_cigar[1])):
                        read_map_status.append([query_pos,sequence[query_pos],ref_pos+1,ref_sequence[ref_pos],query_quality[query_pos]])
                        ref_pos+=1
                        query_pos+=1
            #CIGAR parsing end
            #Determine whether read covers SNV from CIGAR
            readCover=read_cover_single_snv(pos,read_map_status)
            if not readCover:continue#Read does not cover SNV position or SNV position quality too low
#            if readid=="E100047980L1C004R00203310736_163":print(read_map_status)
#            readCover_left=read_cover_flank_snv(pos-1,read_map_status)
#            readCover_right=read_cover_flank_snv(pos+1,read_map_status)
#            if not readCover_left and readCover_right:
#                continue
            fragment_direction=parser_alignment_direction(eachread)
            if not fragment_direction:
                continue
            #Extract all from parsed CIGAR results
            key=0# not support
            for each_position in read_map_status:#[35, '-', 151856125, 'A'], 0 is 0-based, 2 is 1-based
#                if each_position[3]=="-": #insertion
#                elif each_position[1]=='-':#deletion
#                elif each_position[1]!=each_position[3]:#sn 
                
                if pos == each_position[2] and refbase==each_position[3] and altbase==each_position[1]:# single base snv support reads
                    readCover_left=read_cover_LR_single_snv(pos-1,read_map_status)
                    readCover_right=read_cover_LR_single_snv(pos+1,read_map_status)
                    if not readCover_left and readCover_right:
                        continue
                    add_item_dict(eachread,chrname,fragment_direction)
                    if simple_readid in cover_snv_readlist and simple_readid not in simple_readid_list:#Mate covers but does not support variant
                        break
                    assert readid not in readid_list
                    readid_list[readid]=None
                    cigar_readid_list[readid]=None
                    if simple_readid not in simple_readid_list:
                        simple_readid_list[simple_readid]=None
                    key=1#support
                    simple_readid_readid[simple_readid]=readid
            if key==0:
                if  simple_readid in simple_readid_list:#Mate read was added before, current read does not support variant
                    del simple_readid_list[simple_readid]
                    del readid_list[simple_readid_readid[simple_readid]]
                    del cigar_readid_list[simple_readid_readid[simple_readid]]

            cover_snv_readlist[simple_readid]=None
        snv_support_dict[eachsnv]=[list(readid_list.keys()),list(simple_readid_list.keys()),ratio,list(cigar_readid_list.keys())]
    return snv_support_dict


def indel_snv(bamfile,snv_indel_dict,ref_dict):

    #Process indel
    snv_indel_support={}
    pattern=re.compile('((\d+)([SHMXNDI]))')#Pattern for extracting CIGAR information
    
    for eachsnv in list(snv_indel_dict.keys()):
        rearrangement_read=[]
        readid_list={}#support readid with flag 
        simple_readid_list={}#support readid without flag
        cigar_readid_list={}
        (chrname,pos,refbase,altbase)=snv_indel_dict[eachsnv]
        if altbase=="-" or len(refbase) > len(altbase):
            mutation_type="Deletion"
        elif refbase=="-" or len(refbase) < len(altbase):
            mutation_type="Insertion"
        else:
            raise Exception
        might_support_reads=bamfile.get_pos_reads(chrname,pos-1,pos+len(refbase)+1)#Ins and del are different, ins covers 3bp, but the difference is minimal, ignore
        
        for eachread in might_support_reads:
            readid=eachread.query_name+'_'+str(eachread.flag)
            simple_readid=eachread.query_name
            cigar=eachread.cigarstring
            SAtag=eachread.has_tag("SA")
            READ_CIGAR_DICT[readid]=[cigar,SAtag]
            if not cigar:continue
            if 'H' in cigar:continue
            #Parse CIGAR
            ref_pos=eachread.reference_start # Same as initial ref_start, but ref_pos is subsequently used as ref pointer
            ref_start=eachread.reference_start
            sequence=eachread.query_sequence# Sequence containing softclip
            query_quality=eachread.query_qualities
            assert len(query_quality)==len(sequence)
            m_chrname=bamfile.pysam_bam.get_reference_name(eachread.reference_id)
            try:
                assert chrname == m_chrname
            except:
                print("some trouble with read "+simple_readid+" ,maybe there is secondary aligned for this read?")
            ref_sequence=ref_dict[chrname]
            cigar_parse=re.findall(pattern,cigar)
            query_pos=0#Used as pointer to mark current base position relative to read start
            read_map_status=[]
            for each_cigar in cigar_parse:#Process each element in CIGAR, first element of each sublist is not needed [['',10,'S'],['',33,'M'],['',2,'I'],['',15,'M'],['',6,'D'],['',29,'M'],['',1,'D'],['',11,'M']]
            #S part, read sequence will be displayed in bam, start position is actually somewhere in the middle of the read
                if each_cigar[2]=='S':#Skip, don't modify anything, but read coordinates need to change, ref coordinates unchanged
                    query_pos+=int(each_cigar[1])
                    #----------
                    #5S10M
                    #Pointer at first base: [A]AAAATTTTTTTTTT
                    #After modifying pointer: AAAAA[T]TTTTTTTTT
                    #-----------
                #H has no effect on bam parsing, bam does not display H part sequence
                if each_cigar[2]=='H': #Skip directly
                    pass
                #I consumes read bases, ref position unchanged, so only read pointer needs to change
                if each_cigar[2]=='I':#Consumes reads, ref unchanged
                    read_map_status.append([query_pos,sequence[query_pos:query_pos+int(each_cigar[1])],ref_pos+1,"-",query_quality[query_pos-1:query_pos+int(each_cigar[1])]])
                    query_pos+=int(each_cigar[1])
                    #query_pos and query_pos+int(each[1]) are start and end positions of insertion on reads
                    #ref_pos+1 is insertion position, should be inserted between ref_pos and ref_pos+1
                    #"-", indicates ref type is '-'
                    #-----------------------
                #D consumes reference bases, read base position unchanged, so only ref pointer needs to change
                if each_cigar[2]=='D':#Consumes ref, reads unchanged
                    read_map_status.append([query_pos,"-",ref_pos+1,ref_sequence[ref_pos:ref_pos+int(each_cigar[1])],[40]])
                    ref_pos+=int(each_cigar[1])
                    
                #For match, both ref and alt pointers need to change, ignore mismatch case for now
                if each_cigar[2]=='M':
                    for i in range(int(each_cigar[1])):
                        read_map_status.append([query_pos,sequence[query_pos],ref_pos+1,ref_sequence[ref_pos],query_quality[query_pos-1]])
                        ref_pos+=1
                        query_pos+=1
            #CIGAR parsing end
            #Extract all indels from parsed CIGAR results
            #14S87M
            #[[14, 'A', 55242480, 'A', 37], [15, 'A', 55242481, 'A', 37], [16, 'C', 55242482, 'C', 37], [17, 'A', 55242483, 'A', 37], [18, 'T', 55242484, 'T', 37], [19, 'C', 55242485, 'C', 37], [20, 'T', 55242486, 'T', 37], [21, 'C', 55242487, 'C', 37], [22, 'C', 55242488, 'C', 37], [23, 'G', 55242489, 'G', 37], [24, 'A', 55242490, 'A', 37], [25, 'A', 55242491, 'A', 37], [26, 'A', 55242492, 'A', 37], [27, 'G', 55242493, 'G', 37], [28, 'C', 55242494, 'C', 37], [29, 'C', 55242495, 'C', 37], [30, 'A', 55242496, 'A', 37], [31, 'A', 55242497, 'A', 37], [32, 'C', 55242498, 'C', 37], [33, 'A', 55242499, 'A', 37], [34, 'A', 55242500, 'A', 37], [35, 'G', 55242501, 'G', 37], [36, 'G', 55242502, 'G', 37], [37, 'A', 55242503, 'A', 37], [38, 'A', 55242504, 'A', 37], [39, 'A', 55242505, 'A', 37], [40, 'T', 55242506, 'T', 37], [41, 'C', 55242507, 'C', 37], [42, 'C', 55242508, 'C', 37], [43, 'T', 55242509, 'T', 37], [44, 'C', 55242510, 'C', 37], [45, 'G', 55242511, 'G', 37], [46, 'A', 55242512, 'A', 37], [47, 'T', 55242513, 'T', 37], [48, 'G', 55242514, 'G', 37], [49, 'T', 55242515, 'T', 37], [50, 'G', 55242516, 'G', 37], [51, 'A', 55242517, 'A', 37], [52, 'G', 55242518, 'G', 37], [53, 'T', 55242519, 'T', 37], [54, 'T', 55242520, 'T', 37], [55, 'T', 55242521, 'T', 37], [56, 'C', 55242522, 'C', 37], [57, 'T', 55242523, 'T', 37], [58, 'G', 55242524, 'G', 37], [59, 'C', 55242525, 'C', 37], [60, 'T', 55242526, 'T', 37], [61, 'T', 55242527, 'T', 25], [62, 'T', 55242528, 'T', 37], [63, 'G', 55242529, 'G', 37], [64, 'C', 55242530, 'C', 37], [65, 'T', 55242531, 'T', 37], [66, 'G', 55242532, 'G', 37], [67, 'T', 55242533, 'T', 37], [68, 'G',55242534, 'G', 37], [69, 'T', 55242535, 'T', 37], [70, 'G', 55242536, 'G', 37], [71, 'G', 55242537, 'G', 37], [72, 'G', 55242538, 'G', 37], [73, 'G', 55242539, 'G', 37], [74, 'G', 55242540, 'G', 37], [75, 'T', 55242541, 'T', 37], [76, 'C', 55242542, 'C', 37], [77, 'C', 55242543, 'C', 37], [78, 'A', 55242544, 'A', 37], [79, 'T', 55242545, 'T', 37], [80, 'G', 55242546, 'G', 37], [81, 'G', 55242547, 'G', 37], [82, 'C', 55242548, 'C', 37], [83, 'T', 55242549, 'T', 37], [84, 'C', 55242550, 'C', 37], [85, 'T', 55242551, 'T', 37], [86, 'G', 55242552, 'G', 37], [87, 'A', 55242553, 'A', 37], [88, 'A', 55242554, 'A', 37], [89, 'C', 55242555, 'C', 37], [90, 'C', 55242556, 'C', 37], [91, 'T', 55242557, 'T', 37], [92, 'C', 55242558, 'C', 37], [93, 'A', 55242559, 'A', 37], [94, 'G', 55242560, 'G', 37], [95, 'G', 55242561, 'G', 37], [96, 'C', 55242562, 'C', 37], [97, 'C', 55242563, 'C', 37], [98, 'C', 55242564, 'C', 37], [99, 'A', 55242565, 'A', 37], [100, 'C', 55242566, 'C', 37]]

            support_tag=False
            for each_position in read_map_status:#[35, '-', 151856125, 'A']
                if each_position[3]=="-" or each_position[1]=='-':#ins or del
                    if pos == each_position[2] and refbase==each_position[3] and altbase==each_position[1]:
                        support_tag=True
                        fragment_direction=parser_alignment_direction(eachread)
                        if not fragment_direction:
                            continue
                        if simple_readid not in simple_readid_list:
                            simple_readid_list[simple_readid]=None
                        add_item_dict(eachread,chrname,fragment_direction)
                        assert readid not in readid_list
                        readid_list[readid]=None
                        cigar_readid_list[readid]=None
                else:
                    pass 
            #rearrangement
            #Reconstruct reference sequence
            if mutation_type=="Insertion":
                readCover1=read_cover_single_pos(pos-1,read_map_status)
                readCover2=read_cover_single_pos(pos,read_map_status)
                
                new_ref_sequence_upstream=ref_sequence[pos-100:pos]+altbase[1:]
                new_ref_sequence_downstream=altbase+ref_sequence[pos-1:pos+99]

            elif mutation_type=="Deletion":
                readCover1=read_cover_single_pos(pos-2,read_map_status)
                readCover2=read_cover_single_pos(pos+len(refbase)+1,read_map_status)
                new_ref_sequence_upstream=ref_sequence[pos-101:pos-1]
                new_ref_sequence_downstream=ref_sequence[pos+len(refbase)-1:pos+len(refbase)+99]
            else:
                raise Exception
            if not support_tag:
#                if cigar_parse[0][2] in "IS" and readCover1+readCover2==1:#Upstream end is insert or softclip
                if cigar_parse[0][2] in "IS" and readCover1+readCover2==1:#Upstream end is insert or softclip
                    read_part_seq=sequence[:int(cigar_parse[0][1])]
                    # print('up')
                    # print(readid)
                    # print(read_part_seq)
                    # print(new_ref_sequence_upstream)
                    if new_ref_sequence_upstream.endswith(read_part_seq):
                        support_tag=True
                if cigar_parse[-1][2] in "IS" and readCover1+readCover2==1:#Downstream end is insert or softclip
                    read_part_seq=sequence[-int(cigar_parse[-1][1]):]
                    # print('down')
                    # print(readid)
                    # print(read_part_seq)
                    # print(new_ref_sequence_downstream)
                    if new_ref_sequence_downstream.startswith(read_part_seq):
                        support_tag=True
                if support_tag:
                    rearrangement_read.append(readid)
                    # print(readid+'\t'+cigar+'\trearrangement')
                    fragment_direction=parser_alignment_direction(eachread)
                    if not fragment_direction:
                        continue
                    if simple_readid not in simple_readid_list:
                        simple_readid_list[simple_readid]=None
                    add_item_dict(eachread,chrname,fragment_direction)
                    assert readid not in readid_list
                    readid_list[readid]=None

        ratio=pos_softclip_ratio(ref_dict,bamfile,chrname,pos-1,pos+len(refbase)+1,rearrangement_read)
        snv_indel_support[eachsnv]=[list(readid_list.keys()),list(simple_readid_list.keys()),ratio,rearrangement_read,list(cigar_readid_list.keys())]
    # print("step1 done")
    return snv_indel_support
                
def complex_snv(bamfile,snv_complex_dict,ref_dict):
    #Process complex mutations
    snv_complex_support_dict={}
    pattern=re.compile('((\d+)([SHMXNDI]))')#Pattern for extracting CIGAR information
    for eachsnv in list(snv_complex_dict.keys()):
        rearrangement_read=[]
        snv_complex_support_dict[eachsnv]=[[],[]]
        (chrname,pos,refbase,altbase)=snv_complex_dict[eachsnv]
        leftpos=pos-1#Position left of SNV
        rightpos=pos+len(refbase)#Position right of SNV
        
        might_support_reads=bamfile.get_pos_reads(chrname,leftpos-1,rightpos+1)
        cover_snv_reads={}
        simple_readid_readid={}
        cigar_support_list={}
    

        for eachread in might_support_reads:            
            readid=eachread.query_name+'_'+str(eachread.flag)
            readid_simple=eachread.query_name
            query_quality=eachread.query_qualities
            fragment_direction=parser_alignment_direction(eachread)
            if not fragment_direction:
                continue
            cigar=eachread.cigarstring
            SAtag=eachread.has_tag("SA")
            READ_CIGAR_DICT[readid]=[cigar,SAtag]
            if not cigar:continue
            #Parse CIGAR
            ref_pos=eachread.reference_start # Same as initial ref_start, but ref_pos is subsequently used as ref pointer
            ref_start=eachread.reference_start
            sequence=eachread.query_sequence# Sequence containing softclip
            m_chrname=bamfile.pysam_bam.get_reference_name(eachread.reference_id)
            try:
                assert chrname == m_chrname
            except:
                 raise Exception("some trouble with read "+readid_simple+" ,maybe there is secondary aligned for this read?")
            ref_sequence=ref_dict[chrname]            
            cigar_parse=re.findall(pattern,cigar)
            query_quality=eachread.query_qualities
            query_pos=0#Used as pointer to mark current base position relative to read start
            read_map_status=[]
            for each_cigar in cigar_parse:#Process each element in CIGAR, first element of each sublist is not needed [['',10,'S'],['',33,'M'],['',2,'I'],['',15,'M'],['',6,'D'],['',29,'M'],['',1,'D'],['',11,'M']]
                #S part, read sequence will be displayed in bam, start position is actually somewhere in the middle of the read
                if each_cigar[2]=='S':#Skip, don't modify anything, but read coordinates need to change, ref coordinates unchanged
                    query_pos+=int(each_cigar[1])
                    #----------
                    #5S10M
                    #Pointer at first base: [A]AAAATTTTTTTTTT
                    #After modifying pointer: AAAAA[T]TTTTTTTTT
                    #-----------
                #H has no effect on bam parsing, bam does not display H part sequence
                if each_cigar[2]=='H': #Skip directly
                    pass
                #I consumes read bases, ref position unchanged, so only read pointer needs to change
                if each_cigar[2]=='I':#Consumes reads, ref unchanged
                    read_map_status.append([query_pos,sequence[query_pos:query_pos+int(each_cigar[1])],ref_pos+1,"-",query_quality[query_pos:query_pos+int(each_cigar[1])]])
                    query_pos+=int(each_cigar[1])
                    #query_pos and query_pos+int(each[1]) are start and end positions of insertion on reads
                    #ref_pos+1 is insertion position, should be inserted between ref_pos and ref_pos+1
                    #"-", indicates ref type is '-'
                    #-----------------------
                #D consumes reference bases, read base position unchanged, so only ref pointer needs to change
                if each_cigar[2]=='D':#Consumes ref, reads unchanged
                    read_map_status.append([query_pos,"-",ref_pos+1,ref_sequence[ref_pos:ref_pos+int(each_cigar[1])]])
                    ref_pos+=int(each_cigar[1])
                    
                #For match, both ref and alt pointers need to change, ignore mismatch case for now
                if each_cigar[2]=='M':
                    for i in range(int(each_cigar[1])):
                        read_map_status.append([query_pos,sequence[query_pos],ref_pos+1,ref_sequence[ref_pos],query_quality[query_pos]])
                        ref_pos+=1
                        query_pos+=1
                 #CIGAR parsing end
            
            #Determine whether read covers SNV from CIGAR
            #read_cover_flank_snv takes the place of read_cover_single_snv, because of difference between single snv and complex snv
            readCover_left=read_cover_flank_snv(leftpos,read_map_status)
            readCover_right=read_cover_flank_snv(rightpos,read_map_status)
#            if not readCover_left or not readCover_right:#Must cover bases on both sides
#                continue
            
            #cover required refbase==readbase
            query_ref_pileup_left=bamfile.new_get_read_pos(leftpos,eachread,read_map_status)#[readpos,readbase,refpos,refbase], if not covered or if this position is a deletion, returns empty
            query_ref_pileup_right=bamfile.new_get_read_pos(rightpos,eachread,read_map_status)#[readpos,readbase,refpos,refbase])
            k=0#not cover
            if query_ref_pileup_left and query_ref_pileup_right:
                k=1#cover
            
            #Need refbase and altbase to be the same
            if query_ref_pileup_right and query_ref_pileup_left:
                query_seq=bamfile.get_read_seq(eachread,query_ref_pileup_left[0]+1,query_ref_pileup_right[0]-1)
                query_quality_snv_part=query_quality[query_ref_pileup_left[0]:query_ref_pileup_right[0]-1]
            elif query_ref_pileup_right and not query_ref_pileup_left:
                query_seq=bamfile.get_read_seq(eachread,1,query_ref_pileup_right[0]-1)
                query_quality_snv_part=query_quality[0:query_ref_pileup_right[0]-1]
            elif query_ref_pileup_left and not query_ref_pileup_right:
                query_seq=bamfile.get_read_seq(eachread,query_ref_pileup_left[0]+1,"end")
                query_quality_snv_part=query_quality[query_ref_pileup_left[0]:]
            elif "S" in cigar:
                query_seq=""
                pass
            else:
                continue

            if not query_seq and not "S" in cigar:#Only covered positions on both sides of complex, not the complex position itself
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
                if low_freq>0.5:#Treat low quality as not covered
                    continue    
                if readid_simple in cover_snv_reads and readid_simple not in snv_complex_support_dict[eachsnv][1]:
                    continue
                snv_complex_support_dict[eachsnv][0].append(readid)
                cigar_support_list[readid]=None
                if readid_simple not in snv_complex_support_dict[eachsnv][1]:
                    snv_complex_support_dict[eachsnv][1].append(readid_simple)
                add_item_dict(eachread,chrname,fragment_direction)
                simple_readid_readid[readid_simple]=readid
            else:#This read does not support complex variant
                if k==1:#Under the premise that this read covers
                    if readid_simple in cover_snv_reads and readid_simple in snv_complex_support_dict[eachsnv][1]:#If previously mate reads supported variant, delete this fragment
                        del_index=snv_complex_support_dict[eachsnv][1].index(readid_simple)#Index to delete
                        del snv_complex_support_dict[eachsnv][1][del_index]#Delete this fragment
                        del_index2=snv_complex_support_dict[eachsnv][0].index(simple_readid_readid[readid_simple])
                        del snv_complex_support_dict[eachsnv][0][del_index2]
            if k==1:#Regardless of anything, as long as it counts as cover, add to dictionary
                cover_snv_reads[readid_simple]=None
            
            #rearrangement
            readCover1=read_cover_single_pos(leftpos-1,read_map_status)
            readCover2=read_cover_single_pos(rightpos+1,read_map_status)

            if readCover1+readCover2!=1 or not "S" in cigar:
                continue
            
            new_ref_sequence_upstream=ref_sequence[leftpos-100:leftpos]+altbase
            new_ref_sequence_downstream=altbase+ref_sequence[rightpos-1:rightpos+99]
            support_tag=False
            if cigar_parse[0][2] in "IS" and readCover1+readCover2==1:#Upstream end is insert or softclip
                read_part_seq=sequence[:int(cigar_parse[0][1])]
                # print('up')
                # print(readid)
                # print(read_part_seq)
                # print(new_ref_sequence_upstream)
                if new_ref_sequence_upstream.endswith(read_part_seq):
                    support_tag=True
            if cigar_parse[-1][2] in "IS" and readCover1+readCover2==1:#Downstream end is insert or softclip
                read_part_seq=sequence[-int(cigar_parse[-1][1]):]
                # print('down')
                # print(readid)
                # print(read_part_seq)
                # print(new_ref_sequence_downstream)
                if new_ref_sequence_downstream.startswith(read_part_seq):
                    support_tag=True
            
            if support_tag:
                # print(readid+'\trearrangement')
                snv_complex_support_dict[eachsnv][0].append(readid)
                if readid_simple not in snv_complex_support_dict[eachsnv][1]:
                    snv_complex_support_dict[eachsnv][1].append(readid_simple)
                add_item_dict(eachread,chrname,fragment_direction)
                simple_readid_readid[readid_simple]=readid
                rearrangement_read.append(readid)
        cigar_support_num=len(cigar_support_list.keys())
        ratio=pos_softclip_ratio(ref_dict,bamfile,chrname,leftpos,leftpos,rearrangement_read)     
        snv_complex_support_dict[eachsnv].extend([ratio,rearrangement_read,cigar_support_num])

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
    

def read_vcf(infile):
    vcf_dict={}
    with open(infile,encoding='utf-8') as infile_handle:
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
            black_flag_index=header_info.index('black_list_flag')
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
            if black_flag_index:
                black_flag=lineInfo[black_flag_index]
            else:
                black_flag=""
            vcf_dict[key]={"chrname":chrname,"pos":pos,"refbase":refbase,"altbase":altbase,
            "freq":freq,"black_flag":black_flag}
    return vcf_dict

def split_dict(in_dict):#Read SNV file, return 3 dictionaries, 2 lists
    snv_dict={}#{"chr1\t2491442\tC\tG":["chr1","2491442","C","G"],...} Key and value are actually the same, dictionary may not seem necessary, but for possible future modifications, keep consistent with complex variant storage method
    snv_indel_dict={}# Form like T:TA or -:A
    snv_complex_dict={}#Complex may contain multiple mutations {"chr1\t2491442\tCA\tGT":[["chr1","2491442","C","G"],["chr1","2491443","A","T"],...}
    raw_snv_dict={}# {"chr1\t2491442\tCA\tGT":line,...}, Use these four columns as unique identifier for the row, convenient for subsequent output
    raw_snv_list=[]#["chr1\t2491442\tCA\tGT",...] Fixed output order, consistent with original SNV results
    snv_info_dict={}#{"chr1\t2491442\tCA\tGT":[info],...} Store SNV information, freq etc.
    for key in in_dict:
        value=in_dict[key]
        chrname=value["chrname"]
        pos=value["pos"]
        refbase=value["refbase"]
        altbase=value["altbase"]
        freq=value["freq"]
        black_flag=value["black_flag"]
        raw_snv_list.append(key)
        raw_snv_dict[key]=key
        if len(refbase)==1 and len(altbase)==1 and refbase!="-" and altbase!="-":#Single nucleotide mutation
            assert key not in snv_dict
            snv_dict[key]=[chrname,int(pos),refbase,altbase]
        elif (len(refbase)>1 or len(altbase)>1) and refbase[0]==altbase[0]:#A form of indel, e.g., chr1 100 AT A, chr1 101 TTA T.
            pos=str(int(pos)+1)
            assert key not in snv_indel_dict
            refbase=refbase[1:]
            if not refbase:refbase='-'
            altbase=altbase[1:]
            if not altbase:altbase='-'
            snv_indel_dict[key]=[chrname,int(pos),refbase,altbase]
        elif refbase=="-" or altbase=="-":#Indel represented with "-":
            assert key not in snv_indel_dict
            snv_indel_dict[key]=[chrname,int(pos),refbase,altbase]
        else:#Complex mutation
            assert key not in snv_complex_dict
            snv_complex_dict[key]=[chrname,int(pos),refbase,altbase]
        snv_info_dict[key]=[freq,black_flag]
    return snv_dict,snv_indel_dict,snv_complex_dict,raw_snv_dict,raw_snv_list,snv_info_dict

def parser_snv_file(infile):#Read SNV file, return 3 dictionaries, 2 lists
    snv_dict={}#{"chr1\t2491442\tC\tG":["chr1","2491442","C","G"],...} Key and value are actually the same, dictionary may not seem necessary, but for possible future modifications, keep consistent with complex variant storage method
    snv_indel_dict={}# Form like T:TA or -:A
    snv_complex_dict={}#Complex may contain multiple mutations {"chr1\t2491442\tCA\tGT":[["chr1","2491442","C","G"],["chr1","2491443","A","T"],...}
    raw_snv_dict={}# {"chr1\t2491442\tCA\tGT":line,...}, Use these four columns as unique identifier for the row, convenient for subsequent output
    raw_snv_list=[]#["chr1\t2491442\tCA\tGT",...] Fixed output order, consistent with original SNV results
    snv_info_dict={}#{"chr1\t2491442\tCA\tGT":[info],...} Store SNV information, freq etc.
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
            black_flag_index=header_info.index('black_list_flag')
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
#            if freq>30:continue# will not to find support reads of snv freq >30
            #Separate single nucleotide mutations, indels, and complex mutations below because they are processed differently
            if len(refbase)==1 and len(altbase)==1 and refbase!="-" and altbase!="-":#Single nucleotide mutation
                assert key not in snv_dict
                snv_dict[key]=[chrname,int(pos),refbase,altbase]
            elif (len(refbase)>1 or len(altbase)>1) and refbase[0]==altbase[0]:#A form of indel, e.g., chr1 100 AT A, chr1 101 TTA T.
                pos=str(int(pos)+1)
                assert key not in snv_indel_dict
                refbase=refbase[1:]
                if not refbase:refbase='-'
                altbase=altbase[1:]
                if not altbase:altbase='-'
                snv_indel_dict[key]=[chrname,int(pos),refbase,altbase]
            elif refbase=="-" or altbase=="-":#Indel represented with "-":
                assert key not in snv_indel_dict
                snv_indel_dict[key]=[chrname,int(pos),refbase,altbase]
            else:#Complex mutation
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

def support_softclip_ratio(supportreads,exclude_list=[]):
    pattern=re.compile('((\d+)([SHMXNDI]))')#Pattern for extracting CIGAR information
    total=0
    leftnum=0
    rightnum=0
    bothnum=0
    SAtagnum=0
    for read in supportreads.split(','):
        if not read or read=='-':
            break
        total+=1
        cigar=READ_CIGAR_DICT[read][0]
        satag=READ_CIGAR_DICT[read][1]
        if satag:SAtagnum+=1
        cigar_parse=re.findall(pattern,cigar)
        if cigar_parse[0][2]=="S" and int(cigar_parse[0][1])>=5 and read not in exclude_list:
            leftnum+=1
        if cigar_parse[-1][2]=="S" and int(cigar_parse[-1][1])>=5 and read not in exclude_list:
            rightnum+=1
        if cigar_parse[0][2]=="S" and cigar_parse[-1][2]=="S" and int(cigar_parse[0][1])>=5 and int(cigar_parse[-1][1])>=5 and  read not in exclude_list:
            bothnum+=1
    if total!=0:
        left_ratio=round(float(leftnum)/total,6)
        right_ratio=round(float(rightnum)/total,6)
        both_ratio=round(float(bothnum)/total,6)
    else:
        left_ratio="NA"
        right_ratio="NA"
        both_ratio="NA"
    return str(left_ratio),str(right_ratio),str(both_ratio),str(total),str(SAtagnum)

def find_support_reads(raw_snv_list,snv_support_dict,snv_indel_support_dict,snv_complex_support_dict):
    for eachsnv in raw_snv_list:
        if eachsnv in list(snv_support_dict.keys()):
            supportreads=','.join(snv_support_dict[eachsnv][0])
            count_dup=str(len(snv_support_dict[eachsnv][0]))
            count_nondup=str(len(snv_support_dict[eachsnv][1]))
            softnum=int(snv_support_dict[eachsnv][2])
            rearrangement_read_list=[]
            cigar_support_num=str(len(snv_support_dict[eachsnv][3]))
        elif eachsnv in list(snv_complex_support_dict.keys()):
            supportreads=','.join(snv_complex_support_dict[eachsnv][0])
            count_dup=str(len(snv_complex_support_dict[eachsnv][0]))
            count_nondup=str(len(snv_complex_support_dict[eachsnv][1]))
            softnum=int(snv_complex_support_dict[eachsnv][2])
            rearrangement_read_list=snv_complex_support_dict[eachsnv][3]
            cigar_support_num=str(snv_complex_support_dict[eachsnv][4])
        elif eachsnv in list(snv_indel_support_dict.keys()):
            supportreads=','.join(snv_indel_support_dict[eachsnv][0])
            count_dup=str(len(snv_indel_support_dict[eachsnv][0]))
            count_nondup=str(len(snv_indel_support_dict[eachsnv][1]))
            softnum=int(snv_indel_support_dict[eachsnv][2])
            rearrangement_read_list=snv_indel_support_dict[eachsnv][3]
            cigar_support_num=str(len(snv_indel_support_dict[eachsnv][4]))
        else:
            supportreads="-"
            count_dup="-"
            count_nondup="-"
            softnum='-'
            cigar_support_num="-"
            rearrangement_read_list=[]
        left_ratio,right_ratio,both_ratio,traceback_support_reads_num,SA_support_reads_num=support_softclip_ratio(supportreads,rearrangement_read_list)
        ratio = str(round(float(softnum)/(softnum+int(traceback_support_reads_num)),6)) if softnum!='-' else softnum

        outline=raw_snv_dict[eachsnv]+'\t'+'\t'.join([supportreads,count_dup,count_nondup,left_ratio,right_ratio,both_ratio,ratio])
#        print('\t'.join([eachsnv,left_ratio,right_ratio,both_ratio,ratio,cigar_support_num,count_dup,supportreads]))
    return left_ratio,right_ratio,both_ratio,ratio,traceback_support_reads_num,SA_support_reads_num,outline,cigar_support_num


if __name__=="__main__":
    argv=argument_parser()
    bamfile=argv['bamfile']
    fasta=argv['fasta']
    outdir=argv['outdir']
    outfile=argv['outfile']
    snvfile=argv['snv']


    out_prefix=os.path.join(outdir,outfile)

    #Read fasta
    ref_dict=read_fasta(fasta)
    print("read fasta done")

    #Read bam file
    bamfile=bamhunter.bamhunter(bamfile)
#    check_reads(bamfile,readid)
#    for eachread in bamfile.pysam_bam:
#        if eachread.query_name=='E100044124L1C014R00303299624':
#            bamfile.get_alignment(eachread)

    #Parse SNV file
    vcf_dict=read_vcf(snvfile)


    #Read header
    r=open(snvfile)
    snv_header=r.readline().strip()
    r.close()

    (snv_dict,snv_indel_dict,snv_complex_dict,raw_snv_dict,raw_snv_list,snv_info_dict)=split_dict(vcf_dict)
#    (snv_dict,snv_indel_dict,snv_complex_dict,raw_snv_dict,raw_snv_list,snv_header,snv_info_dict)=parser_snv_file(snvfile)
    print("read vcf done")
    #Single nucleotide mutation support reads retrieval
    snv_support_dict=single_base_snv(bamfile,snv_dict,ref_dict)
    READID_SEQUENCE_DICT={}

    snv_indel_support_dict=indel_snv(bamfile,snv_indel_dict,ref_dict)

    #complex snv
    snv_complex_support_dict=complex_snv(bamfile,snv_complex_dict,ref_dict)
    print("parse snv done")


    #Final output
    final_outfile=out_prefix+'.snv_supportReads.txt'
    w=open(final_outfile,'w')
    w.write(snv_header+'\t'+'\t'.join(['support_reads_original','support_reads_num_nondedup','support_reads_num_dedup','left_softclip_ratio','right_softclip_ratio','both_softclip_ratio','neighbor_softclip_ratio'])+'\n')
    left_ratio,right_ratio,both_ratio,ratio,outline,a,b,c=find_support_reads(raw_snv_list,snv_support_dict,snv_indel_support_dict,snv_complex_support_dict)
    
    w.write(outline+'\n')

    w.close()
