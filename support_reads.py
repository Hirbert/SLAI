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

#Read fasta file, return dictionary with chromosome as key and sequence as value. Reading too many tasks simultaneously may cause memory errors. To avoid this, a loop is added here. If fails, pause and retry. No loop termination set because normal reading should succeed.

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
        for eachread in might_support_reads:
            readid=eachread.query_name+'_'+str(eachread.flag)
            simple_readid=eachread.query_name
            cigar=eachread.cigarstring
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
                    if simple_readid not in simple_readid_list:
                        simple_readid_list[simple_readid]=None
                    key=1#support
                    simple_readid_readid[simple_readid]=readid
            if key==0:
                if  simple_readid in simple_readid_list:#Mate read was added before, current read does not support variant
                    del simple_readid_list[simple_readid]
                    del readid_list[simple_readid_readid[simple_readid]]

            cover_snv_readlist[simple_readid]=None
        snv_support_dict[eachsnv]=[list(readid_list.keys()),list(simple_readid_list.keys())]
    return snv_support_dict


def indel_snv(bamfile,snv_indel_dict,ref_dict):

    #Process indel
    snv_indel_support={}
    pattern=re.compile('((\d+)([SHMXNDI]))')#Pattern for extracting CIGAR information
    for eachsnv in list(snv_indel_dict.keys()):
        readid_list={}#support readid with flag
        simple_readid_list={}#support readid without flag
        (chrname,pos,refbase,altbase)=snv_indel_dict[eachsnv]
        might_support_reads=bamfile.get_pos_reads(chrname,pos-1,pos+len(refbase))#Ins and del are different, ins covers 3bp, but the difference is minimal, ignore
        for eachread in might_support_reads:
            readid=eachread.query_name+'_'+str(eachread.flag)
            simple_readid=eachread.query_name
            cigar=eachread.cigarstring
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
                        if each_position[3]=="-":#If indel, need to consider base quality
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
    #Process complex mutations
    snv_complex_support_dict={}
    pattern=re.compile('((\d+)([SHMXNDI]))')#Pattern for extracting CIGAR information
    for eachsnv in list(snv_complex_dict.keys()):
        snv_complex_support_dict[eachsnv]=[[],[]]
        (chrname,pos,refbase,altbase)=snv_complex_dict[eachsnv]
        leftpos=pos-1#Position left of SNV
        rightpos=pos+len(refbase)#Position right of SNV
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
            if not readCover_left or not readCover_right:#Must cover bases on both sides
                continue
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
            else:
                continue
            if not query_seq:#Only covered positions on both sides of complex, not the complex position itself
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

    #Read bam file
    bamfile=bamhunter.bamhunter(bamfile)
#    check_reads(bamfile,readid)
#    for eachread in bamfile.pysam_bam:
#        if eachread.query_name=='E100044124L1C014R00303299624':
#            bamfile.get_alignment(eachread)

    #Parse SNV file
    (snv_dict,snv_indel_dict,snv_complex_dict,raw_snv_dict,raw_snv_list,snv_header,snv_info_dict)=parser_snv_file(snvfile)

    #Single nucleotide mutation support reads retrieval
    snv_support_dict=single_base_snv(bamfile,snv_dict,ref_dict)
    READID_SEQUENCE_DICT={}

    snv_indel_support_dict=indel_snv(bamfile,snv_indel_dict,ref_dict)

    #complex snv
    snv_complex_support_dict=complex_snv(bamfile,snv_complex_dict,ref_dict)



    #Final output
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
