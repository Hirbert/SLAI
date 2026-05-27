#!coding:utf-8
#----------------------------------------------
#Project: Commercial Inspection - Basic Process - Basic Analysis
#Description: Core module of enzyme digestion filtering program - alignment part, can also be used separately to output negative strand alignment results near SNV
#Usage: python3 artifact_identify.py -h
#Author: Luo Lei
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
    Parse external parameters
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
    Input sequence, return the corresponding reverse complementary sequence
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
    Return dictionary {chr1:'AGCATGCTAGCTACGATC...',....}
    '''
    ref_dict={}#{chr:sequence,...}
    fasta_handle=pysam.FastxFile(fasta)
    for each_chr in fasta_handle:
        ref_dict[each_chr.name]=each_chr.sequence.upper()
    return ref_dict

def get_fasta_sequence(ref_dict,chrname,begin,end):
    '''
    Input genomic coordinates, return sequence based on coordinates
    '''
    return ref_dict[chrname][begin-1:end]

def structure_identify(site,ref_dict):
    '''
    Extract the reference genome sequence near the SNV, denote it as A, obtain the reverse complementary sequence of A as B, and align A and B using dynamic programming
    '''
    sites=site.split(':')
    chrname=sites[0]
    pos=int(sites[1])
    ref_base=sites[2]
    alt_base=sites[3]
    refseq=get_fasta_sequence(ref_dict,chrname,pos-flank,pos+flank-1+len(ref_base))#Reference sequence including the variant and its upstream/downstream
    query_1=get_fasta_sequence(ref_dict,chrname,pos-flank,pos-1)#Reference sequence upstream of the variant
    query_2=get_fasta_sequence(ref_dict,chrname,pos+len(ref_base),pos+flank-1+len(ref_base))#Reference sequence downstream of the variant
    query=complement_reverse(query_1+alt_base+query_2)#Reverse complementary sequence of the variant sequence, used for local alignment with reference sequence


    #Local alignment
    (m_score,m_cigar,start1,end1,mseq1,start2,end2,mseq2,x,y,z,similarity)=SmithWaterman.local_alignment(refseq,query,0,flank+1)#Reverse complement of the sequence covering SNV aligns to reference genome sequence
    #-------------------------------------------------------------------------------------------
    #(m_score,m_cigar,start1,end1,mseq1,start2,end2,mseq2,x,y,z)=SmithWaterman.local_alignment(refseq,query,0,0)
    #m_score Alignment score
    #m_cigar Alignment details
    #start1 Alignment start site of extracted forward reference genome sequence
    #end1 Alignment end site of extracted forward reference genome sequence
    #mseq1 Aligned forward reference genome sequence
    #start2 Start site of the reverse complementary sequence of the variant sequence (forward end site)
    #end2 End site of the reverse complementary sequence of the variant sequence (forward start site)
    #mseq2 Reverse complementary sequence of the variant sequence
    #-----------------------------------------------------------------------------


    #Analysis after alignment
    if not int(m_score):#Site without palindrome false positive basis
        return None

    #Map alignment results to reference genome coordinates
    new_start1= pos - flank - 1 +int(start1)#Actual alignment start site on reference genome sequence
    new_end1=pos - flank - 1 + int(end1)#Actual alignment end site on reference genome sequence
    start2_r=2*flank+len(alt_base)-int(end2)+1#Start site of the reverse complementary sequence of the variant sequence in forward orientation
    end2_r=2*flank+len(ref_base)-int(start2)+1#End site of the reverse complementary sequence of the variant sequence in forward orientation
    new_start2=pos-flank-1 +start2_r#Start site of the reverse complementary sequence of the variant sequence on the forward genome strand
    new_end2=pos-flank-1+end2_r#End site of the reverse complementary sequence of the variant sequence on the forward genome strand
    reverse_query_seq=complement_reverse(mseq2)#Reverse complement of the aligned region of the variant sequence (forward orientation)
    return [chrname,pos,ref_base,alt_base,m_score,m_cigar,new_start1,new_end1,new_start2,mseq1,new_end2,reverse_query_seq,mseq2,similarity]

def read_ref_alignment(read,chrname,pos,ref_base,alt_base,ref_dict,flank,keysite):
    refseq=get_fasta_sequence(ref_dict,chrname,pos-flank,pos+flank-1+len(ref_base))#Reference sequence including the variant and its upstream/downstream
    query=complement_reverse(read)#Reverse complementary sequence of the variant sequence, used for local alignment with reference sequence
    keysite=len(read)-keysite+1


    #Local alignment
    (m_score,m_cigar,start1,end1,mseq1,start2,end2,mseq2,x,y,z,similarity)=SmithWaterman.local_alignment(refseq,query,0,keysite)#Reverse complement of the sequence covering SNV aligns to reference genome sequence
    #-------------------------------------------------------------------------------------------
    #(m_score,m_cigar,start1,end1,mseq1,start2,end2,mseq2,x,y,z)=SmithWaterman.local_alignment(refseq,query,0,0)
    #m_score Alignment score
    #m_cigar Alignment details
    #start1 Alignment start site of extracted forward reference genome sequence
    #end1 Alignment end site of extracted forward reference genome sequence
    #mseq1 Aligned forward reference genome sequence
    #start2 Start site of the reverse complementary sequence of the variant sequence (forward end site)
    #end2 End site of the reverse complementary sequence of the variant sequence (forward start site)
    #mseq2 Reverse complementary sequence of the variant sequence
    #-----------------------------------------------------------------------------


    #Analysis after alignment
    if not int(m_score):#Site without palindrome false positive basis
        return None

    #Map alignment results to reference genome coordinates
    new_start1= pos - flank - 1 +int(start1)#Actual alignment start site on reference genome sequence
    new_end1=pos - flank - 1 + int(end1)#Actual alignment end site on reference genome sequence
    new_start2=len(read)-int(end2)+1 #Start site before reverse complement of read
    new_end2=len(read)-int(start2)+1 #End site before reverse complement of read
    reverse_mseq1=complement_reverse(mseq1)
    reverse_mseq2=complement_reverse(mseq2)
    return [chrname,pos,ref_base,alt_base,m_score,m_cigar,new_start1,new_end1,mseq1,reverse_mseq1,new_start2,new_end2,reverse_mseq2,mseq2,similarity]
def reverse_cigar(cigar):
    '''
    Output reversed cigar
    '''
    new_cigar=''
    pattern=re.compile('((\d+)([SDIM]))')#S here means substitution, CIGAR in BAM has no base substitution
    cigar_parse=re.findall(pattern,cigar)
    cigar_parse.reverse()
    for each_cigar in cigar_parse:
        new_cigar+=each_cigar[0]
    return new_cigar

def present_artifact(match_list):
    '''
    Output alignment results. When used as a module, this function is not executed. When the script is used alone, output alignment results to stdout and mark SNV with underscores
    '''
    (chrname,pos,ref_base,alt_base,m_score,m_cigar,new_start1,new_end1,new_start2,mseq1,new_end2,reverse_query_seq,mseq2,similarity)=match_list
    #ref highlight
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
    #query highlight
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

    #Replace ref type with alt type, alt cover snv
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
#   Read fasta
    ref_dict=read_fasta(fasta)

#   Process snv, whether it is a single site or a file
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
