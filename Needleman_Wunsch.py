#!coding:utf-8

#----------------------------------------------
#Project: Commercial Inspection - Basic Process - Basic Analysis
#Description: Core module of enzyme digestion filtering program - dynamic programming global alignment
#Usage: Not used alone
#Author: Luo Lei
#------------------------------------------------

import sys
#seq1=sys.argv[1]
#seq2=sys.argv[2]

global gap,mismatch,match
#gap=-29 #v0.3
gap=-11 #v0.4
#mismatch=-19 #v0.3
mismatch=-9 #v0.4
match=2


zero_site=[]
def local_alignment(seq1,seq2,keysite1=None,keysite2=None):

    len1=len(seq1)
    len2=len(seq2)
    score_matrix=[]

    #Generate initial scoring matrix
    for i in range(len1+1):
        score_matrix.append([])
        for j in range(len2+1):
            score_matrix[i].append(0)

    #Initialize penalty matrix
    for i in range(len(score_matrix[0])):
        if i==0:continue
        score_matrix[0][i]=score_matrix[0][i-1]+gap
    for i in range(len(score_matrix)):
        if i==0:continue
        score_matrix[i][0]=score_matrix[i-1][0]+gap

    
    #Build scoring matrix
    for i in range(len(score_matrix)):
        if i==0:continue
        for j in range(len(score_matrix[i])):
            if j==0:continue
            if seq1[i-1]==seq2[j-1]:
                score_0=score_matrix[i-1][j-1]+match
            else:
                score_0=score_matrix[i-1][j-1]+mismatch
            score_1=score_matrix[i][j-1]+gap
            score_2=score_matrix[i-1][j]+gap
            score_matrix[i][j]=max([score_0,score_1,score_2])
#    print_matrix(score_matrix)

    #Global alignment backtracking
    m_cigar=Needleman_Wunsch(seq1,seq2,score_matrix)
    #Local alignment backtracking
    return m_cigar
    



def print_matrix(score_matrix):
    for j in range(len(score_matrix[0])):
        tmp=[]
        for i in range(len(score_matrix)):
            tmp.append(str(score_matrix[i][j]))
        print('\t'.join(tmp))

def Needleman_Wunsch(seq1,seq2,score_matrix):
    i=len(seq1)
    j=len(seq2)
    recall=[]
    m_cigar=[]
    while i>0 and j>0:
        score=score_matrix[i][j]#Alignment score at this point
        #Backtrack to upper-left
        if seq1[i-1]==seq2[j-1]:
            if score_matrix[i-1][j-1]+match==score:
                recall.append([0,i,j])
                i-=1
                j-=1
                m_cigar.append('M')
                continue
        if seq1[i-1]!=seq2[j-1]:
            if score_matrix[i-1][j-1]+mismatch==score:
                recall.append([0,i,j])
                i-=1
                j-=1
                m_cigar.append('S')
                continue
        #Backtrack left
        if score_matrix[i-1][j]+gap==score:
                recall.append([-1,i,j])
                i-=1
                m_cigar.append('I')
                continue
        #Backtrack up
        if score_matrix[i][j-1]+gap==score:
                recall.append([1,i,j])
                j-=1
                m_cigar.append('D')
                continue
    m1=[]
    m2=[]
    m_cigar.reverse()
    m_cigar=stat_cigar(m_cigar)
    for i in range(len(recall)-1,-1,-1):
        if recall[i][0]==0:
            m1.append(seq1[recall[i][1]-1])
            m2.append(seq2[recall[i][2]-1])
        if recall[i][0]==-1:
            m1.append(seq1[recall[i][1]-1])
            m2.append('-')
        if recall[i][0]==1:
            m1.append('-')
            m2.append(seq2[recall[i][2]-1])
#    print(''.join(m2))
#    print(''.join(m1))
    return m_cigar

def Smith_Waterman(seq1,seq2,score_matrix):
#    i=len(seq1)
#    j=len(seq2)
    #Find the highest scoring position
    max_score=-100000000000000000000000000
    
    for i in range(len(score_matrix)):
        for j in range(len(score_matrix[i])):
            if score_matrix[i][j]>=max_score:
                max_score=score_matrix[i][j]
#                print(max_score)
                end1=i
                end2=j
#    #Only look at the last column, assuming query sequence is fully consumed
#    for j in range(len(score_matrix[-1])):
#        if score_matrix[-1][j] >= max_score:
#            max_score=score_matrix[-1][j]
#            end1=len(score_matrix)-1
#            end2=j
    total=0
    similar=0
    matchsitelist=[]
    i=end1
    j=end2
    matchsitelist.append([i,j])
    recall=[]
    m_cigar=[]
    k=0
    while i>0 and j>0:
        score=score_matrix[i][j]#Alignment score at this point
        last_score=score_matrix[i-1][j-1]#Upper-left alignment score
        if score==0:break
            
        #Backtrack to upper-left
        if seq1[i-1]==seq2[j-1]:
            if score_matrix[i-1][j-1]+match==score:
                recall.append([0,i,j])
                i-=1
                j-=1
                m_cigar.append('M')
                matchsitelist.append([i,j])
                similar+=1
                total+=1
                continue                
        #Incorrect backtrack to upper-left
        if seq1[i-1]!=seq2[j-1] and k==1:
            if score_matrix[i-1][j-1]+mismatch==score:
                recall.append([0,i,j])
                i-=1
                j-=1
                m_cigar.append('S')
                total+=1
                matchsitelist.append([i,j])
                continue    
        k=1                    
        #Backtrack left
        if score_matrix[i-1][j]+gap==score:
                recall.append([-1,i,j])
                i-=1
                m_cigar.append('I')
                total+=1
                matchsitelist.append([i,j])
                continue
        #Backtrack up
        if score_matrix[i][j-1]+gap==score:
                recall.append([1,i,j])
                j-=1
                m_cigar.append('D')
                total+=1
                matchsitelist.append([i,j])
                continue
        #Incorrect backtrack to upper-left
        if seq1[i-1]!=seq2[j-1]:
            if score_matrix[i-1][j-1]+mismatch==score:
                recall.append([0,i,j])
                i-=1
                j-=1
                m_cigar.append('S')
                total+=1
                matchsitelist.append([i,j])
                continue
        #Highest scoring position eliminated
        if score_matrix[i][j-1]==-1 or score_matrix[i-1][j]==-1 or score_matrix[i-1][j-1]==-1:
            return -1,'',0,0,'',0,0,'',end1,end2,[],0
        break
    m_cigar.reverse()
    if not m_cigar:
        return 0,0,0,0,0,0,0,0,0,0,0,0
    m_cigar=stat_cigar(m_cigar)
#    print(m_cigar)
    start1=i
    start2=j
    m1=[]
    m2=[]
    for i in range(len(recall)-1,-1,-1):
        if recall[i][0]==0:
            m1.append(seq1[recall[i][1]-1])
            m2.append(seq2[recall[i][2]-1])
        if recall[i][0]==-1:
            m1.append(seq1[recall[i][1]-1])
            m2.append('-')
        if recall[i][0]==1:
            m1.append('-')
            m2.append(seq2[recall[i][2]-1])
#    print(''.join(m2))
#    print(''.join(m1))
#    print(max_score)
#    print(score_matrix[i][j])
    m_score=max_score - score_matrix[i][j]
    similarity=float(similar)/total
#    print(str(m_score)+'\t'+m_cigar+'\t'+str(start1+1)+'\t'+str(end1)+'\t'+seq1[start1:end1]+'\t'+str(start2+1)+'\t'+str(end2)+'\t'+seq2[start2:end2])
    return m_score,m_cigar,start1+1,end1,seq1[start1:end1],start2+1,end2,seq2[start2:end2],end1,end2,matchsitelist,similarity

def stat_cigar(m_cigar):
    dic={}
    dic[m_cigar[0]]=1
    s=""
    i=1
    while i<len(m_cigar):
        if m_cigar[i] not in dic:
            tmp=str(dic[m_cigar[i-1]])+m_cigar[i-1]
            s+=tmp
            dic={}
            dic[m_cigar[i]]=1
        else:
            dic[m_cigar[i]]+=1

        i+=1
    tmp=str(dic[m_cigar[i-1]])+m_cigar[i-1]
    s+=tmp
    return s

#dtgh(sys.argv[1],sys.argv[2],int(sys.argv[3]),int(sys.argv[4]))


#local_alignment('ATCGGGCGC','ATCGC')
