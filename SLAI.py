#!coding:utf-8

#----------------------------------------------
#Project: Commercial Inspection - Basic Process - Basic Analysis
#Description: Main program for enzyme digestion filtering, filters false positive sites caused by enzyme digestion, main program is responsible for launching multi-processes
#Usage: python3 SLAI.py -h
#Author: Luo Lei
#------------------------------------------------

#----------------
#v0.3
#Changed read counting to fragment counting
#v0.4
#Hook structure reads require snv position and softclip
#v0.6
#Higher tolerance for softclip and hardclip, stable version
#v0.7
#freq limit set to 40, removed hardclip, version for large-scale testing
#v0.8
#add cpu=1
#os.system changed to subprocess.Popen. subprocess.run only work in python3
#----------------
import time
from multiprocessing import Pool
import os
import sys
import argparse
import judge_artifact as judge_artifact
import subprocess as sub
import shutil
import filter_artifact_snv
import support_reads as SR

def argument_parser():
    parser = argparse.ArgumentParser(description="input bam file and snv file to find out which read support this snv")
    parser.add_argument('--bamfile',help='input bam file',required=True)
    #/beegfs/work/commercial_test/cupcake/databases/gatk_bundle/2.8/hg19/ucsc.hg19.noconfig.fasta
    parser.add_argument('--fasta',help='fasta file used',required=True)
    parser.add_argument('--outdir',help='outdir, default current directory',default='./')
    parser.add_argument('--outfile',help='prefix of outfie',required=True)
    parser.add_argument('--snv',help='snv file',required=True)
    parser.add_argument('--flank',help='flank of snv to match',type=int,default=300)
    parser.add_argument('--cpu',help='cpu core number',type=int,default=8)
    parser.add_argument('--final',help='final output',default=False)
#    parser.add_argument('--do_clean',dest="do_clean",help='If given, clean tem folder',default=False,action="store_true")
    argv=vars(parser.parse_args())
    return argv

def run_filter_artifact(fix_conf_dict,snvfile_f,outfile_f):
    '''
    Function called by multi-threading, actually starts multiple python processes, need to modify conf output and input for multi-processing
    '''
    filter_artifact_snv.main(bam=fix_conf_dict['bamfile'],fasta=fix_conf_dict['fasta'],outdir=fix_conf_dict['outdir'],outfile=outfile_f,snvfile=snvfile_f,flank=fix_conf_dict['flank'])
def split_task(my_conf_dict):
    '''
    Split input snv sites into multiple parts, number specified by --cpu
    '''
    cpu=int(my_conf_dict['cpu'])
    snv=my_conf_dict['snv']
    split_snv_list=[]
    for split_item in range(cpu):
        split_snv_list.append([])
    with open(snv) as r:
        header=r.readline().strip("\n")
        k=0#Control which split unit to add to
        while True:
            line=r.readline()
            if not line:break
            lineInfo=line.strip("\n").split('\t')
            key='\t'.join(lineInfo[:4])
            split_snv_list[k%cpu].append(line.strip("\n"))
            k+=1
    return split_snv_list,header

def write_split_input(my_conf_dict,split_snv_list,header):
    '''
    Write split input files
    '''
    outdir=my_conf_dict['outdir']
    outfile=my_conf_dict['outfile']
    outprefix=outdir+'/'+outfile
    sub.Popen('mkdir -p '+outdir,shell=True).wait()
    inputlist=[]
    for i in range(len(split_snv_list)):
        eachsnvList=split_snv_list[i]
        inputfile=outprefix+'.snv.part'+str(i)
        inputlist.append(inputfile)
        w=open(inputfile,'w')
        w.write(header+'\n')
        for eachsnv in eachsnvList:
            w.write(eachsnv+'\n')
        w.close()
    return inputlist




if __name__=='__main__':
    print("Program starts at: "+time.strftime('%Y-%m-%d %H:%M:%S'))
    argv=argument_parser()
    conf_dict=argv.copy()

    #Modify parameters and specify input/output files
    ini_outdir=conf_dict['outdir']
    ini_outfile=conf_dict['outfile']
    finaloutput=conf_dict['final']
    match_out=conf_dict['outdir']+'/'+conf_dict['outfile']+".artifact_match.txt"
    if finaloutput:
        out_final_snv_report=finaloutput
    else:
        out_final_snv_report=conf_dict['outdir']+'/'+conf_dict['outfile']+".final.artifact_filter.txt"
    fix_match_out=conf_dict['outdir']+'/'+conf_dict['outfile']+".fix_artifact_match.txt"
    support_out=conf_dict['outdir']+'/'+conf_dict['outfile']+".snv_supportReads.txt"
    outfile_judge=ini_outdir+'/'+ini_outfile+'.snv.filter_artifact_snv.txt'
#   if not os.path.exists(support_out) or not os.path.exists(match_out):
    threads=argv['cpu']
    if True:
        if threads>1:#Use multiprocess for multi-threading when cpu is not 1
            conf_dict['outdir']=argv['outdir']+'/tem'
            (split_snv_list,header)=split_task(conf_dict)
            inputlist=write_split_input(conf_dict,split_snv_list,header)
            task_list=[]
            pool=Pool(threads)
            for eachinput in inputlist:
                snvfile=eachinput
                outfile=os.path.basename(eachinput)
                task=pool.apply_async(run_filter_artifact,args=(conf_dict,snvfile,outfile))
                task_list.append(task)
            pool.close()
            pool.join()
            for task in task_list:
                 assert not task.get()
#                p = Process(target=run_filter_artifact, args=(conf_dict,))
#                task_list.append(p)
#                p.start()
#                print('job start')
#                time.sleep(3)
#            for each_task in task_list:
#                each_task.join()
            
            print('over')
            merge_support_reads_cmd='''awk 'NR==1||FNR>1' %s/tem/%s.snv.part*.snv_supportReads.txt >%s'''%(ini_outdir,ini_outfile,support_out)
            merge_match_results_cmd='cat %s/tem/%s.snv.part*_artifact_match.txt >%s' %(ini_outdir,ini_outfile,match_out)
            sub.Popen(merge_support_reads_cmd,shell=True).wait()
            sub.Popen(merge_match_results_cmd,shell=True).wait()
        else:#When cpu is 1, run directly without multi-threading
            conf_dict['outdir']=argv['outdir']
            conf_dict['outfile']=argv['outdir']+'/'+argv['outfile']
            conf_dict['snv']=argv['snv']
            snvfile=conf_dict['snv']
            outfile=ini_outdir+'/'+ini_outfile
            run_filter_artifact(conf_dict,snvfile,outfile)
#            sub.Popen('mv '+conf_dict['outfile']+'_artifact_match.txt '+match_out,shell=True).wait()
    
    #judgement
    inputfile_for_judge=match_out
    snvfile_for_judge=support_out
    (init_judge_dict,match_dict)=judge_artifact.snv_init_judge_dict(inputfile_for_judge,snvfile_for_judge)
    judge_dict=judge_artifact.snv_judge_dict(init_judge_dict,match_dict,fix_match_out,"panel15_pro")#Core working function
    judge_artifact.write_output(outfile_judge,judge_dict,snvfile_for_judge)#Output
    judge_artifact.write_final_output(out_final_snv_report,judge_dict,snvfile_for_judge)#Output
#    if argv['do_clean'] and os.path.exists(os.path.join(ini_outdir, 'tem')):
#        shutil.rmtree(os.path.join(ini_outdir, 'tem'))
    print("Program ends at: "+time.strftime('%Y-%m-%d %H:%M:%S'))
