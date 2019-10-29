#!/usr/bin/env python2.7

import os,sys
import json
import argparse
import subprocess
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from distutils.spawn import find_executable

dirpath = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
blast_cmd = 'blastn'
db_prophage = dirpath+'/database/PHASTER_prophage_database.fna'
db_plasmid = dirpath+'/database/Refseq_plasmid_database.fna'
log_file = 'MGE_pipeline.log'

# Parse the input arguments.
def parse_args():
  parser = argparse.ArgumentParser(usage='mge_pipeline.py -i <assembly_path> -o <output_path> -t <threads> -c <minimum_coverage> -p <minimum_ident>')
  parser.add_argument("-i",help="<string>: path to input assembly")
  parser.add_argument("-o",help="<string>: path to output directory")
  parser.add_argument("-t",default=1,help="<string>: threads")
  parser.add_argument("-c",default=50,help="<int>: minimum coverage of each hit for MGE verification, default 50")
  parser.add_argument("-p",default=90,help="<int>: minimum identical percentage of each hit for MGE verification, default 90")
  parser.add_argument("--db_prophage",type=os.path.abspath,help="<string>: path to prophage database")
  parser.add_argument("--db_plasmid",type=os.path.abspath,help="<string>: path to plasmid database")
  parser.add_argument("--check",action="store_true",help="<flag>: use '--check' flag to check the required dependencies")
  parser.add_argument("--rename",action="store_true",help="<flag>: use '--rename' flag to rename the contigs before running MGE detection")
  return parser.parse_args()

# check dependencies
check_dependencies = parse_args().check
dependencies = ['seqtk','prokka','blastn','ProphET_standalone.pl','gff_rewrite.pl','mob_recon','phigaro']
if check_dependencies:
    for item in dependencies:
        ext_path = find_executable(item)
        if ext_path is not None:
            print ("Using "+item+" - "+ext_path)
        else:
            print ("ERROR: can not find "+item+" in PATH")
    sys.exit()
# end of --check

# index database
parse_db_prophage = parse_args().db_prophage
if parse_db_prophage:
  db_prophage = parse_db_prophage
if os.path.exists(db_prophage+'.nal') or os.path.exists(db_prophage+'.nhr'):
  pass
else:
  subprocess.check_call('makeblastdb -in '+db_prophage+' -dbtype nucl',shell=True)

parse_db_plasmid = parse_args().db_plasmid
if parse_db_plasmid:
  db_plasmid = parse_db_plasmid
if os.path.exists(db_plasmid+'.nal') or os.path.exists(db_plasmid+'.nhr'):
  pass
else:
  subprocess.check_call('makeblastdb -in '+db_plasmid+' -dbtype nucl',shell=True)

# Function for genome annotation using PROKKA.
def prokka(input_path,filestr,filename,output_path,threads,prokka_fna,prokka_gff):
  rename_contigs = parse_args().rename
  if rename_contigs:
    ## Rename contigs.
    subprocess.check_call('seqtk rename '+input_path+' '+filename+'_ > '+filestr,shell=True)
  else:
    subprocess.check_call('cp '+input_path+' '+filestr,shell=True)
  subprocess.check_call('prokka --mincontiglen 1 --force --kingdom Bacteria --addgenes --prefix '+filename+' --outdir '+filename+'_prokka --cpus '+threads+' '+filestr+' >> '+log_file+' 2>&1',shell=True)
  subprocess.check_call('cp '+filename+'_prokka/'+filename+'.fna '+prokka_fna,shell=True)
  subprocess.check_call('cp '+filename+'_prokka/'+filename+'.gff '+prokka_gff,shell=True)

# Function for prophage detection using ProphET.
def prophet(filename,prokka_fna,prokka_gff):
  prophet_gff=filename+'_prophet.gff'
  prophet_out=filename+'_prophet'
  ## Change the format of the GFF file.
  subprocess.check_call('gff_rewrite.pl --input '+prokka_gff+' --output '+prophet_gff+' --add_missing_features'+' >> '+log_file+' 2>&1',shell=True)
  subprocess.check_call('ProphET_standalone.pl --fasta_in '+prokka_fna+' --gff_in '+prophet_gff+' --outdir '+prophet_out+' >> '+log_file+' 2>&1',shell=True)
  subprocess.check_call('cp '+prophet_out+'/phages_coords '+prophet_out+'.txt',shell=True)

# Function for prophage detection using Phigaro.
def phigaro(filename,prokka_fna,threads,phigaro_contig):
  fasta_tmp=prokka_fna+'.tmp'
  sequences_to_delete = []
  records_to_save = []
  for record in SeqIO.parse(prokka_fna, "fasta"):
    if len(record) < 20000:
      sequences_to_delete.append(record.id)
    else:
      records_to_save.append(record)
  print('These sequencess are shorter than 20000 bp. Phigaro is running without them.')
  print('\n'.join(sequences_to_delete))
  SeqIO.write(records_to_save, fasta_tmp, "fasta")
  phigaro_out=filename+'_phigaro.txt'
  subprocess.check_call('phigaro --not-open -e txt -f '+fasta_tmp+' -t '+threads+' -o '+phigaro_out+' >> '+log_file+' 2>&1',shell=True)

# Function for plasmid dection using MOB-suite.
def mobsuite(filename,prokka_fna,threads):
  mobsuite_dir=filename+'_mobsuite'
  subprocess.check_call('mob_recon --infile '+prokka_fna+' --outdir '+mobsuite_dir+' -n '+threads+' >> '+log_file+' 2>&1',shell=True)
  subprocess.check_call('cp '+mobsuite_dir+'/contig_report.txt '+mobsuite_dir+'.txt',shell=True)

# Function for extracting prophage sequences.
def extract_seq(fasta_file,contig_id,beg,end):
  genome=SeqIO.index(fasta_file,'fasta')
  sequence = genome[contig_id].seq[int(beg)-1:int(end)]
  fasta='mge_sequence.fasta.tmp'
  handle=open(fasta,'w')
  handle.write('>'+contig_id+'_'+beg+'-'+end+'\n')
  handle.write(str(sequence)+'\n')
  handle.close()
  return fasta

# Function for Blastn.
def blastn(blast_cmd,query,db,blastout,threads,max_hsps,max_target_seqs,perc_qcov,perc_ident):
  evalue=1e-5
  blast_fmt = "'6 qseqid stitle pident length qcovhsp qlen slen qstart qend sstart send evalue bitscore'"
  if db=='Prophage':
    blast_db=db_prophage
    blast_database='PHASTER'
  elif db=='Plasmid':
    blast_db=db_plasmid
    blast_database='Refseq_plasmid'
  blastn_out=NcbiblastnCommandline(cmd=blast_cmd,query=query,db=blast_db,evalue=evalue,outfmt=blast_fmt,out=blastout,num_threads=threads,max_hsps=max_hsps,max_target_seqs=max_target_seqs,qcov_hsp_perc=perc_qcov,perc_identity=perc_ident)
  stdout, stderr = blastn_out()
  if os.path.getsize(blastout) > 0:
    blast_result=open(blastout,'r').readline()
    items_blast=blast_result.split('\t')
    blast_subject=items_blast[1]
    blast_ident=items_blast[2]
    blast_length=items_blast[3]
    blast_cov=items_blast[4]
    blast_slength=items_blast[6]
  else:
    blast_database=blast_subject=blast_ident=blast_length=blast_cov=blast_slength='NA'
  return blast_database,blast_subject,blast_ident,blast_length,blast_cov,blast_slength

def main():
  args = parse_args()
  input_path = os.path.abspath(args.i)
  output_path = os.path.abspath(args.o)
  threads = args.t
  qcov=args.c
  ident=args.p

  subprocess.check_call('mkdir -p '+output_path,shell=True)
  os.chdir(output_path)
  filename=os.path.basename(input_path).rsplit('.',1)[0]
  filestr=os.path.split(input_path)[1]
  prokka_fna=filename+'_prokka.fna'
  prokka_gff=filename+'_prokka.gff'
  output_file=filename+'_MGE_results.txt'
  blastout=filename+'_blastout.tmp'
  handle_results=open(output_file+'.tmp','w')
  handle_results.write('Contig_ID\tContig_length\tStart_pos\tEnd_pos\tMGE_length\tMGE_type\tMethod\tDatabase\tBlast_hit\tHit_identity\tAlignment_length\tQuery_coverage\tHit_length\n')

  ## Run Prokka.
  print 'Running Prokka...'
  prokka(input_path,filestr,filename,output_path,threads,prokka_fna,prokka_gff)
  prokka_fna_file=SeqIO.index(prokka_fna,'fasta')

  ## Run ProphET.
  print 'Running ProphET...'
  prophet(filename,prokka_fna,prokka_gff)
  mgetype='Prophage'
  method='ProphET'
  prophet_out=filename+'_prophet.txt'
  prophet_data=open(prophet_out,'r').readlines()
  for line in prophet_data:
    items=line.strip().split('\t')
    contig_id=items[0]
    beg=items[2]
    end=items[3]
    contig_length=str(len(prokka_fna_file[contig_id].seq))
    query=extract_seq(prokka_fna,contig_id,beg,end)
    blast_database,blast_subject,blast_ident,blast_length,blast_cov,blast_slength=blastn(blast_cmd,query,mgetype,blastout,threads,1,1,qcov,ident)
    handle_results.write(contig_id+'\t'+contig_length+'\t'+beg+'\t'+end+'\t'+str(int(end)-int(beg)+1)+'\t'+mgetype+'\t'+method+'\t'+blast_database+'\t'+blast_subject+'\t'+blast_ident+'\t'+blast_length+'\t'+blast_cov+'\t'+blast_slength+'\n')

  ## Run Phigaro.
  print 'Running Phigaro...'
  phigaro_contig='~/.phigaro/config.yml'
  phigaro(filename,prokka_fna,threads,phigaro_contig)
  mgetype='Prophage'
  method='Phigaro'
  phigaro_out=filename+'_phigaro.txt'
  phigaro_data=open(phigaro_out,'r').readlines()[1:]
  for line in phigaro_data:
    items=line.strip().split('\t')
    contig_id=items[0]
    beg=items[1]
    end=items[2]
    contig_length=str(len(prokka_fna_file[contig_id].seq))
    query=extract_seq(prokka_fna,contig_id,beg,end)
    blast_database,blast_subject,blast_ident,blast_length,blast_cov,blast_slength=blastn(blast_cmd,query,mgetype,blastout,threads,1,1,qcov,ident)
    handle_results.write(contig_id+'\t'+contig_length+'\t'+beg+'\t'+end+'\t'+str(int(end)-int(beg)+1)+'\t'+mgetype+'\t'+method+'\t'+blast_database+'\t'+blast_subject+'\t'+blast_ident+'\t'+blast_length+'\t'+blast_cov+'\t'+blast_slength+'\n')

  ## Run MOB-suite.
  print 'Running MOB-suite...'
  mobsuite(filename,prokka_fna,threads)
  mgetype='Plasmid'
  method='MOB-suite'
  mobsuite_out=filename+'_mobsuite.txt'
  mobsuite_data=[x for x in open(mobsuite_out,'r').readlines()[1:] if x.split('\t')[1]!='chromosome']
  for line in mobsuite_data:
    items=line.strip().split('\t')
    contig_id=items[2].split('|')[-1]
    beg='1'
    end=items[3]
    contig_length=str(len(prokka_fna_file[contig_id].seq))
    query=extract_seq(prokka_fna,contig_id,beg,end)
    blast_database,blast_subject,blast_ident,blast_length,blast_cov,blast_slength=blastn(blast_cmd,query,mgetype,blastout,threads,1,1,qcov,ident)
    handle_results.write(contig_id+'\t'+contig_length+'\t'+beg+'\t'+end+'\t'+str(int(end)-int(beg)+1)+'\t'+mgetype+'\t'+method+'\t'+blast_database+'\t'+blast_subject+'\t'+blast_ident+'\t'+blast_length+'\t'+blast_cov+'\t'+blast_slength+'\n')

  handle_results.close()
  subprocess.check_call('(head -n 1 '+output_file+'.tmp'+' && tail -n +2 '+output_file+'.tmp'+' | sort -V) > '+output_file,shell=True)
  subprocess.check_call('rm *tmp error.log',shell=True)
  subprocess.check_call('rm *_prokka.fna *_prokka.gff *_prophet.gff',shell=True)
  subprocess.check_call('rm -r *_prokka',shell=True)
  print 'Done.'

if __name__ == '__main__':
  main()
