#!/home/deakig/usr/local/bin/python
# coding: utf-8
#HirBin function clusterBinsToSubbins

from hirbin import *
from hirbin.parsers import *
from hirbin.thirdparty import runUclust
import os
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
import argparse
import sys
import re
from multiprocessing.dummy import Pool

def parseArgs():
  ''' 
    Function for parsing arguments 
  '''
  parser = argparse.ArgumentParser(description='Cluster annotated sequences (bins) into sub-bins')
  parser.add_argument('-m', dest='mapping_file',help='The path to the mapping/metadata file (required)',required=True)
  parser.add_argument('--type', dest='type', default='protein', help='Sequence type. "prot" if the reference and annotation are for amino acid sequences. "nucl" if nucleotide sequences. In the latter case the nucleotide sequences will be translated. (default:  %(default)s)')
  parser.add_argument('-o', dest='output_dir', help='The name of the output directory. If not specified, a new directory is created')
  parser.add_argument('-id', '--sequenceIdentity', dest='identity',default=0.7,help='Sequence identity cutoff used for sub-binning (default:  %(default)s)')
  parser.add_argument('-p','--minRepresented', dest='p', default=0.75, help='Require non-zero counts in at least the fraction p of the samples for the bin/sub-bin to be representative. (default:  %(default)s)')
  parser.add_argument('--minMeanCount', dest='minMeanCount', default=3, help='Minimum mean count per sample for the bin/sub-bin to be representative. (default:  %(default)s)')
  parser.add_argument('-n', dest='n',help='number of threads (default = %(default)s)',default=1, type=int)
  parser.add_argument('--reClustering', dest='reClustering',action="store_true",help='Perform reclustering and parsing. Use when files for clustering already exists in the hirbin directory. (e.g. rerun clustering with new parameters).')
  parser.add_argument('--onlyParsing', dest='onlyParsing',action="store_true",help='Perform only parsing. Use when files after clustering already exists in the hirbin directory.')
  parser.add_argument('--onlyClustering', dest='onlyClustering',action="store_true",help='Perform only clustering (e.g. clustering with new parameters).')
  parser.add_argument('-f',dest='force_output_directory',action="store_true",help='Force, if the output directory should be overwritten')
  parser.add_argument('--usearchPath',dest='usearchpath',default='usearch',help='The path to the usearch program (default = %(default)s)')
  arguments=parser.parse_args(sys.argv[1:]) 
  return arguments

  
def extract_sequences_one_sample(args):
  '''
    Function for extracting protein sequences given annotated regions and produce fasta files that can be used for clustering
  '''
  (fasta_path,annotation_path,sample, output_dir,domain)=args
  domainInfo=load_annotation_pfam(annotation_path,domain)
  print "Generating domain fasta sequences for "+sample+" ..."
  from Bio import SeqIO
  from Bio.SeqFeature import SeqFeature, FeatureLocation
  from Bio.SeqRecord import SeqRecord 
  (annot,start,stop,strand,evalue)=domainInfo
  record_dict=index_fasta(fasta_path)
  recordlist=[]
  outfilename=output_dir +'/forClustering/'+ sample +'.fasta'
  outhandle=open(outfilename,'w')
  for domainID in annot.keys():
      for i in range(len(annot[domainID])):
          domain=annot[domainID][i]
          try:
              seq=record_dict[domain]
          except KeyError:
             print "Error: " + domain + " not in fasta file.\n"
             break    
          a=start[domainID][i]
          b=stop[domainID][i]
          seq_strand=strand[domainID][i]
          seq_evalue=evalue[domainID][i]
          if seq_strand in '+':
             domain_feature = SeqFeature(FeatureLocation(a-1, b-1), type="domain", strand=1)
          elif seq_strand in '-':
               domain_feature = SeqFeature(FeatureLocation(a-1, b-1), type="domain", strand=-1)
          feature_seq = domain_feature.extract(seq)
          feature_seq.id=feature_seq.id+' '+domainID+' '+seq_evalue 
          recordlist.append(feature_seq)
  SeqIO.write(recordlist, outhandle, "fasta")
  outhandle.close()
  print "Done"

def getSequencesPerDomain(metadata,dom):
  '''
    Function for creating one fasta file per domain (TIGRFAM/PFAM) that can be used as input to the clustering
  '''
  import mmap

  outputdir=metadata.output_directory
  filelist=[outputdir+'/forClustering/'+samplename+".fasta" for samplename in metadata.reference.keys()]
  domains={}

  with open(dom,'rb') as m:
  	f=mmap.mmap(m.fileno(), 0, prot=mmap.PROT_READ)
   	line = f.readline()
  while line:
  	line=line.rstrip('\n')
  	domain=line
	domains[domain]=[]
	line = f.readline()
  m.close()

  for filename in filelist:
    f1=open(filename)
    sequences = SeqIO.parse(f1, "fasta")
    for s in sequences:
      domain=s.description.split()[1]
#      if domain not in domains.keys():
#        domains[domain]=[]
#        domains[domain].append(s)
#      else:
      domains[domain].append(s)
    f1.close()
  print "Loaded seqeunces\n"
 
  for domain in domains:
    g=open(outputdir+'/forClustering/'+domain+'.fasta','w')
    SeqIO.write(domains[domain], g, "fasta")
    g.close()
  for filename in filelist:
    os.remove(filename)


def extractSequences(metadata,n,force):
  '''
    Function for extracting protein sequences given annotated regions and produce fasta files that can be used for clustering
  '''
  outputdir=metadata.output_directory
  try:
    os.mkdir(outputdir+'/forClustering/')
  except OSError as e:
    if not force:
      print "Error: Output directory already exists, you can use an already existing output directory by including the flag -f"
      raise
  domain=metadata.domain.values()
  fasta_list=metadata.reference.values()
  annot_list=metadata.annotation.values()
  sample_list=metadata.annotation.keys()
  #arglist=['BIGWOOD.pep','BIGWOOD.hmmout', 'BIGWOOD',outputdir]
  arglist=[(i,j,sample,outputdir,dom) for (i,j,sample,dom) in zip(fasta_list,annot_list,sample_list,domain)]
  print "Initiating " + str(n) + " processes\n" + str(arglist)
  p=Pool(n)
  p.map(extract_sequences_one_sample,arglist)
  getSequencesPerDomain(metadata,str(domain[0]))
  
def getSubBins(groups,clustpath,cutoffnumber,minMeanCount,identity,countDict):
  '''
    Creating a dictonary with the structure of the sub-bins, i.e. which contigs that are included in each sub-bin.
    The sub-bins not passing the criteria numberOfSamples > cutoffnumber and meanCount>minMeanCount are excluded.
  '''
  samplelisttotal=groups.keys()
  cutoffnumber=float(cutoffnumber)
  cutoff=round(len(samplelisttotal)*cutoffnumber)
  domainlist=os.listdir(clustpath)
  directory="clust"+str(identity) #loop through each ID cutoff
  print "Getting subBins for " + directory
  countsvec={}#initiate counts vector
  countsvec[directory]={} #initiate counts vector
  for fname in domainlist:
    domain=fname.rstrip('.fasta.uc')
    #countsvec[directory][domain]={} #initiate counts vector
    clusters=getClusterStruct(clustpath+'/'+fname)
    for key in clusters.keys(): #for each cluster, check if the cluster is large enough for statistical testing
      #if n>=cutoff: #if the cluster is large enough, get the counts from results file.
      samplelist=[]
      n=0
      for value in clusters[key]:
        contigID=re.sub('_[0-9]$','',value) #remove trailing number
        try:
          crow=countDict[contigID+'_'+domain]      
          for sample in crow:            
            if sample not in samplelist and countDict[contigID+'_'+domain][sample]>0:
              samplelist.append(sample)
              n+=1
        except KeyError:
          #print "No count data found for "+ contigID+'_'+domain
          xxx=1
       #   raise
      if n>=cutoff:
        countsvec[directory][domain]={} #initiate counts vector
        countsvec[directory][domain][key]={}#initiate counts vector
        totcount=0
        for value in clusters[key]:
          contigID=re.sub('_[0-9]$','',value) #remove trailing number
          try:
            crow=countDict[contigID+'_'+domain]
            for sample in crow:
              if sample in countsvec[directory][domain][key].keys():
                countsvec[directory][domain][key][sample]+=countDict[contigID+'_'+domain][sample]
              else:
                countsvec[directory][domain][key][sample]=countDict[contigID+'_'+domain][sample]
                  #check mean count
          except KeyError:
            print "No count data found for "+ contigID+'_'+domain
            # xxx=1
            # raise
        
        tnum=1
        if float(len(samplelist))>0:
          tnum=float(len(samplelist))
        meancount=sum(countsvec[directory][domain][key].values())/tnum  #float(len(samplelist))
        if float(meancount)<float(minMeanCount): #remove if mean count is less than minMeanCount.
          print "Deleted:"+str(countsvec[directory][domain][key])
          del(countsvec[directory][domain][key])
             
  
  #countsvec[directory][domain][key][sample]=samplecount
  
  print "Writing counts to file for " + directory
  g=open(clustpath+'/../abundance_matrix_subbins_'+directory+'.txt','w')
  samplelisttotal=sorted(samplelisttotal)
  for s in samplelisttotal:
      g.write('\t'+s)
  g.write('\n')  
  for domain in countsvec[directory].keys():
      for key in countsvec[directory][domain].keys():
          #print "Domain_value: "+domain+" "+str(key)
          g.write(domain+'_'+directory+'_'+str(key))
          for s in samplelisttotal:
              if s in countsvec[directory][domain][key]:
                 g.write('\t'+str(countsvec[directory][domain][key][s]))
              else:
                   g.write('\t0')
          g.write('\n')
  g.close()
  print "Done"

def createAbundanceMatrixAll(metadata,p,minMeanCount):
  ''' Function for parsing the output from the mapping in TIGRFAM format and creating an abundance matrix'''
  
  filelist=metadata.counts
  samplelist=metadata.samples
  domains={}
  #create sample names
  for sample in filelist: #loop through files and pick out TIGRFAMs + counts
      f=open(filelist[sample])
      for line in f:
          name=line.split('\t')[0]
          count=int(line.split('\t')[1])
          tigrfam=name.split(':')[0]#.split('_')[-1]
          tigrfam=re.sub(r'(.*[0-9]_)([a-zA-Z].*)',r'\2',tigrfam)
          #if tigrfam in str(range(0,10)): #check if TIGRFAM name ends with _1
          #  tigrfam=name.split(':')[0].split('_')[-2]
          if tigrfam not in domains:
             domains[tigrfam]={}
             for i in samplelist:
                 domains[tigrfam][i]=0
          domains[tigrfam][sample] += count
      f.close()
  #filtereddomains=filterAbundanceMatrix(metadata,domains,p,minMeanCount)
  writeAbundanceMatrix(metadata,domains)
  #writeAbundanceMatrix(metadata,filtereddomains)
  return domains

def writeAbundanceMatrix(metadata,domains):
  '''Function for writing the abundance matrix (for bins) to a file in the hirbin output directory'''
  #save results in the file abundance_matrix.txt
  filelist=metadata.counts
  samplelist=metadata.samples
  output_directory=metadata.output_directory
  g=open(output_directory+'/abundance_matrix_bins.txt','w')
  g.write('domain')
  for sample in samplelist:
      g.write('\t'+sample)
  g.write('\n')
  for domain in domains.keys():
      g.write(domain)
      for sample in samplelist:
          g.write('\t'+str(domains[domain][sample]))
      g.write('\n')
  g.close()


def main(mappingFile,output_directory,type,p,minMeanCount,identity,n,onlyClustering,onlyParsing,reClustering,force,usearchpath):
  #reading metadata file and creating output directory
  metadata=Hirbin_run(output_directory)
  metadata.readMetadata(mappingFile)
  output_directory=metadata.createOutputDirectory(output_directory)
  metadata.output_directory=output_directory
  if not reClustering:
  	extractSequences(metadata,n,force)
  try:
    os.mkdir(output_directory+'/clust'+str(identity)) #create a directory for the clustering files
  except OSError as e:
    if not force:
      print "Output directory already exists, you can use an already existing output directory by including the flag -f"
      raise
  if not onlyParsing:
    runUclust(output_directory+"/forClustering/",identity,usearchpath) #run the clustering step
  if not onlyClustering:
    countDict=getCountStruct(metadata) #read the results from the mapping
    getSubBins(metadata.groups,output_directory+"/clust"+str(identity),p,minMeanCount,identity,countDict) #get the results from clustering
    domains=createAbundanceMatrix(metadata,p,minMeanCount) #create abundance matrix - doesn't work(GD)
    domains=createAbundanceMatrixAll(metadata,p,minMeanCount)


if __name__=='__main__':
  arguments=parseArgs()
  if arguments.onlyParsing:
    arguments.reClustering=True
  if arguments.reClustering==True:
    arguments.force_output_directory=True
  main(arguments.mapping_file,arguments.output_dir,arguments.type,arguments.p,arguments.minMeanCount,arguments.identity,arguments.n,arguments.onlyClustering,arguments.onlyParsing,arguments.reClustering,arguments.force_output_directory,arguments.usearchpath)

