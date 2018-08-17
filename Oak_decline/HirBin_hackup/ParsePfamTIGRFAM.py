from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
import mmap

def load_annotation_pfam(filename,dom):
    '''function for loading annotation and coordinates from PFAM and TIGRfam
    (HMMersearch domtblout format) and saving into dictionaries'''

    #list the annotation files in annotation directory
    print "loading annotation from file " + filename
    annot={}
    start={}
    stop={}
    strand={}
    evalue={}
    
    with open(dom,'rb') as m:
    	f=mmap.mmap(m.fileno(), 0, prot=mmap.PROT_READ)
    	line = f.readline()
   	while line:
   		line=line.rstrip('\n')
   		domain=line
		annot[domain]=[]
		start[domain]=[]
		stop[domain]=[]
		strand[domain]=[]
		evalue[domain]=[]   	
   		line = f.readline()
    m.close()
    print "Read in domains "
    with open(filename, 'rb') as m:
	    #open the file
	    f=mmap.mmap(m.fileno(), 0, prot=mmap.PROT_READ)
	    print "File read into memory"
	    line = f.readline()
	    while line:
		line=line.rstrip('\n')
		vec=line.split()
		domain=vec[0] #domain ID
		#print "domain " + domain
		#if domain not in "NA":
		contigID=vec[1] #sequence id
		#print "contig " + contigID
		domainevalue=vec[5] #e-value
		startcoord=int(vec[2]) #internal coordinates
		stopcoord=int(vec[3]) #internal coordinates
		strandinfo=vec[4]    #all features are on sense strand
		#check if domain exists in dictionary and add fields
		#if domain not in annot.keys():

		#	annot[domain].append(contigID)
		#	start[domain].append(startcoord)
		#	stop[domain].append(stopcoord)
		#	strand[domain].append(strandinfo)
		#	evalue[domain].append(domainevalue)
		#else:
		annot[domain].append(contigID)
		start[domain].append(startcoord)
		stop[domain].append(stopcoord)
		strand[domain].append(strandinfo)
		evalue[domain].append(domainevalue)
		line = f.readline()
    m.close()
    print "Number of motifs loaded: "+str(len(annot.keys()))
    #return annotation information as a tuple.
    output= (annot,start,stop,strand,evalue)
    return(output)