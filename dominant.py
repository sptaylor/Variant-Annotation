#! /usr/bin/python

import os
import sys
import string
import re

file=open("/Users/sptaylor/Dropbox/python/allSamples.withID.nodbsnp.txt")
counter=0

def main():
    result={}
    for line in file.readlines():
        if not line.startswith("#"):
            columns = line.split()
            sample,uploaded_var, location, ref, alt, genotype, consequence, transcript, gene = columns[0], columns[1], columns[2], columns[3], columns[4],columns[5],columns[10], columns[11], columns[12]
            if genotype == '0/1':
                if consequence in ('ESSENTIAL_SPLICE_SITE',
                                   'NMD_TRANSCRIPT,ESSENTIAL_SPLICE_SITE',
                                   'NMD_TRANSCRIPT,NON_SYNONYMOUS_CODING',
                                   'NMD_TRANSCRIPT,SPLICE_SITE,NON_SYNONYMOUS_CODING',
                                   'NON_SYNONYMOUS_CODING',
                                   'SPLICE_SITE,NON_SYNONYMOUS_CODING',
                                   'STOP_GAINED',
                                   'STOP_GAINED,NMD_TRANSCRIPT',
                                   'STOP_LOST',
                                   'WITHIN_MATURE_miRNA',
                                   'WITHIN_NON_CODING_GENE,ESSENTIAL_SPLICE_SITE'):
                    #Assign key value, transcipt=>sample
                    tx=transcript+","+gene
                    if tx not in result:
                        result[tx] = set()
                    result[tx].add(sample)
		    
    #return result
    control='1005'
    print "Transcript\tCase\tControl\tSample_Set"
    candidates = set()
    for tx in result:
        #print len(result[transcript]), transcript
        #if any("1005" in s for s in result[transcript]):
        if control in result[tx]:
            caseCount=len(result[tx])-1
            #print transcript,"\t",len(result[tx])-1,"\t",'1',"\t",",".join(result[tx])
        else:
            caseCount=len(result[tx])
            if caseCount == 5:
                print tx,"\t",len(result[tx]),"\t",'0',"\t",",".join(result[tx])
		transcript = tx.split(',',1)
		candidates.add(transcript[0])
    #return candidates
		
    file2=open("/Users/sptaylor/Dropbox/allSamples.withID.txt")
    out=open('VEP.dominant.5_cases.txt','w')
    for line in file2.readlines():
        if not line.startswith("#"):
            columns = line.split('\t')
            sample,uploaded_var, location, ref, alt, genotype, consequence, transcript, gene, existingvar = columns[0], columns[1], columns[2], columns[3], columns[4],columns[5],columns[10], columns[11], columns[12], columns[14]
            if genotype == '0/1':
                if consequence in ('ESSENTIAL_SPLICE_SITE',
                                   'NMD_TRANSCRIPT,ESSENTIAL_SPLICE_SITE',
                                   'NMD_TRANSCRIPT,NON_SYNONYMOUS_CODING',
                                   'NMD_TRANSCRIPT,SPLICE_SITE,NON_SYNONYMOUS_CODING',
                                   'NON_SYNONYMOUS_CODING',
                                   'SPLICE_SITE,NON_SYNONYMOUS_CODING',
                                   'STOP_GAINED',
                                   'STOP_GAINED,NMD_TRANSCRIPT',
                                   'STOP_LOST',
                                   'WITHIN_MATURE_miRNA',
                                   'WITHIN_NON_CODING_GENE,ESSENTIAL_SPLICE_SITE'):
		    if not existingvar.startswith("rs"):
			if transcript in candidates:
			    print >> out, line

    out.close()
			
			

            
    
        

if __name__ == '__main__':
	main()



                
    

