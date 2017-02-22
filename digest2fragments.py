#!/usr/bin/env python
# -*- coding: utf-8 -*-
## digest2fragments
## 2017 CPCantalapiedra/

import sys

def f_does_overlap(site1, site2):
    overlap = False
    
    if site1[0]!=site2[0]: # if different chromosome, no overlap
        overlap = False
        
    else: # same chromosome
        if long(site2[1]) > long(site1[2]): # if start of site2 > end of site1, no overlap
            overlap = False
        else: # all other cases should be overlap if the BED file was sorted
            overlap = True
    
    return overlap

#

def f_create_fragment(site1, site2):
    fragment = {"chrom":"", "start":-1, "end":-1, "enz1":"", "enz2":""}
    
    if site1:
        chrom = site2[0] # or site1[0], it should be the same
        start = long(site1[2])+1
        end = long(site2[1])-1
        enz1 = site1[3]
        enz2 = site2[3]
    elif site2:
        chrom = site2[0]
        start = 1
        end = long(site2[1])-1
        enz1 = "5'"
        enz2 = site2[3]
    
    fragment["chrom"] = chrom
    fragment["start"] = start
    fragment["end"] = end
    fragment["enz1"] = enz1
    fragment["enz2"] = enz2
    
    return fragment


# Read arguments

sites_file = sys.argv[1]

sys.stderr.write("Digestion sites input file "+sites_file+"\n")

# Print header

sys.stdout.write("\t".join(["chr", "start", "end", "fragtype"])+"\n")

prev_site_data = None
for site in open(sites_file, 'r'):
    
    if site.startswith("names"): continue
    
    site_data = site.strip().split()
    
    if prev_site_data:
        
        if f_does_overlap(prev_site_data, site_data):
            sys.stderr.write("Overlap: "+str(prev_site_data)+" - "+str(site_data)+"\n")
            continue
        
        if site_data[0] == prev_site_data[0]: # if same chromosome
            fragment = f_create_fragment(prev_site_data, site_data)
            
            if fragment["end"]-fragment["start"]>0:
                sys.stdout.write("\t".join([str(x) for x in [fragment["chrom"], fragment["start"], fragment["end"], fragment["enz1"]+":"+fragment["enz2"]]]))
                sys.stdout.write("\n")
        #else:
        #    prev_site_data = None
        #    
        #    fragment = f_create_fragment(prev_site_data, site_data)
        #
        #if fragment["end"]-fragment["start"]>0:
        #    sys.stdout.write("\t".join([str(x) for x in [fragment["chrom"], fragment["start"], fragment["end"], fragment["enz1"]+":"+fragment["enz2"]]]))
        #    sys.stdout.write("\n")
    #else:
    #    fragment = f_create_fragment(prev_site_data, site_data)
    #    
    #    if fragment["end"]-fragment["start"]>0:
    #        sys.stdout.write("\t".join([str(x) for x in [fragment["chrom"], fragment["start"], fragment["end"], fragment["enz1"]+":"+fragment["enz2"]]]))
    #        sys.stdout.write("\n")
    
    prev_site_data = site_data

# It could create fragment from the last digestion point to the end of the chromosome
# but it is just one fragment per chromosome and it would require to know chromosomes sizes
    
sys.stderr.write("finished.\n")

## END