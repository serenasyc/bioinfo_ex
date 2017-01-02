#!/usr/bin/env python
import sys
import argparse
import os
import time
import requests
from csv import DictReader
import pandas as pd
import numpy as np
import re

class FusionEvent:
    def __init__(self, cid, chrm1, chrm2, bkpt1, bkpt2, strd1, strd2, ftype, mid, seq):
        '''
        Attributes correspond to the following columns of deFuse esults:
        cid = cluster_id
        chrm1   = gene_chromosome1
        chrm2   = gene_chromosome2
        bkpt1   = genomic_break_pos1
        bkpt2   = genomic_break_pos2 
        strd1   = genomic_strand1
        strd2   = genomic_strand2
        ftype   = X in { 'deletion', 'inversion', 'eversion', 'interchromosomal' }
                          where value of X is 'Y'
        splitr_seq  = splitr_sequence
        mid     = cluster_id of matching event in other file
        '''
        self.cid  = cid
        self.chrm1 = chrm1
        self.chrm2 = chrm2
        self.bkpt1 = bkpt1
        self.bkpt2 = bkpt2
        self.strd1 = strd1
        self.strd2 = strd2
        self.ftype = ftype
        self.splitr_seq = seq
        self.mid = mid
        self.comment = ''

def ensembl_mapping(asm_one, asm_two, chrm, start_pos, end_pos, strand=0):
    '''
    Uses the ensembl assembly converter: 
    http://rest.ensembl.org/documentation/info/assembly_map
    @param asm_one  - ensembl version of input assembly (ex. GRCh37)
    @param asm_two  - ensembl version of output assembly (ex. GRCh38)
    @param chrm     - chromosome
    @param start_pos - start position in chromosome
    @param end_pos  - end position in chromosome
    @param strand   - strand of region

    @return list    - [mapped_start_pos, mapped_end_pos]
    '''
    if end_pos == None:
        end_pos = start_pos

    args = (chrm, start_pos, end_pos, strand)
    server = 'http://rest.ensembl.org'
    ext = '/map/human/{}/%s:%s..%s:%s/{}?'.format(asm_one, asm_two)%args
    r = requests.get(server+ext, headers={'Content-Type' : 'application/json'})

    if not r.ok:
        r.raise_for_status()
        sys.exit()
    decode = r.json()

    N = len(decode['mappings'])
    if N != 0:      
        if N != 1:
            print('\tmultiple mappings ({}) for %s:%s..%s:%s\t'.format(N)%args)
            temp = ','.join(map(str, [x['mapped']['start'] for x in 
                                      decode['mappings']]))
            print('\t\t{} --> {}'.format(start_pos, temp))
        return [decode['mappings'][0]['mapped']['start'],  
                decode['mappings'][0]['mapped']['end']]
    else:           # no result
        print('\tno mapping found\t%s:%s..%s:%s'%args)
        return [-1,-1]


def same_fusion_event(f,g):
    '''
    @param f - a FusionEvent object
    @param g - a FusionEvent object

    Assuming you know the two FusionEvent objects have the same chromosomes and
    breakpoints, do they have the same strands and fusion type?
    @return bool
    '''
    return (f.strd1 == g.strd1 and f.strd2 == g.strd2 and f.ftype == g.ftype)

def read_tsv(tsv, ver, sample_id):
    '''
    @param tsv - deFuse results.filtered.tsv
    @param ver - Ensembl version used to run deFuse (ex. GRCh38)

    Convert breakpoint coordinates to GRCh38 if in GRCh37
    Returns a dictionary of chromosome_pairs : { breakpoint_pairs : FusionEvent }
    @return dict
    '''
    filename = tsv.split('/')[-1]
    mapped = True if 'mapped' in filename else False
    print('\nReading {} and calling fusion types ... '.format(tsv))
    ddf = {}
    duplicates = []
    N = 0
    with open(tsv, 'r') as tsvfile:
        reader = DictReader(tsvfile, delimiter='\t')
        for r in reader:
            # Unpacking values
            cid = r['cluster_id']
            chrm1 = r['gene_chromosome1']
            chrm2 = r['gene_chromosome2']
            bkpt1 = r['genomic_break_pos1']
            bkpt2 = r['genomic_break_pos2']
            strd1 = r['genomic_strand1']
            strd2 = r['genomic_strand2']
            # Convert breakpoint if from GRCh37
            if ver == 'GRCh37':
                mapped_bkpt1 = ensembl_mapping('GRCh37', 'GRCh38', chrm1, bkpt1,
                                               None,
                                               '-1' if strd1 == '-' else '1')[0]
                mapped_bkpt2 = ensembl_mapping('GRCh37', 'GRCh38', chrm2, bkpt2,
                                               None,
                                               '-1' if strd2 == '-' else '1')[0]
                # if no mapping available
                if (mapped_bkpt1 == -1 or mapped_bkpt2 == -1):
                    continue
                else:
                    bkpt1 = str(mapped_bkpt1)
                    bkpt2 = str(mapped_bkpt2)

            # Determine fusion type
            if mapped:
                ftype = r['fusion_type']
            else:
                for f in ['deletion', 'inversion', 'eversion', 'interchromosomal']:
                    if r[f] == 'Y':
                        ftype = f
                        break
            # Create FusionEvent
            f = FusionEvent(cid, chrm1, chrm2, int(bkpt1), int(bkpt2), strd1, strd2,
                            ftype, '', r['splitr_sequence'])

            # Add FusionEvent to dictionary ddf
            cp = '{},{}'.format(chrm1, chrm2)
            bp = '{},{}'.format(bkpt1, bkpt2)
            if cp not in ddf:
                ddf[cp] = {}
            if bp not in ddf[cp]:
                ddf[cp][bp] = f
                N += 1
            # FusionEvents not unique
            # event can be supported by different spanning reads clusters
            else:
                g = ddf[cp][bp]
                if not same_fusion_event(f, g):
                    print('Contradicting data for {}:{} | {}:{}'.format(chrm1, bkpt1, chrm2, bkpt2))
                    print('\t({},{},{}) vs ({},{},{})'.format(strd1, strd2, ftype, g.strd1, g.strd2, g.ftype))
                else:
                    duplicates.append('{}:{} | {}:{}'.format(chrm1, bkpt1, chrm2, bkpt2))
    if len(duplicates) != 0:
        print('{} duplicate events found:'.format(len(duplicates)))
        for i in duplicates:
            print('\t{}'.format(i))
    # Save mapped coordinates of the GRCh37 deFuse results
    if ver == 'GRCh37':
        mapfile = '{}.hg38mapped.{}'.format(sample_id, filename)
        f = open(mapfile, 'w')
        headers = [ 'cluster_id',
                    'gene_chromosome1', 'genomic_break_pos1', 'genomic_strand1',
                    'gene_chromosome2', 'genomic_break_pos2', 'genomic_strand2',
                    'fusion_type','match_id', 'splitr_sequence', 'comment']
        f.write('\t'.join(headers)+'\n')
        for cp in ddf.keys():
            for bp in ddf[cp].keys():
                event = ddf[cp][bp]
                row = [event.cid,
                       event.chrm1, str(event.bkpt1), event.strd1,
                       event.chrm2, str(event.bkpt2), event.strd2,
                       event.ftype, event.mid, event.splitr_seq, event.comment]
                f.write('\t'.join(row)+'\n')
        f.close()
    print('Found {} unique events.'.format(N))
    return ddf


def align_breakpoints(s,t):
    '''
    @param s - a sequence 
    @param t - another sequence 
    given the 5 nucleotides before and after the breakpoint of each sequence
    @return int
        if n is negative, move t left by n  (move breakpoint right by n)
        if n is positive, move t right by n (move breakpoint left by n)
        if n is 0, the alignment score is too low (< 6)
    '''
    n = len(s)
    m = len(t)
    M = np.zeros([n,m], dtype =int)
    for j in range(m):      # j = column index
        for i in range(n):  # i = row index
            if i > 0 and j > 0:
                M[i][j] = M[i-1][j-1] + 1 if s[i] == t[j] else 0
    maxscore = M.max()  # length of longest match
    if abs(maxscore) < 25:
        return 0
    maxindex = np.argmax(M)
    maxi = int(maxindex/n)
    maxj = maxindex%m
    return maxi-maxj


def get_breakpoint_region(s1, s2, z=25):
    '''
    @param s1   - sequence with a breakpoint ('|')
    @param s2   - another sequence a breakpoint ('|')
    @param z    - length of region on each side of the breakpoint
    @return str,str -   z of each side of breakpoint for s1
                        omitting the '|'
    '''
    p1 = s1.index('|')
    p2 = s2.index('|')
    s = s1[p1-z:p1]+s1[p1+1:p1+z+1]
    t = s2[p2-z:p2]+s2[p2+1:p2+z+1]
    return s,t


def overlap_seq(s1, s2):
    '''
    @param s1 - a sequence
    @param s2 - another sequence
    @return str,str - the overlap region and a comment
    '''
    o = ''
    p1 = s1.index('|')
    p2 = s2.index('|')

    z = 25
    # Check if the region surrounding the breakpoint is identical
    s,t = get_breakpoint_region(s1,s2,25)
    if s != t:
        bkpt_diff = align_breakpoints(s,t)
        # The splitr_sequences are different
        if bkpt_diff == 0:
            return '','NOT A MATCH'
        # The breakpoints are different but the splitr_sequences are identical
        elif s1[:p1]+s1[p1+1:] == s2[:p2]+s2[p2+1:]:
            return s1 , ', identical splitr_seq'.format(bkpt_diff)
        # The breakpoints are different but the splitr sequences are NOT identical but overlapping
        else:
            s2left = s2[:p2]
            s2right = s2[p2+1:]
            # shift the breakpoint by the difference calculated
            if bkpt_diff > 0:
                s2 = s2left[:bkpt_diff] + '|' + s2left[len(s2left)-bkpt_diff:] + s2right
            else:
                s2 = s2left + s2right[:-bkpt_diff] + '|' + s2right[-bkpt_diff:]
    is_exact = True
    # Find the overlapping sequence
    # first the right side
    for i in range(1, min(len(s1)-p1,len(s2)-p2)):
        if s1[p1+i] == s2[p2+i]:
            o = o + s1[p1+i]
        else:
            is_exact = False
            break
    # then the left side
    for i in range(0, min(p1,p2)):
        if s1[p1-i] == s2[p2-i]:
            o = s1[p1-i] + o
        else:
            is_exact = False
            break
    comment = ', identical splitr_seq' if is_exact else ', overlapping splitr_seq'
    return o, comment


def get_opp_strand(s):
    '''
    @param s - a DNA sequence
    @return the reverse strand of s
    '''
    s = ''.join(reversed(s))
    d = {'A':'T', 'T':'A', 'C':'G', 'G':'C', '|':'|'}
    o = ''
    for i in s:
        o += d[i]
    return o    


def dictionary_cleanup(ddfA, ddfB, to_deleteA, to_deleteB=None):
    '''
    @param ddfA - dictionary of fusion events from fileA
    @param ddfB - dictionary of fusion events from fileB
    @param delete - dictionary of chromosome_pairs : breakpoint_pairs to delete from ddfA and ddfB
    '''
    if to_deleteB == None:
        to_deleteB = to_deleteA
    for cp in to_deleteA:
        for bp in to_deleteA[cp]:
            del ddfA[cp][bp]
        if len(ddfA[cp]) == 0:
            del ddfA[cp]
    for cp in to_deleteB:
        for bp in to_deleteB[cp]:
            del ddfB[cp][bp]
        if len(ddfB[cp]) == 0:
            del ddfB[cp]


def overlay(sample_id, fileA, fileB, ddfA, ddfB, matchfile):
    '''
    @param sample_id - identifier
    @param fileA - results.filtered.tsv from deFuse
    @param fileB - another result.filtered.tsv from deFuse
    @param ddfA - dictionary of chromosome_pairs : { breakpoint_pairs :  FusionEvents } from fileA
    @param ddfB - dictionary of chromosome_pairs : { breakpoint_pairs :  FusionEvents } from fileB
    @param matchfile - name of file containing matched events
    @param nonmatchfile - name of file containing non-matched events
    If events match (same chromosome pairs, same breakpoint pairs, same FusionEvents)
        write to matchfile and remove from ddfA and ddfB
    '''
    print('\nOverlaying results ...')
    matches = []
    to_delete = {}  
    chrm_pairs = set(ddfA.keys()).intersection(set(ddfB.keys()))
    for cp in chrm_pairs:
        bkpt_pairs = set(ddfA[cp].keys()).intersection(set(ddfB[cp].keys()))
        for bp in ddfB[cp].keys():
            # if they match
            if bp in bkpt_pairs and same_fusion_event(ddfA[cp][bp], ddfB[cp][bp]): 
                seq1 = ddfA[cp][bp].splitr_seq
                seq2 = ddfB[cp][bp].splitr_seq
                o, comment = overlap_seq(str(seq1), str(seq2))
                if len(o) == 0:
                    print("WARNING: for cid19 = {} and cid38={}, breakpoints match but split reads do not.".format(ddB[cp][bp].cid, ddfA[cp][bp].cid))
                ddfB[cp][bp].comment += ', {}'.format(comment)
                #Split read sequences are usually ~400 bp long
                #It's possible that only one side of breakpoint region is matched?
                if len(o) < 100 or o.index('|') < 25:
                    print("WARNING: for cid19={} and cid38={}, split read sequences are poorly matched, overlap length of {}".format(ddfB[cp][bp].cid, ddfA[cp][bp].cid, len(o)-19))
                ddfB[cp][bp].splitr_seq = o
                ddfB[cp][bp].mid = ddfA[cp][bp].cid
                matches.append(ddfB[cp][bp])
                if cp not in to_delete:
                    to_delete[cp] = []
                to_delete[cp].append(bp)
 
    print('\t# matches:\t{}'.format(len(matches)))

    headers = [ 'cluster_id',
                'gene_chromosome1', 'genomic_break_pos1', 'genomic_strand1', 
                'gene_chromosome2', 'genomic_break_pos2', 'genomic_strand2', 
                'fusion_type','match_id', 'splitr_sequence', 'comment']

    #If matchfile already exists, it's overwritten
    open(matchfile, 'w').close()
    with open(matchfile, 'a') as f:
        f.write('## Overlayed {} and {}\n'.format(fileA, fileB))
        f.write('## Converted coordinates to GRCh38\n')
        f.write('##File generated on {}\n'.format(time.strftime('%B %d, %Y')))
        f.write('\t'.join(headers)+'\n')

        for e in matches:
            fields = [e.cid, e.chrm1, e.bkpt1, e.strd1, e.chrm2,e.bkpt2, 
                      e.strd2, e.ftype, e.mid, e.splitr_seq, e.comment]
            f.write('\t'.join(map(str,fields))+'\n')
    print('Wrote matched events to {}'.format(matchfile))

    # Cleanup ddfA and ddfB (remove matched events)
    dictionary_cleanup(ddfA,ddfB,to_delete, None)


def check_breakpoint_match(f, g, swap):
    '''
    @param f    - FusionEvent
    @param g    - FusionEvent
    @param swap - bool is it a swap?
    @param bool, bool, str - is_match?, is_swapped?, comment
    '''
    fbkpt1,fbkpt2 = int(f.bkpt1),int(f.bkpt2)
    gbkpt1,gbkpt2 = int(g.bkpt1),int(g.bkpt2)
    
    if swap:
        if fbkpt1 == gbkpt2 and fbkpt2 == gbkpt1:
            return True, ', genes swapped'    
        # are the breakpoints off by a few nucleotides and swapped?
        d_bkpt1 = fbkpt1 - gbkpt2
        d_bkpt2 = fbkpt2 - gbkpt1
        if abs(d_bkpt1) < 1000 and abs(d_bkpt2) < 1000:
            comment = ', genes swapped, bkpt1 differs by {}, bkpt2 differs by {}'.format(d_bkpt1, d_bkpt2)
            return  True, comment
    else:
        d_bkpt1 = fbkpt1 - gbkpt1
        d_bkpt2 = fbkpt2 - gbkpt2
        if abs(d_bkpt1) < 1000 and abs(d_bkpt2) < 1000:
            comment = ', bkpt1 differs by {}, bkpt2 differs by {}'.format(d_bkpt1, d_bkpt2)
            return True, comment
    return False, ''


def get_near_matches(ddfB, cp, bp, ddfA, to_deleteA, to_deleteB, swap):
    '''
    @param cp - chromosome pair from ddfB
    @param bp - breakpoint pair from ddfB
    @param ddfA - dictionary of fusion events from file A
    '''
    [bp1,bp2] = bp.split(',')
    # allowing a difference of < 10000 bp to be considered a near-match
    bp1,bp2 = bp1[:-4],bp2[:-4]
    if swap:
        matches = filter(lambda x: bp2 == x[:len(bp2)] and bp1 in x,
                                  list(ddfA[cp].keys()))
    else:
        matches = filter(lambda x: bp1 == x[:len(bp1)] and bp2 in x,
                         list(ddfA[cp].keys()))
    matches = list(matches)
    if len(matches) == 0:
        return
    if cp not in to_deleteA:
        to_deleteA[cp] = []
    if cp not in to_deleteB:
        to_deleteB[cp] = []
    for m in matches:
        flagm,comment = check_breakpoint_match(ddfB[cp][bp], ddfA[cp][m], swap)
        if flagm:
            ddfB[cp][bp].mid = ddfA[cp][m].cid
            ddfB[cp][bp].comment += comment

            seqA = ddfA[cp][m].splitr_seq
            seqB = ddfB[cp][bp].splitr_seq
            if swap:
                seqA = get_opp_strand(seqA)
            o, comment = overlap_seq(seqB, seqA)
            ddfB[cp][bp].splitr_seq = o
            ddfB[cp][bp].comment += comment
            to_deleteA[cp].append(m)
            to_deleteB[cp].append(bp)

                     
def near_match(matchfile, nonmatchfile, ddfA, ddfB):
    '''
    @param matchfile - name of file containing matched events
    @param nonmatchfile - name of file containing not matched events
    @param same_ver - bool whether compared files are of the same version
    ''' 
    # Go through non-matches
    print('\nFinding near matches...')
    to_deleteA = {}
    to_deleteB = {}
    cp_match = set(ddfA.keys()).intersection(set(ddfB.keys()))
    for cp in cp_match:
        for bp in ddfB[cp].keys():
            get_near_matches(ddfB, cp, bp, ddfA, to_deleteA, to_deleteB, True)
            get_near_matches(ddfB, cp, bp, ddfA, to_deleteA, to_deleteB, False)
    N = 0
    for cp in to_deleteA:
        N += len(to_deleteA[cp])
    print('\t# near-matches found:\t{}\n'.format(N))

    headers = [ 'cluster_id','gene_chromosome1', 'genomic_break_pos1', 'genomic_strand1', 
                'gene_chromosome2', 'genomic_break_pos2', 'genomic_strand2', 
                'fusion_type','match_id', 'splitr_sequence', 'comment']

    f = open(matchfile, 'a')
    for cp in to_deleteB:
        for bp in to_deleteB[cp]:
          e = ddfB[cp][bp]
          fields = [e.cid, e.chrm1, e.bkpt1, e.strd1, e.chrm2,e.bkpt2,
                    e.strd2, e.ftype, e.mid, e.splitr_seq, e.comment]
          f.write('\t'.join(map(str,fields))+'\n')
    f.close()
    print('Wrote near-matches to {}'.format(matchfile))

    dictionary_cleanup(ddfA, ddfB, to_deleteA, to_deleteB)
          
    #If nonmatchfile already exists, it's overwritten
    nonmatchesA = 0
    nonmatchesB = 0
    with open(nonmatchfile, 'a') as f:
        f.write('\t'.join(headers)+'\n')
        for cp in ddfA:
            for bp in ddfA[cp]:
                e = ddfA[cp][bp]
                fields = [e.cid, e.chrm1, e.bkpt1, e.strd1, e.chrm2,e.bkpt2, 
                          e.strd2, e.ftype, 'File A only', e.splitr_seq, e.comment]
                f.write('\t'.join(map(str,fields))+'\n')
                nonmatchesA += 1
        for cp in ddfB:
            for bp in ddfB[cp]:
                e = ddfB[cp][bp]
                fields = [e.cid, e.chrm1, e.bkpt1, e.strd1, e.chrm2,e.bkpt2, 
                          e.strd2, e.ftype, 'File B only', e.splitr_seq, e.comment]
                f.write('\t'.join(map(str,fields))+'\n')
                nonmatchesB += 1
    print('Wrote non-matched events to {}'.format(nonmatchfile))
    print('\t# of nonmatches from File A:\t{}'.format(nonmatchesA))
    print('\t# of nonmatches from File B:\t{}'.format(nonmatchesB))


def __main__():
    parser = argparse.ArgumentParser(description="Takes two result files from deFuse and outputs the overlay with fusion types called")
    parser.add_argument('fileA', type=str, help='first result file from deFuse')
    parser.add_argument('verA', type=str, help='ensembl version for fileA (ex. GRCh38)')
    parser.add_argument('fileB', type=str, help='second result file from deFuse')
    parser.add_argument('verB', type=str, help='ensembl version for fileB (ex. GRCh37)')
    parser.add_argument('-i', '--id', type=str, help='identifier for results being compared')
    args = parser.parse_args()

    fileA = args.fileA
    filenameA = fileA.split('/')[-1]
    verA = args.verA
    if verA == 'hg38':  verA = 'GRCh38'
    elif verA == 'hg19':  verA = 'GRCh37'
    fileB = args.fileB
    filenameB = fileB.split('/')[-1]
    verB = args.verB
    if verB == 'hg38':  verB = 'GRCh38'
    elif verB == 'hg19':  verB = 'GRCh37'

    if not os.path.isfile(fileA) and not os.path.isfile(fileB):
        if(os.path.isfile(fileA) == False):  
            print('Cannot find file: {}'.format(fileA))
        if(os.path.isfile(fileB) == False):
            print('Cannot find file: {}'.format(fileB))
        sys.exit(-1)

    print('\nInput:')
    print('File A = {}\tReference used = {}'.format(fileA, verA))
    print('File B = {}\tReference used = {}'.format(fileB, verB))
    
    # if the positions are already mapped, read from mapfile           
    mapfileA = '{}.hg38mapped.{}'.format(args.id, filenameA)        
    if verA == 'GRCh37' and os.path.isfile(mapfileA):
        ddfA = read_tsv(mapfileA, 'GRCh38', args.id)
    else:
        ddfA = read_tsv(fileA, verA, args.id)   

    mapfileB = '{}.hg38mapped.{}'.format(args.id, filenameB)
    if verB == 'GRCh37'and os.path.isfile(mapfileB):
        ddfB = read_tsv(mapfileB, 'GRCh38', args.id)
    else:
        ddfB = read_tsv(fileB, verB, args.id)
    # We want to convert the old version to the new version
    ver = min(verA,verB)
    if ver == verA:
        old_file, old_filename, old_ddf = fileA, filenameA, ddfA
        new_file, new_filename, new_ddf = fileB, filenameB, ddfB
    else:
        old_file, old_filename, old_ddf = fileB, filenameB, ddfB
        new_file, new_filename, new_ddf = fileA, filenameA, ddfA

    matchfile = '{}.{}.matched.{}'.format(args.id,  ver, old_filename)
    nonmatchfile = '{}.not_matched.{}'.format(args.id, old_filename)
    with open(nonmatchfile, 'w') as f:
        f.write('## Unmatched events from {} and {}\n'.format(fileA, fileB))
        f.write('## Converted coordinates to GRCh38\n')
        f.write('##File generated on {}\n'.format(time.strftime('%B %d, %Y')))
    overlay(args.id, new_file, old_file, new_ddf, old_ddf, matchfile)
    near_match(matchfile, nonmatchfile, new_ddf, old_ddf)

if __name__ == "__main__":
    __main__()
