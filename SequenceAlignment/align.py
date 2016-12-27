#! /usr/bin/python
import argparse
import numpy as np
import matplotlib.pyplot as plt

def d(s, t, matchscore=2, mismatchscore=-1, gapscore=-2):
    '''
    used in place of a substitution matrix
    s,t     STR     sequences to align
    '''
    if s == '-' or t == '-':
        return gapscore
    elif s == t:
        return matchscore
    else:
        return mismatchscore
    
def smithwaterman(s, t, M, B, local=False):
    '''
    s,t     STR     two strings to align
    M       ARR     cost matrix
    B       ARR     backtracking matrix
    local   BOOL    if True, finds local alignments (default global)
    '''
    s = '-'+s
    t = '-'+t
    n = len(s)
    m = len(t)
    for i in range(n):
        for j in range(m):
            deletion,insertion,pairing = -100,-100,-100
            direction = ""
            if i == 0 and j == 0:
                B[i,j] = []
                continue
            else:
                if i > 0 and j > 0:
                    pairing = M[i-1,j-1] + d(s[i], t[j])
                    direction = "\\"
                if i > 0:
                    deletion = M[i-1,j] + d(s[i], "-")
                    direction = "-"
                if j > 0:
                    insertion = M[i,j-1] +  d("-", t[j])
                    direction = "|"
                if local:
                    score = max(pairing, deletion, insertion,0)
                else:
                    score = max(pairing, deletion, insertion)
            M[i,j] = score

            direction = []
            if pairing == score:
                direction.append("\\")
            if deletion == score:
                direction.append("-")
            if insertion == score:
                direction.append("|")
            B[i,j] = direction

def backtrack(s, t, i, j, s_align, t_align, B):
    # if you've finished backtracking (reached beginning of alignment)
    if i == 0 and j == 0:
        print(s_align)
        print(t_align+"\n")
    # if you still need to backtrack
    for k in B[i,j]:
        if k == "-":    # a deletion
            backtrack(s, t, i-1, j, s[i-1]+s_align,'-'+t_align,B)
        if k == "|":    # an insertion
            backtrack(s, t, i, j-1, '-'+s_align, t[j-1]+t_align,B)
        if k == "\\":    # a pairing
            backtrack(s, t, i-1, j-1, s[i-1]+s_align, t[j-1]+t_align,B)

def backtrack_local(s, t, i, j, s_align, t_align, M, B):
    # if you've finished backtracking (reached beginning of alignment)
    if i < 0 or j < 0 or M[i,j] == 0:
        print(s_align)
        print(t_align+"\n")
    # if you still need to backtrack
    else:
        for k in B[i,j]:
            if k == "-":    # a deletion
                backtrack_local(s, t, i-1, j, s[i-1]+s_align,'-'+t_align,M,B)
            if k == "|":    # an insertion
                backtrack_local(s, t, i, j-1, '-'+s_align, t[j-1]+t_align,M,B)
            if k == "\\":    # a pairing
                backtrack_local(s, t, i-1, j-1, s[i-1]+s_align, t[j-1]+t_align,M,B)

def print_matrix(s,t,M,B):
    '''
    combines matrices M and B to display the DP table
    '''
    n = len(s)+1
    m = len(t)+1
    d = 0.7

    plt.axis([-1,n,-0.5,m])
    ax = plt.axes()
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    
    
    # header sequence t
    for j in range(m-1):
        ax.text(-0.5, j, t[m-2-j].upper(), horizontalalignment='center',verticalalignment='center')
    # header sequence s
    for i in range(n-1):
        ax.text(i+1, m-0.5, s[i].upper(), horizontalalignment='center',verticalalignment='center')
    # arrows
    gray = "#D3D3D3"
    for i,v in enumerate(B):
        for j,u in enumerate(v):
            for t in u:
                if t == "|":
                    ax.arrow(i, m-j, 0, -d, head_width=0.1, head_length=0.1, fc=gray, ec=gray)
                if t == "-":
                    ax.arrow(i-1, m-j-1, d, 0, head_width=0.1, head_length=0.1, fc=gray, ec=gray)
                if t == "\\":
                    ax.arrow(i-1, m-j, d, -d, head_width=0.1, head_length=0.1, fc=gray, ec=gray)
    # scores
    for i,v in enumerate(M):
        for j,u in enumerate(v):
            ax.text(i, m-j-1, u, horizontalalignment='center',verticalalignment='center')
    plt.show()
    
def __main__():
    parser = argparse.ArgumentParser(description="Given two sequences, find the optimal global or local alignments.\n Returns the alignment scores, optimal alignments, and dynamic programming table.\n(Scoring: match = 2, mismatch = -1, indel = -2).")
    parser.add_argument('s', type=str, help='A string')
    parser.add_argument('t', type=str, help='Another string')
    parser.add_argument('-l', dest='local', action='store_true', help='if specified, find optimal local alignment (default=global)')

    args = parser.parse_args()
    s = args.s
    t = args.t
    print("\nGiven sequences \"{}\" and \"{}\"".format(s,t))
    
    M = np.zeros([len(s)+1,len(t)+1], dtype=int)
    B = np.empty([len(s)+1,len(t)+1], dtype=list)
    smithwaterman(s, t, M, B, args.local)
    if args.local:
        print("Local alignment returns a max score of {}\n".format(M.max()))
    else:
        print("Global alignment returns a score of {}\n".format(M[-1,-1]))
        
    print("The optimal alignments are:")
    if args.local:
        maxscore = M.max()
        for i,row in enumerate(M):
            for j, col in enumerate(row):
                if col == maxscore:
                    partial_s = s[:i]
                    partial_t = t[:j]
                    backtrack_local(partial_s, partial_t, i, j, '', '', M, B)
    else:
        backtrack(s,t,len(s), len(t), '', '', B)
    print_matrix(s,t,M,B)
    
if __name__ == '__main__':
    __main__()
