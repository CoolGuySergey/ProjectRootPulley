#! /usr/bin/env python3
'''
Niema Moshiri 2016
ECE 286 Homework 2

Generate sequences on a phylogenetic tree using the GTR model
'''
import argparse, dendropy
from common import genString,gtr2matrix,roll
from scipy.linalg import expm
try:
    import Queue as Q  # ver. < 3.0
except ImportError:
    import queue as Q

# evolve given sequence to given time using given rate matrix R
def EvolveSeq(seq, TestTree, R):
    
    Bases = ['A','C','G','T']
    Pmat = expm(R*t)
    P = {'A':{}, 'C':{}, 'G':{}, 'T':{}}
    
    for i in range(4):
        for j in range(4):
            P[Bases[i]][Bases[j]] = Pmat[i,j]
            
    return ''.join([roll(P[c]) for c in seq])

# perform sequence simulation process
def SimulateSeqs(tree, rootseq, pi, gtr):
    
    R = gtr2matrix(gtr,pi)
    root = tree.seed_node
    root.seq = rootseq
    todo = Q.Queue()
    
    for child in root.child_nodes():
        todo.put(child)
        
    while not todo.empty():
        
        node = todo.get()
        node.seq = evolveSeq(node.parent_node.seq, node.edge_length, R)
        
        for child in node.child_nodes():
            todo.put(child)
            
    seqs = {}
    
    for leaf in tree.leaf_node_iter():
        
        seqs[leaf.taxon.label] = leaf.seq
        
    return seqs

# parse arguments
def parseArgs():

    pi, gtr = [[float(i) for i in line.split()] for line in args.gtrparams]
    args.pi = {'A':pi[0], 'C':pi[1], 'G':pi[2], 'T':pi[3]}
    args.gtr = {'CT':gtr[0], 'AT':gtr[1], 'GT':gtr[2], 'AC':gtr[3], 'CG':gtr[4], 'AG':gtr[5]}
    if args.rootseq is None and args.rootseqlen is None:
        print("ERROR: Did not specify a root sequence (-i) nor a root sequence length (-r)")
        exit(-1)
    elif args.rootseq is not None and args.rootseqlen is not None:
        print("User specified both root sequence (-i) and root sequence length (-r)")
        print("Root sequence (-i) will be used, and root sequence length (-r) will be ignored")
        args.rootseqlen = None
    elif args.rootseq is None:
        if args.rootseqlen < 1:
            print("ERROR: Root Sequence Length must be >= 1")
            exit(-1)
        args.rootseq = genString(args.rootseqlen, args.pi)
    else:
        args.rootseq = ''.join([line.strip() for line in args.rootseq][1:])
    return args

# main function
if __name__ == "__main__":
    args = parseArgs()
    tree = dendropy.Tree.get(file=args.tree, schema='newick')
    seqs = simulateSeqs(tree, args.rootseq, args.pi, args.gtr)
    for ID in seqs:
        args.out.write('>' + ID + '\n' + seqs[ID] + '\n')
