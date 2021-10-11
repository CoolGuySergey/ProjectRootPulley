
# Description: GTR

#========================================================================


import random

# Define weights
ProbDict = {'A' : 0.2,
            'T' : 0.2,
            'C' : 0.3,
            'G' : 0.3,
}

def MultiRoll(ProbDict, NumberOfRolls):
    
    """
    Generate a seq of length int(NumberOfRolls)
    """
    
    CandidateBases = sorted(ProbDict.keys())
    Biases = [ProbDict[Base] for Base in CandidateBases]
    
    # Describe a biased distribution NOT uniform distribution.
    # i.e. Choose a number in [0,1] but divide 4 parts
    # depending on the bias wanted
    # e.g. in the test case with 40% AT and 60% GC (sorted A, C, G, T):
    # [0, 0.2],[0.2, 0.5], [0.5, 0.7], [0.7, 1]
    cdf = [Biases[0]]
    for i in range(1,len(Biases)):
        cdf.append(Biases[i] + cdf[i-1])

    AllDraws = []

    for _ in range(NumberOfRolls):
        RandomNum = random.uniform(0, 1)
        
        # In which bin does RandomDraw fall?
        i = 0
        while cdf[i] < RandomNum:
            i += 1
        SelectedBase = CandidateBases[i]
        AllDraws.append(SelectedBase)
        
    return ''.join(AllDraws)


#========================================================================


import numpy as np

# e.g.
# pi = ProbDict, stationary vector
ProbDict = {'A' : 0.2,
            'T' : 0.2,
            'C' : 0.3,
            'G' : 0.3,
}

# GTR parameter dictionary
GTRParamDict = {'AC' : 0.25,
                'AG' : 0.25,
                'AT' : 0.25,
                'CG' : 0.25,
                'CT' : 0.25,
                'GT' : 0.25,
}

# numpy rate matrix (indexed: 0 = A, 1 = C, 2 = G, 3 = T)
def gtr2matrix(ProbDict, GTRParamDict):

    """
    Generate a seq of length int(NumberOfRolls)
    """
    
    Bases = ['A','C','G','T']
    R = [[None,None,None,None] for _ in range(4)]
    
    # row A: outgoing & incoming A's
    R[0][1] = GTRParamDict['AC']*ProbDict['C']
    R[0][2] = GTRParamDict['AG']*ProbDict['G']
    R[0][3] = GTRParamDict['AT']*ProbDict['T']
    R[0][0] = -1 * (R[0][1] + R[0][2] + R[0][3]) # main diagonal, pos AA
    
    # row C: outgoing & incoming C's
    R[1][0] = GTRParamDict['AC']*ProbDict['A']
    R[1][2] = GTRParamDict['CG']*ProbDict['G']
    R[1][3] = GTRParamDict['CT']*ProbDict['T']
    R[1][1] = -1 * (R[1][0] + R[1][2] + R[1][3]) # main diagonal, pos CC
    
    # row G: outgoing & incoming G's
    R[2][0] = GTRParamDict['AG']*ProbDict['A']
    R[2][1] = GTRParamDict['CG']*ProbDict['C']
    R[2][3] = GTRParamDict['GT']*ProbDict['T']
    R[2][2] = -1 * (R[2][0] + R[2][1] + R[2][3]) # main diagonal, pos GG
    
    # row T: outgoing & incoming T's
    R[3][0] = GTRParamDict['AT']*ProbDict['A']
    R[3][1] = GTRParamDict['CT']*ProbDict['C']
    R[3][2] = GTRParamDict['GT']*ProbDict['G']
    R[3][3] = -1 * (R[3][0] + R[3][1] + R[3][2]) # main diagonal, pos TT
    
    # normalize such that sum(vi*pi) = 1
    norm = -1 * sum([R[i][i] * ProbDict[Bases[i]] for i in range(4)])
    #print(f"norm is {norm}.")
    
    for i in range(4):
        for j in range(4):
            R[i][j] = R[i][j] / norm
            
    return R, norm

R, norm = gtr2matrix(ProbDict, GTRParamDict)


#========================================================================

R = [[-0.2, 0.075, 0.075, 0.05],
     [0.05, -0.175, 0.075, 0.05],
     [0.05, 0.075, -0.175, 0.05],
     [0.05, 0.075, 0.075, -0.2],
]

norm = 0.185

ProbDict = {'A' : 0.2,
            'T' : 0.2,
            'C' : 0.3,
            'G' : 0.3,
}

# Silly way to recover norm above
def matrix2gtr(R, norm,  ProbDict):

    """
    Generate a seq of length int(NumberOfRolls)
    """
    
    Bases = ['A','C','G','T']

    for i in range(4):
        for j in range(4):
            R[i][j] = R[i][j] * norm
            
    FinalParamDict = {
        'AC':float(R[0][1])/ProbDict['C'],
        'AG':float(R[0][2])/ProbDict['G'],
        'AT':float(R[0][3])/ProbDict['T'],
        'CG':float(R[1][2])/ProbDict['G'],
        'CT':float(R[1][3])/ProbDict['T'],
        'GT':float(R[2][3])/ProbDict['T'],
    }

    return FinalParamDict

GTRParamDict2 = matrix2gtr(R, norm, ProbDict)


#========================================================================


from dendropy import Tree
from decimal import *
from math import log
from scipy.linalg import expm

StarterStr = "(Alice:0.1,Ben:0.2,(Charles:0.3,Diana:0.4)Eva:0.5)Fred;"
TestTree = Tree.get(data = StarterStr,
             schema = "newick")

SeqDict = {'Alice' : 'GACCCGGGTT',
           'Ben' : 'AACCCAGGTT',
           'Charles' : 'GGGTTAACCT',
           'Diana' : 'GGGTCAACCC',
           }

ProbDict = {'A' : 0.2,
            'T' : 0.2,
            'C' : 0.3,
            'G' : 0.3,
}

R = [[-0.2, 0.075, 0.075, 0.05],
     [0.05, -0.175, 0.075, 0.05],
     [0.05, 0.075, -0.175, 0.05],
     [0.05, 0.075, 0.075, -0.2],
]

# compute maximum likelihood of tree
# given tree, sequence data, stationary comp, and GTR parameters
def L(TestTree, SeqDict, ProbDict, R):

    """
    Generate a seq of length int(NumberOfRolls)
    """
    
    Bases = ['A','C','G','T']
    SeqLength = len(list(SeqDict.values())[0])
  
    # Felsenstein pruning algorithm
    for Node in TestTree.postorder_node_iter():
        
        if Node.is_leaf():
            
            LeafLabel = str(Node.taxon).strip("'").strip('"')
            LeafSeq = SeqDict[LeafLabel]
            print(Node.L)
            
            Node.L = [{
                'A' : Decimal(0),
                'C' : Decimal(0),
                'G' : Decimal(0),
                'T' : Decimal(0)} for i in range(SeqLength)]
            
            for i in range(SeqLength):
                Node.L[i][LeafSeq[i]] = Decimal(1)
    
        else:
            
            Node.L =[{
                'A' : Decimal(1),
                'C' : Decimal(1),
                'G' : Decimal(1),
                'T' : Decimal(1)} for i in range(SeqLength)]
            
            for Child in Node.child_node_iter():
                
                Pmat = expm(R * Child.edge_length)
                P = {'A':{}, 'C':{}, 'G':{}, 'T':{}}
                
                for i in range(4):
                    for j in range(4):
                        P[Bases[i]][Bases[j]] = Decimal(Pmat[i, j])
                        
                for i in range(SeqLength):
                    for j in Bases:
                        Node.L[i][j] *= sum([P[j][x] * Child.L[i][b] for b in Bases])

    ans = sum([sum([Decimal(ProbDict[b])*tree.seed_node.L[i][b] for b in Bases]).ln() for i in range(SeqLength)])
    
    return float(ans)
