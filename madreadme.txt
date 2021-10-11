MAD phylogenetic rooting
       - Please cite DOI:10.1038/s41559-017-0193

Files in this folder:

    mad     : Linux executable
    mad.exe : Windows executable
    mad.osx : OS X executable

The executables are standalone wrapers of:
    mad.py  : Python 3 script

Also included:
    mad.m   : Matlab R2015+ function
    mad.R   : R function
    
All 6:
    + Spit a usage message when run with no parameters.
    + Input and output are NEWICK strings.
    - Differences are in sanity checks and options.
      
The sanest and versatilest is the python 
implementation and its executables:
   
Usage:  
        mad filename [-mnsptfgvh]
    
    Where the file contains tree(s) in NEWICK format.
    Trees can be binary or polytomous, unrooted or rooted
    (in which case the pre-exisiting root is summarily ignored).

    Rooted tree(s) will be written to 'filename.rooted',
    rooting statistics to screen.

    By default, the NEWICK output format is the leanest,
    with internal node-labels/bootstrap values stripped,
    and all polytomies arbitrarily resolved. 
    This format should never crash a downstream program.
      
    Flags:

       -m: Multiple trees in input file, output file is verbose and
           contains '>' (info) and '<' (error) lines in addition to 
           NEWICK strings. 
           Default allows exactly one tree, and output is pure NEWICK.
                
           For large datasets, this option speeds up things considerably.
           Also usefull with a single tree to generate the verbose output.
           One drawback, though, is that the additional context lines will
           crash downstream programs if not removed.

       -n: Like -m, but pure newick. Only one rooted tree per input tree,
           errors as a lone ';' in a line.
           
       -s: Statistics - append MAD statistics to output NEWICK strings.

           Again, very useful for automated pipelines, but may crash downstream
           programs if not stripped.
       
       -p: Polytomies - flat polytomies in rooted NEWICK (if present).         
           Default is a binary tree with polytomies arbitrarily resolved
           into zero-length branches.       

       -t: Tiny - retain branch lengths smaller than 10^-6.                
           Default contracts these to 0.0, thereby creating polytomies.  

           The zeroing of extremely short branches is needed to deal with 
           some tree programs that report an epsilon value (usually 10^-8)
           for unresolved or trivial branches, typically coming in from 
           identical sequences. It is also reasonable.

  -f | -g: Figtree - rooted trees formatted for viewing AD scores in figtree.
           Important: after loading in figtree, please check the option 
           'Appearance>Gradient' manually, otherwise the branch colors will
           be misleading. -f reports ancestor deviations (AD) only for nodes,
           while -g reports also within-branch maximizing positions and AD values.
                      
       -v: Version.       
       -h: Help.
                 
Please report bugs to giddy.landan@gmail.com .
            
