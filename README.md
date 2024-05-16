# BridgeRNADesigner
A command line interface for designing Bridge RNAs for Bridge Editing. Accompanies the paper
["Bridge RNAs direct modular and programmable recombination of target and 
donor DNA"](https://www.biorxiv.org/content/10.1101/2024.01.24.577089v1) by Durrant & Perry et al. 2024.

## Dependencies
* Tested on python 3.8
* Tested on Linux and Mac

## Installation
```bash
pip install bridgernadesigner
```

## Usage
Command to design a Bridge RNA for a given target and donor sequence.
```bash
brna-design --target ATCGGGCCTACGCA --donor ACAGTATCTTGTAT
```

Example output:
```bash
# STOCKHOLM 1.0
BridgeRNA_tgt_ATCGGGCCTACGCA_dnr_ACAGTATCTTGTAT     AGTGCAGAGAAAATCGGCCAGTTTTCTCTGCCTGCAGTCCGCATGCCGTATCGGGCCTTGGGTTCTAACCTGTTGCGTAGATTTATGCAGCGGACTGCCTTTCTCCCAAAGTGATAAACCGGACAGTATCATGGACCGGTTTTCCCGGTAATCCGTATTTACAAGGCTGGTTTCACT
#=GC bRNA_template                                  AGTGCAGAGAAAATCGGCCAGTTTTCTCTGCCTGCAGTCCGCATGCCGTNNNNNNNNNTGGGTTCTAACCTGTNNNNNNNNNTTATGCAGCGGACTGCCTTTCTCCCAAAGTGATAAACCGGNNNNNNNNATGGACCGGTTTTCCCGGTAATCCGTNNTTNNNNNNNTGGTTTCACT
#=GC guides                                         .................................................LLLLLLLCC...............RRRRRCCHH........................................lllllllc..........................rr..rrrcchh..........
#=GC SS                                             ((.(((((((((((......)))))))))))))(((((((((.(((.............<(((.<>.)))>...............))))))))))))...........((((...<(((..........<.(((((((.....)))))...))....>.........)))>.))))
//
```
This should run very quickly on a Mac or Linux machine. The code can also automatically generate annealing oligos in 
FASTA format:
```bash
brna-design --target ATCGGGCCTACGCA --donor ACAGTATCTTGTAT -of fasta --include-annealing-oligos
```
Output:
```bash
>BridgeRNA_tgt_ATCGGGCCTACGCA_dnr_ACAGTATCTTGTAT
AGTGCAGAGAAAATCGGCCAGTTTTCTCTGCCTGCAGTCCGCATGCCGTATCGGGCCTTGGGTTCTAACCTGTTGCGTAGATTTATGCAGCGGACTGCCTTTCTCCCAAAGTGATAAACCGGACAGTATCATGGACCGGTTTTCCCGGTAATCCGTATTTACAAGGCTGGTTTCACT
>BridgeRNA_tgt_ATCGGGCCTACGCA_dnr_ACAGTATCTTGTAT_oligo_anneal_fiveprime_stem_loop_top
TAGCAGTGCAGAGAAAATCGGCCAGTTTTCTCTGCCTGCAGTCCGCATG
>BridgeRNA_tgt_ATCGGGCCTACGCA_dnr_ACAGTATCTTGTAT_oligo_anneal_fiveprime_stem_loop_btm
ACGGCATGCGGACTGCAGGCAGAGAAAACTGGCCGATTTTCTCTGCACT
>BridgeRNA_tgt_ATCGGGCCTACGCA_dnr_ACAGTATCTTGTAT_oligo_anneal_tbl_top
CCGTATCGGGCCTTGGGTTCTAACCTGTTGCGTAGATTTATGCAGCGGACTGCCTTTCTC
>BridgeRNA_tgt_ATCGGGCCTACGCA_dnr_ACAGTATCTTGTAT_oligo_anneal_tbl_btm
TTGGGAGAAAGGCAGTCCGCTGCATAAATCTACGCAACAGGTTAGAACCCAAGGCCCGAT
>BridgeRNA_tgt_ATCGGGCCTACGCA_dnr_ACAGTATCTTGTAT_oligo_anneal_dbl_top
CCAAAGTGATAAACCGGACAGTATCATGGACCGGTTTTCCCGGTAATCCGTATTTACAAGGCTGGTTTCACT
>BridgeRNA_tgt_ATCGGGCCTACGCA_dnr_ACAGTATCTTGTAT_oligo_anneal_dbl_btm
GGCCAGTGAAACCAGCCTTGTAAATACGGATTACCGGGAAAACCGGTCCATGATACTGTCCGGTTTATCACT
```

In the `#=GC guides` track, the `L` character indicates the programmable nucleotides of the LTG, the `R` character
indicates the programmable nucleotides of the RTG, the `C` character indicates the guide nucleotides of the 
target core, the `H` character indicates the position of the programmable handshake guide nucleotides, the 
`l` character indicates the programmable nucleotides of the LDG, the `r` character
indicates the programmable nucleotides of the RDG, the `c` character indicates the guide nucleotides of the 
donor core, and the `h` character indicates the position of the programmable handshake guide nucleotides.

# Citation
Please cite [our paper](https://www.biorxiv.org/content/10.1101/2024.01.24.577089v1) if you use any aspect of this 
code in your work.