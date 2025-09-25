# BridgeRNADesigner
A command line interface for designing bridge RNAs for bridge editing. Accompanies the papers
["Bridge RNAs direct modular and programmable recombination of target and 
donor DNA"](https://doi.org/10.1038/s41586-024-07552-4), Durrant & Perry et al. (2024), and ["Megabase-scale human genome rearrangement with programmable bridge recombinases"](https://doi.org/10.1126/science.adz0276), Perry et al. (2025).

## Dependencies
* Tested on python 3.8
* Tested on Linux and Mac

## Installation
```bash
pip install bridgernadesigner
```

## Usage
Command to design a bridge RNA for a given target and donor sequence, and a bridge RNA scaffold.

```bash
brna-design --target ATCAGGCCTACGTC --donor ACAGTATCTTGTAT --scaffold ISCro4_enhanced
```

Example output:
```bash
>BridgeRNA_tgt_ATCAGGCCTACGTC_dnr_ACAGTATCTTGTAT_scaffold_ISCro4_enhanced
AGTGCAGGGAGAACCGGCCAGTTCTCTCTGCCATGCGGTCCGCATGCCGTATCAGGCCTTGGGCTAATAACCCGTGACGTAGATTGGCAGCGGACCGCGCCGTTCTCCACAAGTGACAAACCGGACAGTATCATGGACCGGTTTTCCCGGTAATCCGCATTCACAAGGCTGGTCTCACTTGTGGAGAACG
>BridgeRNA_tgt_ATCAGGCCTACGTC_dnr_ACAGTATCTTGTAT_scaffold_ISCro4_enhanced_split_tbl
AGTGCAGGGAGAACCGGCCAGTTCTCTCTGCCATGCGGTCCGCATGCCGTATCAGGCCTTGGGCTAATAACCCGTGACGTAGATTGGCAGCGGACCGC
>BridgeRNA_tgt_ATCAGGCCTACGTC_dnr_ACAGTATCTTGTAT_scaffold_ISCro4_enhanced_split_dbl
CGTTCTCCACAAGTGACAAACCGGACAGTATCATGGACCGGTTTTCCCGGTAATCCGCATTCACAAGGCTGGTCTCACTTGTGGAGAACG
```

This should run very quickly on a Mac or Linux machine. There is also an option to output the sequence in Stockholm format with multiple tracks.

```bash
brna-design --target ATCAGGCCTACGTC --donor ACAGTATCTTGTAT --scaffold ISCro4_enhanced -of stockholm
```

```bash
# STOCKHOLM 1.0
BridgeRNA_tgt_ATCAGGCCTACGTC_dnr_ACAGTATCTTGTAT_scaffold_ISCro4_enhanced     AGTGCAGGGAGAACCGGCCAGTTCTCTCTGCCATGCGGTCCGCATGCCGTATCAGGCCTTGGGCTAATAACCCGTGACGTAGATTGGCAGCGGACCGCGCCGTTCTCCACAAGTGACAAACCGGACAGTATCATGGACCGGTTTTCCCGGTAATCCGCATTCACAAGGCTGGTCTCACTTGTGGAGAACG
#=GC bRNA_template                                                           AGTGCAGGGAGAACCGGCCAGTTCTCTCTGCCATGCGGTCCGCATGCCGTNNNNNNNNNTGGGCTAATAACCCGTNNNNNNNNNTGGCAGCGGACCGCGCCGTTCTCCACAAGTGACAAACCGGNNNNNNNNATGGACCGGTTTTCCCGGTAATCCGCNNTCNNNNNNNTGGTCTCACTTGTGGAGAACG
#=GC guides                                                                  ..................................................LLLLLLLCC................RRRRRCCHH........................................lllllllc..........................rr..rrrcchh.....................
//
```

In the `#=GC guides` track, the `L` character indicates the programmable nucleotides of the LTG, the `R` character
indicates the programmable nucleotides of the RTG, the `C` character indicates the guide nucleotides of the 
target core, the `H` character indicates the position of the programmable handshake guide nucleotides, the 
`l` character indicates the programmable nucleotides of the LDG, the `r` character
indicates the programmable nucleotides of the RDG, the `c` character indicates the guide nucleotides of the 
donor core, and the `h` character indicates the position of the programmable handshake guide nucleotides.

# Citation
Please cite the linked publications above if you use any aspect of this code in your work.
