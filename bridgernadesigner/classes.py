from Bio.Seq import reverse_complement
from bridgernadesigner import errors
import re

class WTBridgeRNA177nt:

    TEMPLATE = "AGTGCAGAGAAAATCGGCCAGTTTTCTCTGCCTGCAGTCCGCATGCCGTNNNNNNNNNTGGGTTCTAACCTGTNNNNNNNNNTTATGCAGCGGACTGCCTTTCTCCCAAAGTGATAAACCGGNNNNNNNNATGGACCGGTTTTCCCGGTAATCCGTNNTTNNNNNNNTGGTTTCACT"
    DOT_STRUC_P1 = "...(((((((((((......))))))))))).((((((((((.(((............(((((....)))))..............))))).)))))))).........((((....(((.............((((((.....)))))...)...............)))..))))"
    DOT_STRUC_P2 = "((.(((((((((((......)))))))))))))(((((((((.(((.............<(((.<>.)))>...............))))))))))))...........((((...<(((..........<.(((((((.....)))))...))....>.........)))>.))))"
    GUIDES = ".................................................LLLLLLLCC...............RRRRRCCHH........................................lllllllc..........................rr..rrrcchh.........."

    LTG_REGION = (49, 58)
    RTG_REGION = (73, 80)
    LDG_REGION = (122, 130)
    RDG1_REGION = (160, 165)
    RDG2_REGION = (156, 158)
    TBL_HSG_REGION = (80, 82)
    DBL_HSG_REGION = (165, 167)

    def __init__(self):
        self.bridge_sequence = str(self.TEMPLATE)
        self.target = None
        self.donor = None

    def update_target(self, target):

        self.target = target

        core = target[7:9]
        if core not in ["CT", "GT"]:
            errors.NotCTorGTCoreWarning()
        if core not in ["CT", "GT", "AT", "TT"]:
            errors.NotNTCoreWarning()

        bridgeseq = list(self.bridge_sequence)

        bridgeseq[self.LTG_REGION[0]:self.LTG_REGION[1]] = list(target[:9])
        bridgeseq[self.RTG_REGION[0]:self.RTG_REGION[1]] = list(reverse_complement(target[7:]))

        self.bridge_sequence = "".join(bridgeseq)

    def update_donor(self, donor):
        self.donor = donor

        if self.has_donor_p7C():
            errors.DonorP7CWarning()

        bridgeseq = list(self.bridge_sequence)
        bridgeseq[self.LDG_REGION[0]:self.LDG_REGION[1]] = list(donor[:8])
        bridgeseq[self.RDG1_REGION[0]:self.RDG1_REGION[1]] = list(reverse_complement(donor[7:-2]))
        bridgeseq[self.RDG2_REGION[0]:self.RDG2_REGION[1]] = list(reverse_complement(donor[-2:]))

        self.bridge_sequence = "".join(bridgeseq)

    def update_hsg(self):
        bridgeseq = list(self.bridge_sequence)

        target_p6p7, donor_p6p7 = self.get_p6p7()

        if target_p6p7 != donor_p6p7:
            bridgeseq[self.TBL_HSG_REGION[0]:self.TBL_HSG_REGION[1]] = list(reverse_complement(donor_p6p7))
            bridgeseq[self.DBL_HSG_REGION[0]:self.DBL_HSG_REGION[1]] = list(reverse_complement(target_p6p7))
        else:
            errors.MatchingP6P7Warning()
            bridgeseq[self.TBL_HSG_REGION[0]:self.TBL_HSG_REGION[1]] = list(reverse_complement(donor_p6p7))
            bridgeseq[self.DBL_HSG_REGION[0]:self.DBL_HSG_REGION[1]] = list(target_p6p7)

        self.bridge_sequence = "".join(bridgeseq)

    @staticmethod
    def check_target_length(target):
        if len(target) != 14:
            raise errors.TargetLengthError()

    @staticmethod
    def check_donor_length(donor):
        if len(donor) != 14:
            raise errors.DonorLengthError()

    @staticmethod
    def check_core_match(target, donor):
        if target[7:9] != donor[7:9]:
            raise errors.CoreMismatchError()

    @staticmethod
    def check_donor_is_dna(donor):
        if re.match("^[ATCG]*$", donor) is None:
            raise errors.DonorNotDNAError()

    # @staticmethod
    def check_target_is_dna(target):
        if re.match("^[ATCG]*$", target) is None:
            raise errors.TargetNotDNAError()

    def get_p6p7(self):
        return self.target[5:7], self.donor[5:7]

    def has_donor_p7C(self):
        if self.donor[6] == "C":
            return True
        else:
            return False

    def has_matching_p6p7(self):
        if self.target[5:7] == self.donor[5:7]:
            return True
        else:
            return False

    def annealing_oligos(self, lh_overhang="TAGC", rh_overhang="GGCC"):
        fiveprime_stem_loop_top = lh_overhang+self.bridge_sequence[:45]
        fiveprime_stem_loop_btm = reverse_complement(self.bridge_sequence[:49])
        tbl_top = self.bridge_sequence[45:105]
        tbl_btm = reverse_complement(self.bridge_sequence[49:109])
        dbl_top = self.bridge_sequence[105:177]
        dbl_btm = rh_overhang+reverse_complement(self.bridge_sequence[109:177])

        return {
            'fiveprime_stem_loop_top': fiveprime_stem_loop_top,
            'fiveprime_stem_loop_btm': fiveprime_stem_loop_btm,
            'tbl_top': tbl_top, 'tbl_btm': tbl_btm,
            'dbl_top': dbl_top, 'dbl_btm': dbl_btm
        }

    def format_fasta(self, include_annealing_oligos=False,
                     lh_overhang="TAGC", rh_overhang="GGCC"):
        out = ">BridgeRNA_tgt_{}_dnr_{}\n".format(self.target, self.donor)
        out += self.bridge_sequence
        if include_annealing_oligos:
            oligos = self.annealing_oligos(lh_overhang, rh_overhang)
            out += '\n>BridgeRNA_tgt_{}_dnr_{}_oligo_anneal_fiveprime_stem_loop_top\n{}\n'.format(
                self.target, self.donor, oligos['fiveprime_stem_loop_top']
            )
            out += '>BridgeRNA_tgt_{}_dnr_{}_oligo_anneal_fiveprime_stem_loop_btm\n{}\n'.format(
                self.target, self.donor, oligos['fiveprime_stem_loop_btm']
            )
            out += '>BridgeRNA_tgt_{}_dnr_{}_oligo_anneal_tbl_top\n{}\n'.format(
                self.target, self.donor, oligos['tbl_top']
            )
            out += '>BridgeRNA_tgt_{}_dnr_{}_oligo_anneal_tbl_btm\n{}\n'.format(
                self.target, self.donor, oligos['tbl_btm']
            )
            out += '>BridgeRNA_tgt_{}_dnr_{}_oligo_anneal_dbl_top\n{}\n'.format(
                self.target, self.donor, oligos['dbl_top']
            )
            out += '>BridgeRNA_tgt_{}_dnr_{}_oligo_anneal_dbl_btm\n{}'.format(
                self.target, self.donor, oligos['dbl_btm']
            )

        return out

    def format_stockholm(self, whitespaces=5):
        seqname = "BridgeRNA_tgt_{}_dnr_{}".format(self.target, self.donor)
        leader_cols = len(seqname)+whitespaces
        out = "# STOCKHOLM 1.0\n"
        out += seqname+" "*(leader_cols-len(seqname))+self.bridge_sequence+"\n"
        template_feat = "#=GC bRNA_template"
        out += template_feat + " "*((leader_cols-len(template_feat))) + self.TEMPLATE + "\n"
        guide_feat = "#=GC guides"
        out += guide_feat + " "*((leader_cols-len(guide_feat))) + self.GUIDES + "\n"
        ss_feat = "#=GC SS"
        out += ss_feat + " " * ((leader_cols - len(ss_feat))) + self.DOT_STRUC_P2 + "\n"

        donor_p7 = self.donor[6]
        target_p6p7 = self.target[5:7]
        donor_p6p7 = self.donor[5:7]
        core = self.target[7:9]

        warning_feat = "#=GF WARNING"
        if donor_p7 == 'C':
            out += warning_feat + " Donor P7 is C, this was found to be very inefficient in a screen.\n"
        if target_p6p7 == donor_p6p7:
            out += warning_feat + " Target P6-P7 and Donor P6-P7 match, efficiency is unclear.\n"
        if core != 'CT' and core != 'GT':
            out += warning_feat + " Core is not CT or GT, efficiency is unclear.\n"
        if core != 'CT' and core != 'GT' and core != 'AT' and core != 'TT':
            out += warning_feat + " Core does not follow the expected NT format, likely inefficient.\n"
        out += "//"
        return out

    # Custom print function
    def __str__(self):
        out = "<WTBridgeRNA177nt>\n"
        out += '- Programmed target: {}\n'.format(self.target)
        out += '- Programmed donor: {}\n'.format(self.donor)
        out += '- Bridge RNA sequence: {}'.format(self.bridge_sequence)
        return out