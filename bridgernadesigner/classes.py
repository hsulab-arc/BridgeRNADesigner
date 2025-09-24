from abc import ABC, abstractmethod
from Bio.Seq import reverse_complement
from bridgernadesigner import errors
import re

class BridgeRNAScaffold(ABC):
    """
    Abstract base class for Bridge RNA scaffolds.
    """
    @property
    @abstractmethod
    def TEMPLATE(self):
        pass

    @property
    @abstractmethod
    def GUIDES(self):
        pass

    @property
    @abstractmethod
    def LTG_REGION(self):
        pass

    @property
    @abstractmethod
    def RTG_REGION(self):
        pass

    @property
    @abstractmethod
    def LDG_REGION(self):
        pass

    @property
    @abstractmethod
    def RDG1_REGION(self):
        pass

    @property
    @abstractmethod
    def RDG2_REGION(self):
        pass

    @property
    @abstractmethod
    def TBL_HSG_REGION(self):
        pass

    @property
    @abstractmethod
    def DBL_HSG_REGION(self):
        pass

    @property
    @abstractmethod
    def SPLIT_TBL(self):
        pass

    @property
    @abstractmethod
    def SPLIT_DBL(self):
        pass

    def __init__(self):
        self.bridge_sequence = str(self.TEMPLATE)
        self.target = None
        self.donor = None

    def update_target(self, target):

        self.target = target

        core = target[7:9]
        if core not in ["CT", "GT", "AT", "TT"]:
            errors.NotNTCoreWarning()

        bridgeseq = list(self.bridge_sequence)

        bridgeseq[self.LTG_REGION[0]:self.LTG_REGION[1]] = list(target[:9])
        bridgeseq[self.RTG_REGION[0]:self.RTG_REGION[1]] = list(reverse_complement(target[7:]))

        self.bridge_sequence = "".join(bridgeseq)

    def update_donor(self, donor):
        self.donor = donor

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
    def check_donor_is_dna(donor):
        if re.match("^[ATCG]*$", donor) is None:
            raise errors.DonorNotDNAError()

    @staticmethod
    def check_target_is_dna(target):
        if re.match("^[ATCG]*$", target) is None:
            raise errors.TargetNotDNAError()

    @staticmethod
    def check_core_mismatch(target, donor):
        if target[7:9] != donor[7:9]:
            errors.CoreMismatchWarning()

    @staticmethod
    def check_p6p7_match(target, donor):
        if target[5:7] == donor[5:7]:
            errors.P6P7Warning()
        elif target[6] == donor[6]:
            errors.P7Warning()
        elif target[5] == donor[5]:
            errors.P6Warning()

    def get_p6p7(self):
        return self.target[5:7], self.donor[5:7]

    def format_fasta(self):
        out = ">BridgeRNA_tgt_{}_dnr_{}\n".format(self.target, self.donor)
        out += self.bridge_sequence + "\n"

        # Include split loops if they are defined for this scaffold
        if self.SPLIT_TBL is not None and self.SPLIT_DBL is not None:
            out += ">BridgeRNA_tgt_{}_dnr_{}_split_tbl\n".format(self.target, self.donor)
            out += self.bridge_sequence[self.SPLIT_TBL[0]:self.SPLIT_TBL[1]] + "\n"
            out += ">BridgeRNA_tgt_{}_dnr_{}_split_dbl\n".format(self.target, self.donor)
            out += self.bridge_sequence[self.SPLIT_DBL[0]:self.SPLIT_DBL[1]] + "\n"

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

        target_p6p7 = self.target[5:7]
        donor_p6p7 = self.donor[5:7]
        target_core = self.target[7:9]
        donor_core = self.donor[7:9]

        warning_feat = "#=GF WARNING"
        if target_core != 'CT' and target_core != 'GT' and target_core != 'AT' and target_core != 'TT':
            out += warning_feat + " Core does not follow the expected NT format, recombination may not work.\n"
        if target_core != donor_core:
            out += warning_feat + " The target and donor cores should be identical.\n"
        if target_p6p7 == donor_p6p7:
            out += warning_feat + " Target P6-P7 and Donor P6-P7 match, efficiency may be dramatically reduced.\n"
        elif target_p6p7[1] == donor_p6p7[1]:
            out += warning_feat + " Target P7 and Donor P7 match, efficiency may be significantly reduced.\n"
        elif target_p6p7[0] == donor_p6p7[0]:
            out += warning_feat + " Target P6 and Donor P6 match, efficiency may be reduced.\n"
        out += "//"
        return out

    # Custom print function
    def __str__(self):
        out = "<WTBridgeRNA177nt>\n"
        out += '- Programmed target: {}\n'.format(self.target)
        out += '- Programmed donor: {}\n'.format(self.donor)
        out += '- Bridge RNA sequence: {}'.format(self.bridge_sequence)
        return out


class IS621_WT(BridgeRNAScaffold):
    TEMPLATE = "AGTGCAGAGAAAATCGGCCAGTTTTCTCTGCCTGCAGTCCGCATGCCGTNNNNNNNNNTGGGTTCTAACCTGTNNNNNNNNNTTATGCAGCGGACTGCCTTTCTCCCAAAGTGATAAACCGGNNNNNNNNATGGACCGGTTTTCCCGGTAATCCGTNNTTNNNNNNNTGGTTTCACT"
    GUIDES   = ".................................................LLLLLLLCC...............RRRRRCCHH........................................lllllllc..........................rr..rrrcchh.........."

    LTG_REGION = (49, 58)
    RTG_REGION = (73, 80)
    LDG_REGION = (122, 130)
    RDG1_REGION = (160, 165)
    RDG2_REGION = (156, 158)
    TBL_HSG_REGION = (80, 82)
    DBL_HSG_REGION = (165, 167)
    SPLIT_TBL = None
    SPLIT_DBL = None


class IS622_WT(BridgeRNAScaffold):
    TEMPLATE = "AGTGCAGGGAGAACCGGCCAGTTCTCTCTGCCATGCGGTCCGCATGCCGTNNNNNNNNNCAGGCTAATAACCTGTNNNNNNNNNTTATGCAGCGGACCGCCGTTCTCCACAAGTGACAAACCGGNNNNNNNNATGGACCGGTTTTCCCGGTAATCCGCNNTCNNNNNNNTGGTCTCACT"
    GUIDES   = "..................................................LLLLLLLCC................RRRRRCCHH........................................lllllllc..........................rr..rrrcchh.........."

    LTG_REGION = (50, 59)
    RTG_REGION = (75, 82)
    LDG_REGION = (124, 132)
    RDG1_REGION = (162, 167)
    RDG2_REGION = (158, 160)
    TBL_HSG_REGION = (82, 84)
    DBL_HSG_REGION = (167, 169)
    SPLIT_TBL = (0, 100)
    SPLIT_DBL = (111, 179)


class IS622_enhanced(BridgeRNAScaffold):
    TEMPLATE = "AGTGCAGGGAGAACCGGCCAGTTCTCTCTGCCATGCGGTCCGCATGCCGTNNNNNNNNNTGGGCTAATAACCCGTNNNNNNNNNTGGCAGCGGACCGCGCCGTTCTCCACAAGTGACAAACCGGNNNNNNNNATGGACCGGTTTTCCCGGTAATCCGCNNTCNNNNNNNTGGTCTCACTTGTGGAGAACG"
    GUIDES   = "..................................................LLLLLLLCC................RRRRRCCHH........................................lllllllc..........................rr..rrrcchh....................."

    LTG_REGION = (50, 59)
    RTG_REGION = (75, 82)
    LDG_REGION = (124, 132)
    RDG1_REGION = (162, 167)
    RDG2_REGION = (158, 160)
    TBL_HSG_REGION = (82, 84)
    DBL_HSG_REGION = (167, 169)
    SPLIT_TBL = (0, 98)
    SPLIT_DBL = (100, 190)


SCAFFOLD_NAME_TO_CLASS = {
    'IS621_bRNA': IS621_WT,
    'IS622_bRNA_WT': IS622_WT,
    'IS622_bRNA_enhanced': IS622_enhanced,
}
