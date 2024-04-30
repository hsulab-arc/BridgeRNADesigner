from bridgernadesigner.classes import WTBridgeRNA177nt
from bridgernadesigner import errors

def design_bridge_rna(target, donor):

    target = target.upper()
    donor = donor.upper()

    WTBridgeRNA177nt.check_target_length(target)
    WTBridgeRNA177nt.check_donor_length(donor)
    WTBridgeRNA177nt.check_core_match(target, donor)
    WTBridgeRNA177nt.check_target_is_dna(target)
    WTBridgeRNA177nt.check_donor_is_dna(donor)

    brna = WTBridgeRNA177nt()
    brna.update_target(target)
    brna.update_donor(donor)
    brna.update_hsg()

    return brna

