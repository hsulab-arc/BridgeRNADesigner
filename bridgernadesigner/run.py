from bridgernadesigner.classes import SCAFFOLD_NAME_TO_CLASS


def design_bridge_rna(target, donor, scaffold):
    target = target.upper()
    donor = donor.upper()
    brna_scaffold = SCAFFOLD_NAME_TO_CLASS[scaffold]

    # Hard checks, will raise errors if not met
    brna_scaffold.check_target_length(target)
    brna_scaffold.check_donor_length(donor)
    brna_scaffold.check_target_is_dna(target)
    brna_scaffold.check_donor_is_dna(donor)

    # Warning if not met
    brna_scaffold.check_core_mismatch(target, donor)
    brna_scaffold.check_p6p7_match(target, donor)

    brna = brna_scaffold()
    brna.update_target(target)
    brna.update_donor(donor)
    brna.update_hsg()

    return brna
