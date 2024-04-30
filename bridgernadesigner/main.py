import click
import sys
from bridgernadesigner import errors, run
from bridgernadesigner.classes import WTBridgeRNA177nt

@click.group(invoke_without_command=True)
@click.option('--target', '-t', required=True, help="14nt target sequence")
@click.option('--donor', '-d', required=True, help="14nt donor sequence")
@click.option('--output-format', '-of', default='stockholm', type=click.Choice(['stockholm', 'fasta']))
@click.option('--include-annealing-oligos/--no-include-annealing-oligos', default=True,
              help="Include annealing oligos in the output. Only available for FASTA output format. default=True")
@click.option('--annealing-oligos-lh-overhang', default="TAGC", help="5' overhang for annealing oligos. default=TAGC")
@click.option('--annealing-oligos-rh-overhang', default="GGCC", help="3' overhang for annealing oligos. default=GGCC")
@click.pass_context
def cli(ctx, target, donor, output_format, include_annealing_oligos,
        annealing_oligos_lh_overhang, annealing_oligos_rh_overhang):
    if ctx.invoked_subcommand is None:
        ctx.invoke(default_command, target=target, donor=donor, output_format=output_format,
                   include_annealing_oligos=include_annealing_oligos,
                   annealing_oligos_lh_overhang=annealing_oligos_lh_overhang,
                   annealing_oligos_rh_overhang=annealing_oligos_rh_overhang)

@cli.command()
@click.pass_context
def default_command(ctx, target, donor, output_format,
                    include_annealing_oligos, annealing_oligos_lh_overhang, annealing_oligos_rh_overhang):

    target = target.upper()
    donor = donor.upper()

    try:
        WTBridgeRNA177nt.check_target_length(target)
    except errors.TargetLengthError as e:
        print("ERROR:", e)
        print("Exiting...")
        sys.exit()

    try:
        WTBridgeRNA177nt.check_donor_length(donor)
    except errors.DonorLengthError as e:
        print("ERROR:", e)
        print("Exiting...")
        sys.exit()

    try:
        WTBridgeRNA177nt.check_core_match(target, donor)
    except errors.CoreMismatchError as e:
        print("ERROR:", e)
        print("Exiting...")
        sys.exit()

    try:
        WTBridgeRNA177nt.check_target_is_dna(target)
    except errors.TargetNotDNAError as e:
        print("ERROR:", e)
        print("Exiting...")
        sys.exit()

    try:
        WTBridgeRNA177nt.check_donor_is_dna(donor)
    except errors.DonorNotDNAError as e:
        print("ERROR:", e)
        print("Exiting...")
        sys.exit()

    brna = run.design_bridge_rna(target, donor)
    if output_format == "fasta":
        print(brna.format_fasta(include_annealing_oligos=include_annealing_oligos,
                                lh_overhang=annealing_oligos_lh_overhang,
                                rh_overhang=annealing_oligos_rh_overhang))
    elif output_format == "stockholm":
        print(brna.format_stockholm())


if __name__ == '__main__':
    cli()
