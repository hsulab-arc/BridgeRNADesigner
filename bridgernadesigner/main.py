import click
import sys
from bridgernadesigner import errors, run
from bridgernadesigner.classes import SCAFFOLD_NAME_TO_CLASS


@click.group(invoke_without_command=True)
@click.option('--target', '-t', required=True, help="14nt target sequence")
@click.option('--donor', '-d', required=True, help="14nt donor sequence")
@click.option('--scaffold', '-s', required=True, type=click.Choice(['IS621_bRNA', 'IS622_bRNA_WT', 'IS622_bRNA_enhanced']), help="Template bridge RNA scaffold")
@click.option('--output-format', '-of', default='fasta', type=click.Choice(['fasta', 'stockholm']))
@click.pass_context
def cli(ctx, target, donor, scaffold, output_format):
    if ctx.invoked_subcommand is None:
        ctx.invoke(default_command, target=target, donor=donor, scaffold=scaffold, output_format=output_format)

@cli.command()
@click.pass_context
def default_command(ctx, target, donor, scaffold, output_format):

    target = target.upper()
    donor = donor.upper()
    brna_scaffold = SCAFFOLD_NAME_TO_CLASS[scaffold]

    try:
        brna_scaffold.check_target_length(target)
    except errors.TargetLengthError as e:
        print("ERROR:", e)
        print("Exiting...")
        sys.exit()

    try:
        brna_scaffold.check_donor_length(donor)
    except errors.DonorLengthError as e:
        print("ERROR:", e)
        print("Exiting...")
        sys.exit()

    try:
        brna_scaffold.check_target_is_dna(target)
    except errors.TargetNotDNAError as e:
        print("ERROR:", e)
        print("Exiting...")
        sys.exit()

    try:
        brna_scaffold.check_donor_is_dna(donor)
    except errors.DonorNotDNAError as e:
        print("ERROR:", e)
        print("Exiting...")
        sys.exit()

    brna = run.design_bridge_rna(target, donor, scaffold)
    if output_format == "fasta":
        print(brna.format_fasta())
    elif output_format == "stockholm":
        print(brna.format_stockholm())


if __name__ == '__main__':
    cli()
