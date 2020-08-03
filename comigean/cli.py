#! /usr/bin/env python
"""
Usage:
    comigean [--version] [--help] <command> [<args>...]

General options:
    -h, --help       Help screen
    -v, --version    Get version number

The most commonly used comigean commands are:
    install-db       Installs required databases
    get-genomes      Download genomes (& call genes) from taxa of interest
    genome-stats     Get genome stats
    find-markers     Find marker genes, individually align & concatenate

See 'comigean <command> --help' for more information on a specific command.
"""

import os
from docopt import docopt
from subprocess import call

import comigean.dwnld_genomes
import comigean.get_taxa_db
import comigean.find_markers
import comigean.profile_genomes


if __name__ == '__main__':
    args = docopt(__doc__,
                  version='comigean Beta',
                  options_first=True)
    argv = [args['<command>']] + args['<args>']

    if args['<command>'] == 'get-genomes':
        import comigean.commands.getgenomes_command
        args = docopt(comigean.commands.getgenomes_command.__doc__,
                      argv=argv)

        if args["--USER_GENOMES"]:
            if not os.path.exists(args["--USER_GENOMES"]):
                raise Exception("Path to --USER_GENOMES does not exist.")

        if not os.path.exists(args["<OUTDIR>"]) and not args["--count"]:
            os.makedirs(args["<OUTDIR>"])

        if not os.path.exists(args["<REF_DIR>"]):
            raise Exception(f'{args["<REF_DIR>"]} does not contain the required databases. Did you run the install-db command?')

        if args["--REF_TAXA"] is None and args["--OUTGROUP_TAXA"] is None and args["--USER_GENOMES"] is None:
            raise Exception("Must provide at least one of the following: --USER_GENOMES, --REF_TAXA, or --OUTGROUP_TAXA.")

        comigean.dwnld_genomes.get_genomes(args["--REF_TAXA"],
                                              args["--OUTGROUP_TAXA"],
                                              args["--USER_GENOMES"],
                                              args["<OUTDIR>"],
                                              args["<REF_DIR>"],
                                              args["--count"],
                                              args["--refseq_category"],
                                              args["--assembly_level"],
                                              args["--code"])

    elif args['<command>'] == 'genome-stats':
        import comigean.commands.profile_command
        args = docopt(comigean.commands.profile_command.__doc__,
                      argv=argv)

        comigean.profile_genomes.profile(args["<REF_DIR>"], args["<OUTDIR>"])

    elif args['<command>'] == 'install-db':
        import comigean.commands.install_db_command
        args = docopt(comigean.commands.install_db_command.__doc__, argv=argv)

        if not os.path.exists(args["<DIR>"]):
            os.makedirs(args["<DIR>"])
        elif len(os.listdir(args["<DIR>"])) == 0:
            pass
        else:
            raise Exception(f'{args["<DIR>"]} already exists')

        comigean.get_taxa_db.get_db(args["<DIR>"])

    elif args['<command>'] == 'find-markers':
        import comigean.commands.hmm_command
        args = docopt(comigean.commands.hmm_command.__doc__, argv=argv)

        if not os.path.exists(args["<PROTEOME_DIR>"]):
            raise Exception(f'{args["<PROTEOME_DIR>"]} does not exist.')

        if not os.path.exists(args["<REF_DIR>"]):
            raise Exception(f'{args["<REF_DIR>"]} does not exist.')

        comigean.find_markers.phylogenomics(args["<PROTEOME_DIR>"], args["<REF_DIR>"])

    else:
        print(f"{args['<command>']} is not a comigean command. See 'comigean --help'.")
