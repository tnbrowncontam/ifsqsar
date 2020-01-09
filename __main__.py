"""

"""

if __name__ == "__main__":
    from . import ifsqsar
    from . import models
    import argparse

    # define options for commandline interface
    argparser = argparse.ArgumentParser(prog='ifsqsar',
                                        usage='%(prog)s [options]',
                                        description='IFSQSAR Command Line Interface',
                                        epilog='Invoke ifsqsar with no options start the GUI')
    # input must be either from file OR from a list of SMILES
    inputgroup = argparser.add_mutually_exclusive_group()
    inputgroup.add_argument('-i',
                            '--infile',
                            metavar='infile',
                            action='store',
                            type=str,
                            help='path and name of input file')
    inputgroup.add_argument('-s',
                            '--smiles',
                            metavar='smiles',
                            action='store',
                            type=str,
                            help='comma-separated list of SMILES')
    # output
    argparser.add_argument('-o',
                           '--outfile',
                           metavar='outfile',
                           action='store',
                           type=str,
                           help='path and name of output file, if not specified prints to standard output')
    # input format specification; only used if input is from file
    # output format specification
    # model selection
    argparser.add_argument('-q',
                           '--qsars',
                           metavar='qsars',
                           action='store',
                           type=str,
                           default=str(models.qsarlist).replace("['", '').replace("', '", ',').replace("']", ''),
                           help='list of qsars to apply, if not specified all are applied')
    # parse the options passed then decide actions
    args = argparser.parse_args()
    # no input specified, start the GUI
    if args.infile is None and args.smiles is None:
        ifsqsar.main()
    else:
        # load QSARs
        qsarmodels = []
        for q in args.qsars.split(','):
            qsarmodels.append(getattr(models, q))
        # infile name has been passed, apply models to file
        if args.infile is not None:
            filename = args.infile
            smiles = None
        # SMILES list has been passed, apply models to list
        elif args.smiles is not None:
            filename = None
            smiles = args.smiles.split(',')
        # apply models
        result = ifsqsar.apply_qsars_to_molecule_list(qsarmodels,
                                                      infilename=filename,
                                                      smileslist=smiles,
                                                      outfilename=args.outfile,
                                                      )
        # print to screen if no outfile
        if args.outfile is None:
            print(result)


