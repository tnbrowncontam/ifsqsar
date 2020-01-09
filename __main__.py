"""

"""

if __name__ == "__main__":
    from . import ifsqsar
    from . import models
    import argparse
    import os

    thispath = os.path.dirname(os.path.abspath(__file__))

    # define options for commandline interface
    argparser = argparse.ArgumentParser(prog='ifsqsar',
                                        usage='%(prog)s [options]',
                                        description='IFSQSAR Command Line Interface',
                                        epilog='Invoke ifsqsar with no options start the GUI')
    # output full docs
    argparser.add_argument('-d',
                           '--docs',
                           action='store_true',
                           help='Prints the full documentation to standard output')
    # input must be either from file OR from a list of SMILES
    inputgroup = argparser.add_mutually_exclusive_group()
    inputgroup.add_argument('-i',
                            '--infile',
                            metavar='infile',
                            action='store',
                            type=str,
                            help='Path and name of input file')
    inputgroup.add_argument('-s',
                            '--smiles',
                            metavar='smiles',
                            action='store',
                            type=str,
                            help='Comma-separated list of SMILES')
    # output
    argparser.add_argument('-o',
                           '--outfile',
                           metavar='outfile',
                           action='store',
                           type=str,
                           help='Path and name of output file, if not specified print to standard output')
    # input format specification; only used if input is from file
    # inheaderrows = 1,  # number of header lines
    # inheadertargetrow = 1,  # header row to select from
    # inheadersmiles = 'smiles',  # header value indicating SMILES
    # inseparator = '\t',  # any string
    # inendline = '\n',  # any string
    # output format specification
    # outkeepdata = True,  # also output all of the input file contents
    # outformat = 'rows',  # 'dict', 'columns', 'rows'
    # outheader = True,  # True or False
    # outseparator = '\t',  # any string
    # outendline = '\n',  # any string
    # model selection
    argparser.add_argument('-q',
                           '--qsars',
                           metavar='qsars',
                           action='store',
                           type=str,
                           default=str(models.qsarlist).replace("['", '').replace("', '", ',').replace("']", ''),
                           help='Comma-separated list of qsars to apply, if not specified all are applied. Full list: '
                                 'fhlb, hhlb, hhlt, dsm, tm, E, S, A, B, L, V. See full docs for explanation'
                           )
    # output value selection
    argparser.add_argument('-v',
                           '--values',
                           metavar='values',
                           action='store',
                           type=str,
                           default='insmi,normsmi,sminote,units,qsarpred,UL,error,ULnote',
                           help='Comma-separated list of values to return. Full list: '
                                 'insmi, normsmi, sminote, units, qsarpred, UL, error, ULnote. See full docs for explanation'
                           )
    # parse the options passed then decide actions
    args = argparser.parse_args()
    # no input specified, start the GUI
    if args.docs:
        with open(os.path.join(thispath, 'README.md'), 'r') as readmefile:
            infostring = readmefile.read()
        infostring = infostring.replace('<pre>', '').replace('</pre>', '').replace('\[', '[')
        print(infostring)
    elif args.infile is None and args.smiles is None:
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
                                                      values=args.values.split(','),
                                                      outfilename=args.outfile,
                                                      )
        # print to screen if no outfile
        if args.outfile is None:
            print(result)


