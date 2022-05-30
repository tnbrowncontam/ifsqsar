********************************************************************************
**IFSQSAR - A python package for applying QSARs**  
https://github.com/tnbrowncontam/ifsqsar  
Created and maintained by Trevor N. Brown  
Version 1.0.0

IFSQSAR is free to use and redistribute, but is provided "as is" with no implied
warranties or guarantees. The user accepts responsibility for using the software
properly, as outlined in this user guide.

IFSQSAR was coded using only open source python modules, and the QSARs have all
been published in peer-reviewed literature. See below for details on open source
licensing.

**Funding and Acknowledgements**  
Prof. Dr. Frank Wania, Prof. Dr. Kai-Uwe Goss, and Dr. Jon A. Arnot are
acknowledged as coauthors on publications and for supporting this work as 
Ph.D. supervisor, Postdoc supervisor, and P.I., respectively.

Funding sources that have supported relevant publications include:  
- CEFIC-LRI (ECO13, ECO30, B22)
- ACC-LRI
- Alexander von Humboldt Foundation Fellowship 

********************************************************************************
**CONTENTS**
********************************************************************************

 1. IFSQSAR PACKAGE AND SMILES OVERVIEW
 2. COMMAND LINE USAGE
 3. GRAPHICAL USER INTERFACE USAGE
 4. PYTHON USAGE
 5. QSAR DESCRIPTIONS AND INTERPRETATION
 6. INTERPRETING QSAR APPLICABILITY DOMAIN INFORMATION
 7. DEPENDENCIES AND LICENCES
 8. CHANGE LOG
 9. KNOWN BUGS AND PLANNED FEATURES
10. REFERENCES

********************************************************************************
**1. IFSQSAR PACKAGE AND SMILES OVERVIEW**
********************************************************************************

In this section:  
- What is the IFSQSAR package?
- How to use the IFSQSAR package
- Information about SMILES strings
- Custom SMILES specification for mixtures

**What is the IFSQSAR package?**

IFSQSAR is a tool that can be used to apply Quantitative Structure Activity
Relationships (QSARs) developed by Trevor N. Brown, called Iterative Fragment
Selection (IFS) QSARs. There are IFS QSARs implemented to predict chemical
properties relevant for chemical risk assessment, such as biotransformation
half-lives in animals, and Abraham LSER solute descriptors which can be used
to predict environmental partitioning. Select QSPRs from the literature have
also been implemented in IFSQSAR to make predictions for key chemical
properties. IFSQSAR is a provided as a python package, it is open source and
free to use. The source code can be obtained from the GitHub page at:
https://github.com/tnbrowncontam/ifsqsar

**How to use the IFSQSAR package**

To use the IFSQSAR package you must have a python interpreter (3.4 or greater)
installed on your computer. Two dependencies not in the python standard library
are also required: numpy and openbabel. If python is already installed you can
copy the code from GitHub to a folder named "ifsqsar" where python can find it,
typically the Lib\site-packages folder in the python install directory. There
are three ways to use IFSQSAR to make predictions: (1) Import IFSQSAR as a
python package and use it in your own python code. (2) Run IFSQSAR from the
command line to apply the IFS QSARs and handle input and output without writing
any python code. (3) A simple Graphical User Interface (GUI) is provided to
perform single and batch calculations, it can be started from the command line.
The only thing required to make QSAR predictions for a chemical is the structure
as a SMILES string.

**Information about SMILES strings**

Simplified Molecular Input Line Entry System (SMILES) is a way of representing
molecular structures as a string, first described in Weininger 1988 [1]. The
chemistry in IFSQSAR is handled by the open source program Open Babel [2], which
implements the OpenSMILES specification of SMILES [3]. A number of tools exist
for creating or finding SMILES. Many open access database websites such as
PubChem or chemspider.com provide searchable databases with SMILES strings
included. Note that the EPISuite database [4] popular among environmental
chemists contains some SMILES which do not conform to the standard format, these
structures will not be processed correctly if they are entered into IFSQSAR. The
two most common violations are hydrogens outside of brackets, e.g. NH instead of
N[H] or [NH], and counter ions without brackets or disconnected structures,
eg. CC(=O)ONa instead of CC(=O)[O-].[Na+]. In the second example the ONa group
is interpreted as a covalent bond instead of an ion pair, which is incorrect. A
newer alternative also released by the USEPA is the Chemistry Dashboard (found
at https://comptox.epa.gov/dashboard).

IFSQSAR normalizes passed SMILES strings to a canonical format, and attempts
to correct some common issues that would cause the QSAR predictions to fail. The
current QSARs do not make predictions for charged molecules, so IFSQSAR attempts
to neutralize any charged molecules passed to it and removing any counter ions.
Structures are standardized by "inchifying" them [11], which selects a canonical
tautomer, and converts dative bonds (e.g. nitro groups C\[N+](-\[O-])=O to
CN(=O)=O ). If permanent charges are present (e.g. quaternary amines) then no
predictions are made and an error message is printed. Some aromatic structures
may result in an error. In most cases this is because there is some error in
the SMILES and atoms have been flagged as aromatic when they should not have
been, or a valid aromatic structure has been invalidated by other modifications
such as adding counter-ions onto the charged atoms. An error message is printed
in cases of aromatic structure errors and no predictions are made. The modified
SMILES to which the QSARs are applied is provided in the output and any
modifications to structures are noted.

**Custom SMILES specification for mixtures**

Most QSARs implemented in IFSQSAR accept only one SMILES structure as input and
make predictions for one chemical at a time. However, one property (log Ksa) 
makes predictions for solute-solvent pairs (i.e. two chemicals) and future work
will add predictions for the properties of chemicals in mixtures. A custom
SMILES specifications must be used to input chemicals structures for solute-
solvent pairs and mixtures. For solute-solvent pairs the format is:

> {solute}[solute SMILES]{solvent}[solvent SMILES]

Where [solute SMILES] and [solvent SMILES] are the SMILES of the solute and
solvent as they would normally be input for single chemical QSAR calculations.
The solute is specified with the "{solute}" tag and the solvent is specified with
the "{solvent}" tag, all concatenated together with no spaces.

********************************************************************************
**2. COMMAND LINE USAGE**
********************************************************************************

In this section:
- Invoking the command line  
- Command line interface options
- Usage examples 

**Invoking the command line**

To use the ifsqsar as a Command Line Interface (CLI) the package must be run as
a python script from the command line with the syntax:

> python -m ifsqsar [ifsqsar options]

The "-m" option to python runs the package as a script and must be included,
otherwise the models will not load properly and imports will fail. Python must
be installed and environment variables set up properly to call python from a
command prompt. The ifsqsar folder must be placed in a location where python
can find it when invoked, e.g. typically the Lib\site-packages folder in the
python install directory, or the full path to the ifsqsar folder must be given
every time.

**Command line interface options**

Most of the options of the underlying python functions can be accessed from the
CLI. Structures as SMILES must be read from file, or input as a single SMILES or
a comma-separated list of SMILES at the command line. Output can be to a file or
printed to the standard output. Options are available to specify the formatting
of input files and formatting of the output. Some, all or none of the QSARs can
be applied and, in addition to the QSARs predictions, values related to the
uncertainty and details about the SMILES can be output.

Usage:  
> python -m ifsqsar [ifsqsar options]

IFSQSAR Options:  
- -h, --help           : Show help message and exit  
- -d, --docs           : Prints the full documentation to standard output  
- -s, --smiles         : Comma-separated list of SMILES  
- -i, --infile         : Path and name of input file  
- -r, --inheaderrows   : Input file format: number of header rows, default = 1
- -t, --inheadtrgtrow  : Input file format: 1-indexed row number in which to
                         look for SMILES label, default = 1  
- -c, --inheadersmiles : Input file format: label in header that indicates
                         column containing SMILES, default = SMILES
- -p, --inseparator    : Input file format: character used to separate columns,
                         default = \t (tab)
- -e, --inendline      : Input file format: character used to separate rows,
                         default = \n (newline)
- -o, --outfile        : Path and name of output file, if not specified print
                         to standard output
- -u, --outtossdata    : Output format: toss data from input file, default =
                         include data from input file in the output
- -a, --outspprsshead  : Output format: suppress header in output, default =
                         include header in the output
- -f, --outformat      : Output format: chemical structures in rows, or in
                         columns, default = rows
- -m, --outseparator   : Output format: character used to separate columns,
                         default = \\t (tab)
- -n, --outendline     : Output format: character used to separate rows,
                         default = \\n (newline)
- -q, --qsars          : Comma-separated list of QSARs to apply, if option is
                         not invoked all are applied. If option is invoked with
                         no list then none are applied. See section 5 to read
                         about the details of each QSAR. Full list:  
                         fhlb : biotransformation half-life in fish  
                         hhlb : biotransformation half-life in humans  
                         hhlt : total elimination half-life in humans  
                         HLbiodeg : aqueous aerobic biodegradation half-life  
                         dsm : entropy of fusion  
                         tm : melting point (QSPR)  
                         tmpplfer : melting point (PPLFER)  
                         tmconsensus: melting point (mean of QSPR+PPLFER)  
                         tbpplfer : boiling point (PPLFER)  
                         E : Abraham PPLFER solute descriptor  
                         S : Abraham PPLFER solute descriptor  
                         A : Abraham PPLFER solute descriptor  
                         B : Abraham PPLFER solute descriptor  
                         L : Abraham PPLFER solute descriptor  
                         V : Abraham PPLFER solute descriptor  
                         s : Abraham/Goss PPLFER system parameter  
                         a : Abraham/Goss PPLFER system parameter  
                         b : Abraham/Goss PPLFER system parameter  
                         v : Abraham/Goss PPLFER system parameter  
                         l : Abraham/Goss PPLFER system parameter  
                         c : Abraham/Goss PPLFER system parameter  
                         logKow : log Kow (PPLFER)  
                         logKoa : log Koa (PPLFER)  
                         logKaw : log Kaw (PPLFER)  
                         logVPliquid : log Vapor Pressure of liquids (PPLFER)  
                         logSwliquid : log Water Solubil. of liquids (PPLFER)  
                         logSoliquid : log Octanol Solubil. of liquids (PPLFER)  
                         MVliquid : Molar Volume of liquids  
                         MVsolid : Molar Volume of solids  
                         densityliquid : density of liquids  
                         densitysolid : density of solids  
                         MW : Molecular Weight  
                         state : chemical state (gas, liquid or solid)  
                         logKsa : solvent-air partitioning (solute-solvent pair)  
                         pure : all QSARs except logKsa  
                         mixture: only logKsa  
- -v, --values         : Comma-separated list of values to return. If option is
                         not invoked all are included. If option is invoked with
                         no list then none are included. Full list:  
                         insmi    : the SMILES input by the user  
                         normsmi  : the normalized SMILES, see section 1  
                         sminote  : notes about SMILES normalization  
                         qsarpred : the value predicted by the QSAR  
                         endpoint : description of model endpoint  
                         units    : units of QSAR prediction  
                         UL       : Uncertainty Level, see section 6  
                         ULnote   : notes about the Uncertainty Level  
                         error    : estimated prediction uncertainty  
                         citation : literature to cite for QSAR prediction  

**Usage examples**

- command line input:  

  > python -m ifsqsar  
  
  result: start the GUI, see section 3

- command line input:  

  > python -m ifsqsar -s [SMILES] -q pure -f columns 
  
  result: given a [SMILES], print all QSAR results in a column, reproducing the
  behaviour of the single calculation mode in the GUI, see section 3

- command line input:  

  > python -m ifsqsar -i [infile] -p , -o [outfile] -m ,  
  
  result: given a csv input file [infile], output all QSAR results in rows to
  a csv output file [outfile], reproducing the behaviour of the batch
  calculation mode in the GUI, see section 3. Any data in [infile] will also be
  included in [outfile]

- command line input:  

  > python -m ifsqsar -i [infile] -r 3 -t 2 -c [smhe] -o [outfile] -m , -u  
  
  result: given a tab-delimited text file [infile] with three header rows, and
  a custom SMILES column name [smhe] in header row 2, output all QSAR results
  in rows to a csv output file, do not include the other data in [infile].

- command line input:  

  > python -m ifsqsar -s [SMILES] -q -f columns -a -m , -n -v normsmi  
  
  result: given a comma-separated list of [SMILES], print a comma-separated
  list of the normalized SMILES

- command line input:  

  > python -m ifsqsar -s [SMILES] -q tm -a -m , -v qsarpred,error  
  
  result: given a comma-separated list of [SMILES], on one line for each
  structure print the melting point and the estimated prediction uncertainty
  separated by a comma

- command line input:  

  > python -m ifsqsar -s {solute}[SMILES]{solvent}[SMILES] -q mixture  
  
  result: given a SMILES specifying a solute-solvent pair apply all QSARs which
  accept mixtures as inputs (only logKsa in version 1.0.0)

********************************************************************************
**3. GRAPHICAL USER INTERFACE USAGE**
********************************************************************************

In this section:  
- Overview of IFSQSAR single mode interface  
- Overview of IFSQSAR batch mode interface  

**Overview of IFSQSAR single mode interface**

The GUI can be started from the CLI. In single mode QSARs can be applied to a
single structure. The structure must be entered into the interface as a SMILES
string (see below), and then the "Apply IFS QSARs" button must be clicked. Radio
toggle buttons allow the user to choose if QSARs applicable only to pure
chemicals or if those applicable to mixtures are applied. For mixtures the
SMILES specification for mixtures describe in Section 1 must be used. The
results are displayed in the results window as tab-delimited text with headers
and row labels. For further details on interpreting the results see section 2.
INTERPRETING QSAR RESULTS. Note: a new structure can be entered at any time, but
the results will not be displayed until the "Apply IFS QSARs" button is clicked.
Any results in the results window will be lost when a new calculation is
performed. The results cannot be edited in the results window, this is
intentional. If the user switches to batch mode any SMILES structures or results
in the interface will be lost.
<pre>
                            ______________________________________
Enter SMILES string here  ->| Enter a SMILES [_________________] |
Select pure or mixture    ->|                 o Pure  o Mixture  |
                            |                [                 ] |
QSAR results window       ->| Model Results  [                 ] |
                            |                [_________________] |
Switch to batch mode      ->|[ Batch mode ]                      |
Show this info as a popup ->|[    Info    ]   [Apply IFS QSARs]  |
                            |____________________________________|
                             Begin single calculation ^
</pre>
**Overview of IFSQSAR batch mode interface**

In batch mode QSARs can be applied to many structures in series loaded from one
file, and the results are then written to a second file. The file names must be
entered into the interface, and then the "Apply IFS QSARs" button must be
clicked. Radio toggle buttons allow the user to choose if QSARs applicable only
to pure chemicals or if those applicable to mixtures are applied. For mixtures
the SMILES specification for mixtures describe in Section 1 must be used.
For further details on interpreting the results see section 2. INTERPRETING QSAR
RESULTS. Note: selected files can be changed at any time but the QSARs are not
applied and written to file until the "Apply IFS QSARs" button is clicked. Files
must be selected by clicking the relevant select button, they cannot be entered
in the fields on the interface. Warning: selecting an existing file to save the
results to will overwrite the file, its current contents will be lost and
irrecoverable! The output file cannot be the same file as the input file.

The input file must be structured into columns and rows. Each line is a row, and
it is assumed that each row represents the data from a single chemical. SMILES
must be organized into a single column, columns are delimited by a separator
that is specific to the file format (see below). If the file contains a header
(column labels) then the column with the SMILES to be used must be labelled
"SMILES". If the file does not contain a header then it is assumed that the very
first column contains the SMILES. If there are multiple columns containing
SMILES only the first will be used to apply the QSARs. Any additional data in
the file will also be copied to the output file, with the results appended to
the end of each row.

Two file formats are recognized:
1. tab-delimited text files: files with a .txt file extension. Columns are
   separated by tabs.
2. csv files: files with a .csv file extension. Columns are separated by commas.
<pre>
                            __________________________________________
Select input file here    ->|[Select Input File ][_________________] |
Select pure or mixture    ->|                     o Pure  o Mixture  |
Select output file here   ->|[Select Output File][_________________] |
Switch to single mode     ->|[    Single mode   ]                    |
Show this info as a popup ->|[       Info       ] [Apply IFS QSARs]  |
                            |________________________________________|
                                 Begin batch calculations ^
</pre>

********************************************************************************
**4. PYTHON USAGE**
********************************************************************************

In this section:  
- General overview of python usage
- Stability of the API

**General overview of python usage**

Python users should import the ifsqsar package with the syntax:

> from ifsqsar import ifsqsar

The main python interface is the ifsqsar.apply_qsars_to_molecule_list function.
The CLI is built on top of this function and provides access to most of its
functionality. Options are available to specify the format of input files, the
format of output, and what will be output. See the function full documentation.
Users can also access the ifsqsar.apply_qsars_to_molecule function which
applies IFS QSARs to single molecules, but without any file IO options.

Using the ifsqsar.apply_qsars_to_molecule_list function from python provides
access to two additional options not available from the CLI. In addition to
printing the results or writing them to file the results can be returned as
a dict so that the values are directly accessible. In this format, in addition
to the normalized SMILES and model outputs the generated openbabel OBMol object
can be included in the output.

The QSARs are stored in the models subpackage, which can be imported using the
syntax:

> from ifsqsar import models

A list of QSAR objects must be passed to ifsqsar.apply_qsars_to_molecule_list
as the first positional argument. The get_qsar_list function in the models
subpackage creates a list of QSARs from a list of strings of QSAR names. Calling
get_qsar_list without any arguments will return the full list of QSARs,
excluding older version of the QSARs. Older versions of the QSARs can be
accessed by specifying the version number in the optional version argument. As
an example, the overall python code to make a prediction of log Kow for
cyclohexane would look like this:

> from ifsqsar import ifsqsar  
> from ifsqsar import models  
> qsarlist = models.get_qsar_list(['logKow'])  
> results = ifsqsar.apply_qsars_to_molecule(qsarlist, 'C1CCCCC1')

Using the ifsqsar package directly from python will be faster than accessing
the same functionality from the CLI, because every time the CLI is invoked the
models must be loaded, whereas the models only need to be loaded once when the
package is imported as a python module.

**Stability of the API**

The options and naming of the ifsqsar.apply_qsars_to_molecule_list function
will remain stable, as the primary interface to the ifsqsar package. Additional
options may be added based on user feedback, but backwards compatibility will
be maintained across all minor updates, and likely major updates too. The only
API element of the models subpackage that is guaranteed to be stable across
different versions is the models.get_qsar_list function; more options may be
added in the future based on user feedback.

********************************************************************************
**5. QSAR DESCRIPTIONS AND INTERPRETATION**
********************************************************************************

In this section:  
- General description of IFS QSARs  
- fhlb - QSAR for fish biotransformation half-life  
- hhlb - QSAR for human biotransformation half-life  
- hhlt - QSAR for human total elimination half-life  
- E, S, A, B, V, L - QSPRs for Abraham PPLFER solute descriptors
- s, a, b, v, l, c - QSPRs for Abraham/Goss PPLFER system parameters
- dsm - QSPR for Entropy of Melting  
- tm - QSPR for Melting Point  
- General description of other literature QSARs  
- MVliquid, densityliquid - molar volume and density of liquids  
- MVsolid, densitysolid - molar volume and density of solids  
- HLbiodeg - biodegradation half-life based on BIOWIN3/4 from EPISuite  
- General description of Meta QSARs  
- logKow, logKoa, logKaw - commonly used partition coefficients  
- tmpplfer, tbpplfer - melting and boiling point from PPLFER equations  
- tmconsensus - melting point, consensus of the other two QSARs  
- state - physical state of the chemical: gas, liquid or solid  
- logKsa - solvent-air partitioning of a solute-solvent pair  
- logVPliquid, logSwliquid, logSoliquid - pure phase chemical properties  

**General description of IFS QSARs**

IFSQSAR QSARs are developed using a group contribution model, in which specific
substructures of a molecule are counted and their contributions are determined
by coefficients that have been calculated by multiple linear regression. The
specific methodology used is called Iterative Fragment Selection (IFS) and is
described by the author in detail in references [5] and [6]. The QSARS
conform to best practices as outlined by OECD [7], and includes internal and
external model validation. All models have a defined domain of applicability,
which is described in section 6. INTERPRETING QSAR DOMAIN INFORMATION.

**fhlb - QSAR for fish biotransformation half-life**

The model predicts the base-10 logarithm of the whole-body biotransformation
half-life of chemicals in fish (fhlb) measured in days for a reference 10g fish
at a temperature of 288K. In IFSQSAR the model outputs have been converted to
linear (as opposed to log) scale half-lives in hours, and the uncertainty is
provided as a confidence factor (Cf). The QSAR is trained on a dataset derived
from bioaccumulation tests with biotransformation back-calculated with a one-
compartment PK model. The training and validation datasets contain 412 and 207
chemicals, respectively. The validation statistics are r-sq[external] = 0.748
and RMSE of predictions = 0.584, with the data spanning about 5 log units. Full
details are available in reference [5].

**hhlb - QSAR for human biotransformation half-life**

The model predicts the base-10 logarithm of the whole-body biotransformation
half-life of chemicals in humans, (hhlb) measured in hours for a generic 70kg
human. In IFSQSAR the model outputs have been converted to linear (as opposed to
log) scale half-lives in hours, and the uncertainty is provided as a confidence
factor (Cf). As with the fhlb, the QSAR is trained on biotransformation half-
lives backed out of whole body elimination data with a one compartment PK model.
The data is mostly pharmaceuticals measured in humans, with the remainder
environmental contaminants measured in humans. The training and validation
datasets both contain 470 chemicals. The validation statistics are
r-sq[external] = 0.73 and RMSE of predictions = 0.75, with the data spanning
about 7.5 log units. Full details are available in reference [8].

**hhlt - QSAR for human total elimination half-life**

The model predicts the base-10 logarithm of the whole-body total elimination
half-life of chemicals in humans, (hhlt) measured in hours for a generic 70kg
human. In IFSQSAR the model outputs have been converted to linear (as opposed to
log) scale half-lives in hours, and the uncertainty is provided as a confidence
factor (Cf). The QSAR is trained with the same dataset as the hhlb, except as
whole-body elimination half-lives, with some additional data that could not be
transformed into biotransformation half-lives. The training and validation
datasets contain 552 and 553 chemicals, respectively. The validation statistics
are r-sq[external] = 0.72 and RMSE of predictions = 0.70, with the data spanning
about 7.5 log units. Full details are available in reference [8].

**E, S, A, B, V, L - QSPRs for Abraham PPLFER solute descriptors**
A separate IFS QSPR has been developed for each descriptor, excluding V which is
calculated using a simple QSPR from the literature [12]. Preliminary versions of
the QSPRs have been made available on the UFZ LSER database [9], and are
included in IFSQSAR as version 1 of the solute descriptor QSPRs. Results for
version 1 of the QSPRs may differ from the results calculated on the UFZ LSER
database in some cases because of differences in the implementation details.
Version 2 of the PPLFER solute descriptors were trained and externally validated
on a common dataset of 3316 and 1658 chemicals with solute descriptors which are
considered reliable. Two of the QSPRs, A and B, have been bounded to values 
greater than or equal to zero, because values less than zero are not meaningful.
The validation statistics are good, with r-sq[external] > 0.80 in all cases. The
RMSE for the external validation dataset is less than 0.30 for E, S, A, and B;
and is 0.38 for L which has a greater range of values than the other
descriptors. Full details are available in reference [13].

**s, a, b, v, l, c - QSPRs for Abraham/Goss PPLFER system parameters**
The QSPRs are specifically for the system parameters of PPLFER equations for
log Ksa values, ie. partition coefficients of solvent-air systems. The QSPRs
have been developed for each system parameter using a common training and
external validation dataset containing 706 and 353 chemicals, respectively.
Three of the system parameters s, a, and b have been bounded to be greater than
zero. The data are a mix of experimentally calibrated system parameters, and
values which were empirically predicted from reliable solute descriptors using
the methods from reference [14]. Experimental data was available for only 90
solvent-air systems, so adding the empirical data was critical for robust QSPR
development. The validation statistics are adequate with r-sq[external] > 0.70
and < 0.80. The RMSE values are less than 0.30, except for a which is 0.64.
However, the range of the system descriptor values is smaller than for system
parameters than for solute descriptors. Full details are available in reference
[13].

**dsm - QSPR for Entropy of Melting**

The model predicts the entropy of melting (dSm) on a linear scale in units of
kJ/mol. The QSPR is trained on the Jain 2004 dSm data set as described in
reference [10]. The training and validation data sets contain 1056 and 528
chemicals, respectively. The validation statistics are q-sq = 0.792 and std.err.
of prediction = 11.82 kJ/mol, with the data spanning about 175 kJ/mol. Full
details are found in reference [10].

**tm - QSPR for Melting Point**

The model predicts the melting point (Tm) on a linear scale in units of K. The
QSPR is trained on the Bradley data set as described in reference [10]. The
training and validation data sets contain 1922 and 964 chemicals, respectively.
The validation statistics are r-sq[external] = 0.658 and standard error of
prediction = 46.9 K, with the data spanning about 500 K. Full details are found
in reference [10].

**General description of other literature QSARs**

Group contribution QSARs from literature sources have been included in IFSQSAR
to fill in some key chemical properties. SMARTS substructure search strings have
been manually created to match the substructures described in the relevant
publications. Some models use different algorithms for counting the number of
occurrences of fragments in chemical structures. IFS counts all occurrences of a
fragment, including overlapping occurrences, so an atom may be included in
multiple fragments. Other algorithms allow each atom to be a member of only one
fragment, referred to here as exclusive fragments. Another common algorithm is
referred to as structural alerts, in which the fragments only have a count of
zero or one (absence or presence). Each of these algorithms for counting
fragment occurrences has been implemented in IFSQSAR to facilitate implementing
QSARs from the literature.

**MVliquid, densityliquid - molar volume and density of liquids**  
**MVsolid, densitysolid - molar volume and density of solids**  

The QSPR for molar volume (MV), and by combining with MW also the density, is
from Kotomin and Kozlov 2006 [15]. The QSPR is specifically for the MV of solids
with generic correction factors for the MV of liquids applied based on chemical
class. Various fragment types are defined which requires separate model files
for IFS-style fragments, exclusive fragments, and structural alerts. A Meta QSAR
sums the results of these internal models to produce a prediction. See below for
a description of Meta QSARs. No AD has been defined for this QSPR.

**HLbiodeg - biodegradation half-life based on BIOWIN3/4 from EPISuite**  

QSARs that reproduce the BIOWIN3 and BIOWIN4 models in EPI Suite [4] have been
created based on the fragment descriptions and methodology in the original
publication [16]. The dataset was based on ultimate and primary half-life
classifications for 200 chemicals from a survey of experts. Exclusive fragments
and structural alerts are defined in two model files for each QSAR. The IFSQSAR
models reproduce the BIOWIN models with +/- 0.01 for 960 out of a validation
dataset of 1053 diverse chemicals. For chemicals where the predictions do not
exactly match, the average discrepancy is 4.4% of the expected values.
The predictions could not be made to match exactly because of differences in
aromatic ring detection between EPI Suite and OpenBabel, and apparent
inconsistencies in how fragments are detected and counted in EPI Suite.
Biodegradation half-lives are calculated in a Meta QSAR using regressions
between between the log of half-lives and BIOWIN outputs calibrated by Arnot
et al 2005 [17], and taking the average of the log-scale values. See below for
a description of Meta QSARs. Outputs are converted to linear scale hours, and
the uncertainty is provided as confidence factor, same as for the fish and human
elimination half-lives. An estimate of the overall uncertainties was made based
on the variability and uncertainty of the data, QSARs and regressions.

**General description of Meta QSARs**

Meta QSARs take predictions from IFS QSARs, literature QSARs, and other Meta
QSARs and combine them to make a new prediction. The most common application is
to take the solute descriptor QSPR predictions and use them in a PPLFER equation
from the literature to predict partitioning or some other chemical property with
a calibrated equation. Where possible the UL and prediction uncertainties are
calculated by propagation of uncertainty as outlined in reference [13]. In some
cases this was insufficient and UL and prediction uncertainty are calculated
using different methods.

**logKow, logKoa, logKaw - commonly used partition coefficients**

Base-10 logarithm of the partition coefficients (partition ratios) of octanol-
water, octanol-air and air-water systems. These are "dry" systems, i.e. in the
case of octanol-water the two phases are assumed to be pure, with no mutual
solubility as would be found if they were in direct contact during experimental
determination. The system parameters are taken from reference [14] and the UL
and prediction uncertainty are calculated by propagation of uncertainty.

**tmpplfer, tbpplfer - melting and boiling point from PPLFER equations**

Temperatures of melting and boiling as predicted by the PPLFER equations
calibrated in reference [13]. Propagation of uncertainty was found to be
insufficient, underestimating the observed discrepancies in the external
validation dataset. UL is instead set using the leverage vs. training datasets
of solute descriptors and cut-off values which are calculated the same way as
is described in section 6. Prediction uncertainties were calculated for each UL
with the standard error of prediction for the external validation datasets.

**tmconsensus - melting point, consensus of the other two QSARs**

Consensus prediction of the temperature of melting, calculated by taking the
mean of the direct QSPR prediction (tm) and the PPLFER prediction (tmpplfer).
The UL and prediction uncertainties are calculated by propagation of uncertainty
as described in reference [13].

**state - physical state of the chemical: gas, liquid or solid**

This Meta QSAR attempts to determine the physical state of a chemical at room
temperature, i.e. gas, liquid or solid. No numerical outputs are returned, a
text description of the state is return in the ULnote. Two complementary methods
are used; first, the dataset of liquid solvents used to train and validate the
QSPRs for system parameters s, a, b, v, l, c is used to calculate a leverage
with the predicted solute descriptors. If the leverage indicates the chemical is
in the AD then this is considered evidence that it is also a liquid. The
predicted tm (tmconsensus) and tb (tbpplfer) are also checked to see if the
chemical is a liquid, i.e. tm < room temperature < tb. If the two methods agree
then a "likely" adjective is assigned to the assigned state, otherwise a "maybe"
adjective is assigned.

**logKsa - solvent-air partitioning of a solute-solvent pair**

The model predicts the base-10 logarithm of the partition coefficient (partition
ratio) of a solute partitioning in a solvent-air system. The solute and the
solvent must be both be input using the SMILES specification for mixtures
defined in section 1. The model applies Henry's Law assumptions, i.e. that the
solute is at infinite dilution and the mole fraction of the solvent in the
liquid phase is equal to unity. The entire workflow for assigning system
parameters and AD information for the solvent is implemented as described in
reference [13]. The methods were developed specifically for liquid solvents, so
applying them to gases or solids is much more uncertain. The state Meta QSAR is
called as a solvent dependency to determine physical date of the solvent, and
the prediction uncertainty is multiplied by an additional uncertainty scaling
factor if the solvent is not a liquid.

**logVPliquid, logSwliquid, logSoliquid - pure phase chemical properties**

The logVPliquid MetaQSAR is identical to the logKsa Meta QSAR, except that in
this case the solute is also the solvent, so the SMILES specification for
mixtures is not required as input. The output is the base-10 logarithm of the
vapor pressure. The molar volume Meta QSAR (MVliquid) is used to do a unit
conversion to VP in Pa. The VP is assumed to be for either the liquid state or
supercooled liquid state. The solubility of the chemical in water (logSwliquid)
and the solubility of the chemical in octanol (logSoliquid) are calculated from
the VP by applying thermodynamic cycles using the Meta QSARs for logKaw and
logKoa.

********************************************************************************
**6. INTERPRETING QSAR APPLICABILITY DOMAIN INFORMATION**
********************************************************************************

In this section:  
- What is the applicability domain?  
- Applicability domain of for the IFSQSAR QSARs  
- Interpreting IFSQSAR applicability domain outputs  
- Applicability domain of Meta QSARs  

**What is the applicability domain?**

Having a defined applicability domain (AD) is one of the five OECD guidelines
for developing QSARs. The OECD guidance document [7] gives more information on
AD in general, and references [5] and [6] describe the definition of the AD for
IFS QSARs specifically. The basic problem is that the chemical space to which a
QSAR might be applied for predictions is extremely large, but datasets available
 for training QSARs can only ever hope to cover a fraction of this chemical
space. Therefore, applying QSARs means extrapolating beyond the training set in
many cases. The AD tries to quantify how far from the training dataset the model
can make reliable extrapolations.

There is no universally accepted standard method for defining the AD, but most
methods follow a few basic principles. Robust models trained on good data have a
larger AD than uncertain models trained on noisy data. Interpolating between
data points is more reliable than extrapolating outside of the data range.
Making predictions for chemicals similar to those in the training data is more
robust than making predictions for chemicals that are less similar. In addition
the limitations of the model used is a factor. The IFS QSARs are based on the 2D
structure of chemicals, i.e. only atom connectivity is considered, the three
dimensional shape of molecules is not included in the model. Therefore, in cases
where the 3D shape of the molecule might be important, e.g. due to steric
effects, the models are expected to give less accurate results.

**Applicability domain of the IFSQSAR QSARs**

Two different methods of defining the AD are applied to the IFSQSAR QSARs in
conjunction:
1. Chemical Similarity Score (CSS), as defined by the author in [5] and [6].
2. Leverage, as discussed by the author in [6].

The CSS simultaneously measures the similarity of a chemical to those in the
training dataset and how well the model fits those chemicals in the training
dataset. It accounts for both the robustness of the model and the similarity
to chemicals in the training dataset, but does not discriminate strongly between
interpolating and extrapolating from the dataset. Leverage is a measure of
structural similarity to chemicals in the training dataset, it discriminates
strongly against extrapolation outside of the training dataset but does not
account for how well the model fits the data.

Based on cut-off values a chemical is classified on a scale of 0-2 for CSS, and
a scale of 0-3 for leverage, with the higher of these two taken as the
"uncertainty level" (UL) of the prediction. From the validation dataset the
standard error of prediction of each UL is estimated from the prediction
residuals of chemicals with that UL.

**Interpreting IFSQSAR applicability domain outputs**

For each property prediction three domain fields are included:
- error: the estimated standard error of prediction (uncertainty) as described
  above
- UL: the uncertainty level, higher numbers indicate less reliable predictions
- note: explanation for UL
    
Error (Uncertainty):  
For the fhlb, hhlb, hhlt, and HLbiodeg QSARs the original estimated standard
errors of prediction on the log scale have been converted to factor errors, also
called dispersion. The error is multiplied by 1.96 and then transformed to
linear scale. This gives the factor X in the equation:  
<pre>
     P ( prediction/X < actual < prediction*X ) = 0.95
</pre>  
where the P is the probablility of the actual value being in the interval
defined by prediction/X to prediction*X, i.e. the 95% prediction interval of the
QSAR. This is analogous to the confidence factor Cf.

For several of the PPLFER QSPRs (E, S, A, B, s, a, b) the intercept of the model
has been set to zero, and many chemicals do not have any fragment overlap with
the training dataset (warning level 4 below). These are in fact good
predictions, and the low estimated standard prediction errors of prediction
reflect this. These LSER descriptors are essentially normalized to aliphatic
hydrocarbons and so predictions for these types of chemicals will typically have
a warning level 4.

Uncertainty Level (UL) and Notes:  
<pre>
Uncertain.   Note           Explanation
 Level
   0                        Molecule is considered to be within the domain by
                            both leverage and CSS.

   1    low similarity      CSS flags molecule as a "borderline" case more prone
                            to prediction errors than UL 0 molecules

   1    high leverage       Leverage flags the molecule as a "borderline" case
                            more prone to prediction errors than UL 0 molecules.

   2    out of domain       CSS flags the molecule as out of the domain of
                            applicability. The prediction might still be ok, but
                            should be used with caution because the standard
                            error of prediction is typically high.

   2    structural outlier  Leverage flags the molecule as a structural outlier,
                            meaning the prediction involves some extrapolation.
                            The prediction might still be ok, but should be used
                            with caution because standard error of prediction is
                            typically high.

   3    leverage > 1        Leverage is greater than one, indicating that the
                            molecule contains substructures for which neither
                            interpolation or extrapolation is sufficient.
                            These predictions are likely to be erroneous, with
                            very high standard errors of prediction.

   4    no fragment         None of the substructures in the model are in the
        overlap with        molecule. This may be fine if some such chemicals
        training dataset    are in the training dataset, meaning the prediction
                            is simply the intercept. If so the uncertainty is
                            calculated from any such chemicals in the validation
                            dataset. Otherwise the uncertainty is set to the
                            same value as the full validation set.

   5    atom violation,     The simplest domain check is to examine the atoms in
        various messages    the molecule to be predicted and see if it contains
                            unique atom types that are not present in the
                            training dataset. The features of the atom that are
                            unique will be listed, e.g. element type, number of
                            hydrogens attached, etc. The uncertainty is
                            calculated from any such chemicals in the validation
                            dataset. Otherwise, uncertainty is set to the same
                            value as UL 3.

   6    prediction less/    For some datasets predictions outside of the
        greater than        experimental range do not make sense. For example,
        largest/smallest    elimination of chemicals from fish or humans can
        value in training   only proceed so rapidly, so predicting faster
        dataset             elimination is not reasonable. In these "bounded"
                            QSARs predictions outside of the experimental range
                            are set to the value of the closest boundary. The
                            uncertainty is the same as the normal UL 0-3 the
                            molecule belongs to.
</pre>
********************************************************************************
**7. DEPENDENCIES AND LICENCES**
********************************************************************************

IFSQSAR was written with the open source resources listed below. Because two of
these resources use the GNU General Public License the source code for IFSQSAR
is available on GitHub.

Python: GUI and general coding  
- Python Programming Language, http://www.python.org/  
- License: BeOpen.com GPL-compatible license.  

NumPy: Mathematics  
- NumPy, http://www.numpy.org/  
- License: Developer specific license. http://www.numpy.org/license.html  

Open Babel: Chemical structure handling  
- The Open Babel Package, http://openbabel.org  
- License: GNU GPL. http://www.gnu.org/licenses/gpl.html  

********************************************************************************
**8. CHANGE LOG**
********************************************************************************

Version 1.0.0
- First official release
- python package, CLI and GUI functionality implemented
- Chemical properties implemented are:
    - fhlb
    - hhlb, hhlt
    - HLbiodeg
    - dsm
    - tm, tmpplfer, tmconsensus
    - PPLFER solute descriptors: E, S, A, B, L, V
    - PPLFER system parameters (solvent-air system): s, a, b, v, l, c
    - logKow, logKoa, logKaw
    - logVPliquid, logSwliquid, logSoliquid
    - MVliquid, MVsolid, densityliquid, densitysolid, MW
    - state
    - Solvent-air partitioning (solute and solvent input): logKsa

Version 1.0.1
- Added datasets folder to package for disseminating select datasets
- Updated citations for Brown 2022

********************************************************************************
**9. KNOWN BUGS AND PLANNED FEATURES**
********************************************************************************

- SMILES specification for mixtures still to be fully implemented to include
  components, i.e. chemicals which are treated as both solutes and solvents,
  along with input of mole fractions for each component
- Planned properties still to be implemented:
    - mixture vapor-liquid equilibrium (VLE)
    - Human skin permeability (Kp) for pure chemicals and mixtures
    - logVP from PPLFERs for liquids and solids, possible consensus liquid VP
    - logKow, logKoa for "wet", i.e. water saturated, octanol
    - IFS QSAR for HLbiodeg

********************************************************************************
**10. REFERENCES**
********************************************************************************

 1. Weininger D., 1988, J. Chem. Inf. Model. 28 (1): 31-6.  
    DOI:10.1021/ci00057a005

 2. O'Boyle N.M., Banck M., James C.A., Morley C., Vandermeersch T., 
    Hutchison G.R. Open Babel: An open chemical toolbox. J. Cheminf. 2011, 3,  
    33\. DOI: 10.1186/1758-2946-3-33

 3. James, C.A. OpenSMILES Specification. 2016.  
    http://opensmiles.org/

 4. E.P.A., U.S., Estimation Programs Interface (EPI) Suite for Microsoft(R)
    Windows, Ver. 4.1. 2011, U. S. Environmental Protection Agency: Washington,
    D.C.

 5. Brown T.N., Arnot J.A., Wania F. Iterative Fragment Selection: A group
    contribution approach to predicting fish biotransformation half-lives.
    Environmental Science & Technology, 2012, 46, 8253-8260.  
    DOI: 10.1021/es301182a

 6. Brown T.N. Predicting hexadecane-air equilibrium partition coefficients (L)
    with a group contribution approach constructed from high-quality data. SAR
    and QSAR in Environmental Research, 2013.  
    DOI: 10.1080/1062936x.2013.841286

 7. OECD Guidance Document on the Validation of (Quantitative) Structure-
    Activity Relationship [(Q)SAR] Models.; Organisation for Economic
    Cooperation and Development, Environment Directorate: Paris, 2007.

 8. Arnot J., Brown T.N., Wania F.; Estimating screening-level organic chemical
    half-lives in humans. Environmental Science & Technology, 2013.  
    DOI: 10.1021/es4029414

 9. Ulrich, N., Endo, S., Brown, T.N., Watanabe, N., Bronner, G., Abraham, M.H.,
    Goss, K.-U., UFZ-LSER database v 3.2.1 [Internet], Leipzig, Germany,
    Helmholtz Centre for Environmental Research-UFZ. 2017. 
    http://www.ufz.de/lserd

10. Brown T.N., Armitage J.M., Arnot J.A.; Application of an iterative fragment
    selection (IFS) method to estimate entropies of fusion and melting points of
    organic chemicals. Molecular Informatics, 2019. DOI: 10.1002/minf.201800160  
    
11. O’Boyle, N. M., Towards a Universal SMILES representation - A standard
    method to generate canonical SMILES based on the InChI. Journal of
    Cheminformatics 2012, 4 (1), 22. DOI: 10.1186/1758-2946-4-22

12. Mcgowan, J. C., The Estimation of Solubility Parameters and Related
    Properties of Liquids. J Chem Tech Biot A 1984, 34 (1), 38-42.
    DOI: 10.1002/jctb.5040340107

13. Brown, T.N., QSPRs for Predicting Equilibrium Partitioning in Solvent–Air
    Systems from the Chemical Structures of Solutes and Solvents. Journal of
    Solution Chemistry, 2022. DOI: 10.1007/s10953-022-01162-2

14. Brown, T. N., Empirical Regressions between System Parameters and Solute
    Descriptors of Polyparameter Linear Free Energy Relationships (PPLFERs) for
    Predicting Solvent-Air Partitioning. Fluid Phase Equilibria 2021, 113035.
    DOI: 10.1016/j.fluid.2021.113035

15. Kotomin, A. A.; Kozlov, A. S., Calculation of densities of organic compounds
    from contributions of molecular fragments. Russ J Appl Chem 2006, 79 (6),
    957-966. DOI: 10.1134/S1070427206060176

16. Boethling, R. S.;  Howard, P. H.;  Meylan, W.;  Stiteler, W.;  Beauman, J.;
    Tirado, N., Group contribution method for predicting probability and rate of
    aerobic biodegradation. Environmental science & technology 1994, 28 (3),
    459-65. DOI: 10.1021/es00052a018

17. Arnot, J. A.;  Gouin, T.; Mackay, D. Development and Application of Models
    of Chemical Fate in Canada - Practical methods for estimating environmental
    biodegradation rates; 2005.

********************************************************************************