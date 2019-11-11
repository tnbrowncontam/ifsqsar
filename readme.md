********************************************************************************
**IFSAPP - A simple program for applying QSARs**\
Coded by Trevor N. Brown\
Version 0.0, Nov. 2019

IFSAPP is free to use and redistribute, but is provided "as is" with no implied
warranties or guarantees. The user accepts responsibility for using the software
properly, as outlined in this user guide.

IFSAPP was coded using only open source python modules, and the QSARs have all
been published in peer-reviewed literature.
********************************************************************************

CONTENTS

1. SOFTWARE OVERVIEW AND USAGE
2. QSAR DESCRIPTIONS AND INTERPRETATION
3. INTERPRETING QSAR DOMAIN INFORMATION
4. SOFTWARE CODING EXPLANATION
5. CHANGE LOG
6. REFERENCES

********************************************************************************
1. SOFTWARE OVERVIEW AND USAGE
********************************************************************************

In this section:  
-Overview of IFSAPP single mode interface  
-Overview of IFSAPP batch mode interface  
-Information about SMILES strings

**Overview of IFSAPP single mode interface**

In single mode QSARs can be applied to a single structure. The structure must be
entered into the interface as a SMILES string (see below), and then the "Apply
IFS QSARs" button must be clicked. The results are displayed in the results
window as tab-delimited text with headers and row labels. For further details on
interpreting the results see section 2. INTERPRETING QSAR RESULTS. Note: a new
structure can be entered at any time, but the results will not be displayed
until the "Apply IFS QSARs" button is clicked. Any results in the results
window will be lost when a new calculation is performed. The results cannot be
edited in the results window, this is intentional. If the user switches to batch
mode any SMILES structures or results in the interface will be lost.
<pre>
                            ______________________________________
Enter SMILES string here  ->| Enter a SMILES [_________________] |
                            |                [                 ] |
QSAR results window       ->| Model Results  [                 ] |
                            |                [_________________] |
Switch to batch mode      ->|[ Batch mode ]                      |
Show this info as a popup ->|[    Info    ]   [Apply IFS QSARs]  |
                            |____________________________________|
                             Begin single calculation ^
</pre>
**Overview of IFSAPP batch mode interface**

In batch mode QSARs can be applied to many structures in series loaded from one
file, and the results are then written to second file. The structure must be
entered into the interface as a SMILES string (see below), and then the "Apply
IFS QSARs" button must be clicked. For further details on interpreting the
results see section 2. INTERPRETING QSAR RESULTS. Note: selected files can be
changed at any time but the QSARs are not applied and written to file until the
"Apply IFS QSARs" button is clicked. Files must be selected by clicking the
relevant select button, they cannot be entered in the fields on the interface.
Warning: selecting an existing file to save the results to will overwrite the
file, its current contents will be lost and irrecoverable! The output file
cannot be the same file as the input file.

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
Select output file here   ->|[Select Output File][_________________] |
Switch to single mode     ->|[    Single mode   ]                    |
Show this info as a popup ->|[       Info       ] [Apply IFS QSARs]  |
                            |________________________________________|
                                 Begin batch calculations ^
</pre>
**Information about SMILES strings**

Simplified Molecular Input Line Entry System (SMILES) is a way of representing
molecular structures as a string, first described in Weininger 1988 [1]. The
chemistry in IFSAPP is handled by the open source program Open Babel [2], which
implements the OpenSMILES specification of SMILES [3]. A number of tools exist
for creating or finding SMILES. Many open access database websites such as
PubChem or chemspider.com provide searchable databases with SMILES strings
included. Note that the EPISuite database [4] popular among environmental
chemists contains some SMILES which do not conform to the standard format, these
structures will not be processed correctly if they are entered into IFSAPP. The
two most common violations are hydrogens outside of brackets, e.g. NH instead of
N[H] or [NH], and counter ions without brackets or disconnected structures,
eg. CC(=O)ONa instead of CC(=O)[O-].[Na+]. In the second example the ONa group
is interpreted as a covalent bond instead of an ion pair, which is incorrect. A
newer alternative also released by the USEPA is the Chemistry Dashboard (found
at https://comptox.epa.gov/dashboard). If an erroneous SMILES is entered into
IFSAPP QSAR predictions will not be made, only the message "SMILES error" will
be displayed.

********************************************************************************
2. QSAR DESCRIPTIONS AND INTERPRETATION
********************************************************************************

In this section:  
-General description of IFSAPP QSARs  
-FHLB - QSAR for fish biotransformation half-life  
-HHLB - QSAR for human biotransformation half-life  
-HHLT - QSAR for human total elimination half-life  
-dSm - QSPR for Entropy of Melting  
-Tm - QSPR for Melting Point

**General description of IFSAPP QSARs**

QSARs are developed using a group contribution model, in which specific
substructures of a molecule are counted and their contributions are determined
by coefficients that have been calculated by multiple linear regression. The
specific methodology used is called Iterative Fragment Selection (IFS) and is
described by the author in detail in references [5] and [6]. The QSARS
conform to best practices as outlined by OECD [7], and includes internal and
external model validation. All models have a defined domain of applicability,
which is described in section 3. INTERPRETING QSAR DOMAIN INFORMATION.

**FHLB - QSAR for fish biotransformation half-life**

The model predicts the base-10 logarithm of the whole-body biotransformation
half-life of chemicals in fish (FHLB) measured in days for a reference 10g fish
at a temperature of 288K. In IFSAPP the model outputs have been converted to
linear (as opposed to log) scale half-lives in hours. The QSAR is trained on a
dataset derived from bioaccumulation tests with biotransformation back-
calculated with a one-compartment PK model. The training and validation datasets
contain 412 and 207 chemicals, respectively. The validation statistics are
r-sq[external] = 0.748 and RMSE of predictions = 0.584, with the data spanning
about 5 log units. Full details are available in reference [5].

**HHLB - QSAR for human biotransformation half-life**

The model predicts the base-10 logarithm of the whole-body biotransformation
half-life of chemicals in humans, (HHLB) measured in hours for a generic 70kg
human. In IFSAPP the model outputs have been converted to linear (as opposed to
log) scale half-lives in hours. As with the FHLB, the QSAR is trained on
biotransformation half-lives backed out of whole body elimination data with a
one compartment PK model. The data is mostly pharmaceuticals measured in humans,
with the remainder environmental contaminants measured in humans. The training
and validation datasets both contain 470 chemicals. The validation statistics
are r-sq[external] = 0.73 and RMSE of predictions = 0.75, with the data spanning
about 7.5 log units. Full details are available in reference [8].

**HHLT - QSAR for human total elimination half-life**

The model predicts the base-10 logarithm of the whole-body total elimination
half-life of chemicals in humans, (HHLT) measured in hours for a generic 70kg
human. In IFSAPP the model outputs have been converted to linear (as opposed to
log) scale half-lives in hours. The QSAR is trained with the same dataset as the
HHLT, except as whole-body elimination half-lives, with some additional data
that could not be transformed into biotransformation half-lives. The training
and validation datasets contain 552 and 553 chemicals, respectively. The
validation statistics are r-sq[external] = 0.72 and RMSE of predictions = 0.70,
with the data spanning about 7.5 log units. Full details are available in
reference [8].

**dSm - QSPR for Entropy of Melting**

The model predicts the entropy of melting (dSm) on a linear scale in units of
kJ/mol. The QSPR is trained on the Jain 2004 dSm data set as described in
reference [10]. The training and validation data sets contain 1056 and 528
chemicals, respectively. The validation statistics are q-sq = 0.792 and std.err.
of prediction = 11.82 kJ/mol, with the data spanning about 175 kJ/mol. Full
details are found in reference [10].

**Tm - QSPR for Melting Point**

The model prdicts the melting point (Tm) on a linear scale in units of K. The
QSPR is trained on the Bradley data set as described in reference [10]. The
training and validation data sets contain 1922 and 964 chemicals, respectively.
The validation statistics are q-sq = 0.658 and std.err. of prediction = 46.9 K,
with the data spanning about 500 K. Full details are found in reference [10].

********************************************************************************
3. INTERPRETING QSAR DOMAIN INFORMATION
********************************************************************************

In this section:  
-What is the domain of applicability?  
-Domain of applicability for the IFSAPP QSARs  
-Interpreting IFSAPP domain outputs

**What is the domain of applicability?**

Having a defined domain of applicability is one of the five OECD guidelines for
developing QSARs. The OECD guidance document [7] gives more information on
domain of applicability in general, and references [5] and [6] describe the
domain definition of the IFS QSARs specifically. The basic problem is that the
chemical space to which a QSAR might be applied for predictions is extremely
large, but datasets available for training QSARs can only ever hope to cover a
fraction of this chemical space. Therefore, applying QSARs means extrapolating
beyond the training set in many cases. The domain tries to quantify how far from
the training dataset the model can make reliable extrapolations.

There is no universally accepted standard method for defining the domain, but
most methods follow a few basic principles. Robust models trained on good data
have a larger domain than uncertain models trained on noisy data. Interpolating
between data points is more reliable than extrapolating outside of the data
range. Making predictions for chemicals similar to those in the training data
is more robust than making predictions for chemicals that are less similar. In
addition the limitations of the model used is a factor. The IFS QSARs are based
on the 2D structure of chemicals, i.e. only atom connectivity is considered, the
three dimensional shape of molecules is not included in the model. Therefore,
in cases where the 3D shape of the molecule might be important, e.g. due to
steric effects, the models are expected to give less accurate results.

**Domain of applicability for the IFSAPP QSARs**

Two different methods of defining the domain are applied to the IFSAPP QSARs in
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

**Interpreting IFSAPP domain outputs**

For each property prediction three domain fields are included:
-error: the estimated standard error of prediction as described above
-UL: the uncertainty level, higher numbers indicate less reliable predictions
-note: explanation for UL
    
Error:
For the FHLB, HHLB, and HHLT QSARs the original estimated standard errors of
prediction on the log scale have been converted to factor errors, also called
dispersion. The error is multiplied by 1.96 and then transformed to linear 
scale. This gives the factor X in the equation:
     P { prediction/X < actual < prediction*X } = 0.95
where the P is the probablility of the actual value being in the interval
defined by prediction/X to prediciont*X, i.e. the 95% prediction interval of the
QSAR. This is analogous to the confidence factor Cf.

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
                            is simply the intercept. If so the error is
                            calculated from any such chemicals in the validation
                            dataset. Otherwise the error is set to the same
                            value as the full validation set.

   5    atom violation,     The simplest domain check is to examine the atoms in
        various messages    the molecule to be predicted and see if it contains
                            unique atom types that are not present in the
                            training dataset. The features of the atom that are
                            unique will be listed, e.g. element type, number of
                            hydrogens attached, etc. The error is set to the
                            same value as UL 3.

   6    prediction less/    For some datasets predictions outside of the
        greater than        experimental range do not make sense. For example,
        largest/smallest    elimination of chemicals from fish or humans can
        value in training   only proceed so rapidly, so predicting faster
        dataset             elimination is not reasonable. In these "bounded"
                            QSARs predictions outside of the experimental range
                            are set to the value of the closest boundary. The
                            error the same as the normal UL 0-3 the molecule
                            belongs to.
</pre>
********************************************************************************
4. SOFTWARE CODING EXPLANATION
********************************************************************************

IFSAPP was written with the open source resources listed below. Because two of
these resources use the GNU General Public License interested users may request
the source code for IFSAPP.

Python: GUI and general coding  
-Python Programming Language, version 2.7.10, http://www.python.org/  
-License: BeOpen.com GPL-compatible license.  

NumPy: Mathematics  
-NumPy, version 1.16.2, http://www.numpy.org/  
-License: Developer specific license. http://www.numpy.org/license.html  

Open Babel: Chemical structure handling  
-The Open Babel Package, version 2.4.1, http://openbabel.org  
-License: GNU GPL. http://www.gnu.org/licenses/gpl.html  

Pyinstaller: Packaging for .exe distribution  
-Pyinstaller, version 3.3.1, http://www.pyinstaller.org/  
-License: GNU GPL. http://www.pyinstaller.org/license.html  

********************************************************************************
5. CHANGE LOG
********************************************************************************

Version B.0 Completed February 2018.  
-All major features implemented, including:  
 -single mode for single calculations with output to IFSAPP  
 -batch mode for reading and writing files  
 -basic error checking for invalid or truncated SMILES  
 -three QSARs implemented: FHLB, HHLB, HHLT  
 -initial version of this user guide
 
Version B.1 Completed 2018  
-Minor bug fixes to domain of applicability testing and output
 for FHLB, HHLB, and HHLT models  
-Added QSPRs for Abraham LSER descriptors E, S, A, B, V and L  
-Updated this user guide to include description of LSER QSPRs, and added
 section to outline known bugs and planned features
 
Version B.2 Completed 2019  
-Removed QSPRs for Abraham LSER descriptors E, S, A, B, V and L to improve load
 time and run speed  
-Added QSPRs for dSm and Tm  
-Updated warning levels to reflect the new ULs from the dSm&Tm paper  
-Updated this user guide to reflect changes

Version 0.0 Completed Nov. 2019  
-Updated code to run on python 3. Python 3 is now required  
-General code cleanup and updating  
-Migrated code to GitHub  
-Converted IFSAPP docstring to readme.md file  

********************************************************************************
6. KNOWN BUGS AND PLANNED FEATURES
********************************************************************************

Identified in Version B.1:  
-Code needs to be optimized to increase the speed of calculations  
-Updated versions of the LSER QSPRs need to be created and published
  
Identified in Version 0.0:  
-The xtxi matrix for domain checking in the model files needs to be re-encoded
 using python 3


********************************************************************************
7. REFERENCES
********************************************************************************

 1. Weininger D., 1988, J. Chem. Inf. Model. 28 (1): 31-6.  
    DOI:10.1021/ci00057a005

 2. O'Boyle N.M., Banck M., James C.A., Morley C., Vandermeersch T., 
    Hutchison G.R. Open Babel: An open chemical toolbox. J. Cheminf. 2011, 3, 33.  
    DOI:10.1186/1758-2946-3-33

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
    organic chemicals. Molecular Informatics, 2019.  
    DOI: 10.1002/minf.201800160

********************************************************************************