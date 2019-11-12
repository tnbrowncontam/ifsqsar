"""smiles_norm.py
Implements a convert function that takes a chemical structure as a SMILES and converts it to a standardized
representation.
"""

import openbabel as ob
import re


def convert(smiles, obconversion=None):
    """Converts passed SMILES to a normalized form.
    By default returns a converted SMILES string, but can also return the new OBMol.
    An existing OBConversion instance can be passed so that one does not need to instantiated."""

    # instantiate OBMol
    mol = ob.OBMol()
    changes = ''

    # error checking of the smiles string
    if smiles[0] in '()#=-+]1234567890':
        return mol, '', 'error reading SMILES'
    elif smiles.count('(') != smiles.count(')') or smiles.count('[') != smiles.count(']'):
        return mol, '', 'error reading SMILES'

    # instantiate obconversion if necessary
    if obconversion is None:
        obconversion = ob.OBConversion()
        obconversion.SetInAndOutFormats('smi', 'can')

    # add -i option to remove isotopic data if necessary
    if obconversion.IsOption('i', obconversion.OUTOPTIONS) is None:
        obconversion.AddOption('i', obconversion.OUTOPTIONS)

    # initialize smarts for finding atoms with permanent charges
    permsmarts = ob.OBSmartsPattern()
    permsmarts.Init('[#7H0v4+!$([#7+]~[#8-]),#8H0a+,#8H0v3+,#8-$([#8-]~[#8H0v3+]),#15H0v4+!$([#15+]~[#8-]),'
                    '#16H0a+,#16H0v3+!$([#16+]~[#8-]),#6H0v3-$([#6-]#[#7H0v4+])]')
    # initialize smarts for finding atoms with any charge
    anysmarts = ob.OBSmartsPattern()
    anysmarts.Init('[*!+0]')
    # initialize smarts for finding atoms with any charge and also having hydrogen atoms attached
    hydrosmarts = ob.OBSmartsPattern()
    hydrosmarts.Init('[*!+0!H0]')

    # count aromatic atoms using string matching
    aromatch = re.compile('(?<!\[[A-Z])[cnosp]')
    aromaticbeforecount = len(re.findall(aromatch, smiles))

    # read smiles into openbabel
    obconversion.ReadString(mol, smiles)

    # check obmol to see if atoms were added, if not then there was a smiles error
    if mol.NumAtoms() == 0:
        return mol, '', 'error reading SMILES'

    if '.' in smiles:
        if changes == '':
            changes = 'salts stripped'
        else:
            changes += ', salts stripped'

    # remove all but the largest contiguous fragment
    obconversion.AddOption('r', obconversion.GENOPTIONS)
    mol.DoTransformations(obconversion.GetOptions(obconversion.GENOPTIONS), obconversion)
    obconversion.RemoveOption('r', obconversion.GENOPTIONS)

    # get lists of permanently charged atoms
    permsmarts.Match(mol)
    permchargeatoms = [i[0] for i in permsmarts.GetUMapList()]

    # get count of non-permanently charged atoms
    anysmarts.Match(mol)
    countanychargeatoms = len((set([i[0] for i in anysmarts.GetUMapList()]).difference(permchargeatoms)))

    # if there charged atoms that are not permanently charged, then modify the molecule to neutralize them
    newsmiles = None
    if countanychargeatoms:
        # setting begin modify prevents openbabel from re-evaluating
        # the structure as modifications are being made
        mol.BeginModify()

        # convert dative bonds
        obconversion.AddOption('b', obconversion.GENOPTIONS)
        mol.DoTransformations(obconversion.GetOptions(obconversion.GENOPTIONS), obconversion)
        obconversion.RemoveOption('b', obconversion.GENOPTIONS)

        # get new list of charged atoms after converting dative bonds
        anysmarts.Match(mol)
        anychargeatoms = list(set([i[0] for i in anysmarts.GetUMapList()]).difference(permchargeatoms))
        if len(anychargeatoms) < countanychargeatoms:
            if changes == '':
                changes = 'dative bonds converted'
            else:
                changes += ', dative bonds converted'

        # get list of non-permanently charged atoms with hydrogens
        hydrosmarts.Match(mol)
        hydrochargeatoms = list(set([i[0] for i in hydrosmarts.GetUMapList()]).difference(permchargeatoms))

        # if there are non-permanently charged atoms with hydrogens then remove hydrogen atoms,
        # they will be added again after setting charges to zero so there are the correct number
        if len(hydrochargeatoms):
            obconversion.AddOption('d', obconversion.GENOPTIONS)
            mol.DoTransformations(obconversion.GetOptions(obconversion.GENOPTIONS), obconversion)
            obconversion.RemoveOption('d', obconversion.GENOPTIONS)
            # the internal atom ids are reassigned automatically after deleting hydrogens, so the
            # SMARTS need to be applied again to get the new atom ids
            permsmarts.Match(mol)
            permchargeatoms = [i[0] for i in permsmarts.GetUMapList()]
            anysmarts.Match(mol)
            anychargeatoms = list(set([i[0] for i in anysmarts.GetUMapList()]).difference(permchargeatoms))

        # set non-permanently charged atoms to zero charge
        for atom in ob.OBMolAtomIter(mol):
            atomid = atom.GetIdx()
            if atomid in anychargeatoms:
                atom.SetFormalCharge(0)
                if changes == '':
                    changes = 'charged atom(s) neutralized'
                elif 'atoms neutralized' not in changes:
                    changes += ', charged atom(s) neutralized'

        # ending molecule modification causes openbabel to re-evaluate the
        # valence of the atoms with their charges now set to zero
        mol.EndModify()

        # re-add hydrogen atoms with the new valence if they were deleted
        # note this does not necessarily cause them to be shown in the SMILES
        if len(hydrochargeatoms):
            obconversion.AddOption('h', obconversion.GENOPTIONS)
            mol.DoTransformations(obconversion.GetOptions(obconversion.GENOPTIONS), obconversion)
            obconversion.RemoveOption('h', obconversion.GENOPTIONS)

    # check for permanently charged atoms
    if len(permchargeatoms) > 0:
        if changes == '':
            changes = 'structure contains permanently charged atoms'
        elif 'structure contains permanently charged atoms' not in changes:
            changes += ', structure contains permanently charged atoms'

    # count the aromatic atoms to check if aromaticity was broken by manipulations
    newsmiles = obconversion.WriteString(mol).strip()
    aromaticaftercount = len(re.findall(aromatch, newsmiles))
    if aromaticaftercount < aromaticbeforecount:
        changes = 'aromaticity broken'
        newsmiles = ''

    # return results
    return mol, newsmiles, changes


# list for testing
test_list = [
    # (N+) ammonium
    '[NH4+]',
    # (N+) methyl amine
    'C[NH3+]',
    # (N+) trimethyl amine
    'C[NH+](C)C',
    # (N+) nitrobenzene (dative bond representation)
    'c1ccccc1[N+]([O-])=O',
    # (P+) (dative bond representation)
    'c1ccccc1[P+][O-]',
    # (P+) Diethyl hydrogen phosphate (dative bond representation)
    'CCO[P+]([O-])(O)OCC',
    # (S+) Dimethyl sulfite (dative bond representation)
    'CO[S+](-[O-])OC',
    # (n+) 1-Methylnicotinamide
    'C[n+]1cccc(c1)C(=O)N',
    # (N+) tetramethyl quaternary amine
    'C[N+](C)(C)C',
    # (N+) Gentian violet cation- CAS 7438-46-2
    'C[N+](=C1C=CC(=C(c2ccc(cc2)N(C)C)c2ccc(cc2)N(C)C)C=C1)C',
    # (N+) Methyl isocyanide - converting to non-dative representation results
    # in a "quadruple bond" which is not realistic. openbabel will convert and
    # represent the bond as "$" but this bond will not be matched by any
    # smarts strings, resulting in possible errors. should not be converted
    'C[N+]#[C-]',
    # (o+) Pelargonidin (dye) - 134-04-3
    'Oc1ccc(cc1)c1[o+]c2cc(O)cc(c2cc1O)O',
    # (O+) methoxymethylidene(methyl)oxidanium
    'COC=[O+]C',
    # (O+) ozone - converting to non-dative representation results in two double
    # bonds to an oxygen, which is not realistic. should not convert
    '[O-][O+]=O',
    # (s+) 2-Phenyl-1-benzothiopyrylium perchlorate - CAS 22956-28-1
    'c1ccc(cc1)c1ccc2c([s+]1)cccc2.[O-]Cl(=O)(=O)=O',
    # (S+) Dimethylpropiothetin - CAS 7314-30-9
    'C[S+](C)CCC(=O)[O-]',
    # (p+) (P+?) 2-methyl-3-(1-methylphosphinin-1-ium-2-yl)-3,4-dihydropyrazole
    # SMARTS search for [p+] (aromatic phosphorous) in pubchem returns this and
    # other chemicals, but neither pubchem or openbabel display the structures as
    # aromatic
    'CN1C(CC=N1)C2=CC=CC=[P+]2C',
    # (P+) Chlorphonium chloride (pesticide) - CAS 115-78-6
    'CCCC[P+](CCCC)(CCCC)CC1=C(C=C(C=C1)Cl)Cl.[Cl-]',
    # (P+) bis(4-methylanilino)-oxophosphanium
    'CC1=CC=C(C=C1)N[P+](=O)NC2=CC=C(C=C2)C',
    # cis-1,2-dichloroethene
    'Cl/C=C\Cl',
    # false aromaticity
    'c1cccc1',
    # erroneous smiles
    'erroneous smiles',
    'Cc1ncc(N(=O)=O)n1CCO',
    'S(c(ccc(N)c1)c1)c(ccc(N)c2)c2',

]

if __name__ == '__main__':
    for smiles in test_list:
        m, s, n = convert(smiles)
        print(smiles, s, n)
