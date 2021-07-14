"""smiles_norm.py
Implements a convert function that takes a chemical structure as a SMILES and converts it to a standardized
representation.
"""

from openbabel import openbabel as ob
import re

# initialize smarts for finding atoms with any charge
anysmarts = ob.OBSmartsPattern()
anysmarts.Init('[*!+0]')

# initialize smarts for handling silicon implicit hydrogens
silicon4 = ob.OBSmartsPattern()
silicon4.Init('[#14v4]')
silicon3 = ob.OBSmartsPattern()
silicon3.Init('[#14v3]')
silicon2 = ob.OBSmartsPattern()
silicon2.Init('[#14v2]')
silicon1 = ob.OBSmartsPattern()
silicon1.Init('[#14v1]')
silicon0 = ob.OBSmartsPattern()
silicon0.Init('[#14v0]')

# initialize smarts for handling isocyanides
isocyanide = ob.OBSmartsPattern()
isocyanide.Init('[CD1-]#[O,N;+]')

# initialize smarts for organic structure checking
organicsmarts = ob.OBSmartsPattern()
organicsmarts.Init('[*!$([#1,#5,#6,#7,#8,#9,#14,#15,#16,#17,#35,#53])]')
carbonsmarts = ob.OBSmartsPattern()
carbonsmarts.Init('[#6]')

# initialize re for matching aromatic atoms
aromatch = re.compile('(?<!\[[A-Z])[cnosp]')

def convertsmiles(smiles, obconversion=None, standardize=True, neutralize=True, filtertype='organic'):
    """Converts passed SMILES to a normalized form.
    Returns the new OBMol, a normalized SMILES string and notes on changes made.
    An existing OBConversion instance can be passed so that one does not need to instantiated."""

    # instantiate OBMol
    mol = ob.OBMol()
    changes = []

    # error checking of the smiles string
    smiles = smiles.lstrip().rstrip()
    if smiles[0] in '()#=-+]1234567890%':
        return mol, [''], 'error reading SMILES: invalid first character in SMILES'
    elif smiles.count('(') != smiles.count(')') or smiles.count('[') != smiles.count(']'):
        return mol, [''], 'error reading SMILES: missing brackets'

    # instantiate obconversion if necessary
    if obconversion is None:
        obconversion = ob.OBConversion()
        obconversion.SetInAndOutFormats('smi', 'can')

    # read smiles into openbabel and save original smiles
    obconversion.ReadString(mol, smiles)

    # check obmol to see if atoms were added, if not then there was a smiles error
    if mol.NumAtoms() == 0:
        return mol, [''], 'error reading SMILES: no atoms loaded, check input'

    # handle silicon valence, fill up to valence 4 with implicit hydrogens
    silicon0.Match(mol)
    for atom in silicon0.GetUMapList():
        mol.GetAtom(atom[0]).SetImplicitHCount(4)
    silicon1.Match(mol)
    for atom in silicon1.GetUMapList():
        mol.GetAtom(atom[0]).SetImplicitHCount(3)
    silicon2.Match(mol)
    for atom in silicon2.GetUMapList():
        mol.GetAtom(atom[0]).SetImplicitHCount(2)
    silicon3.Match(mol)
    for atom in silicon3.GetUMapList():
        mol.GetAtom(atom[0]).SetImplicitHCount(1)
    silicon4.Match(mol)
    for atom in silicon4.GetUMapList():
        mol.GetAtom(atom[0]).SetImplicitHCount(0)

    # handle isocyanides
    isocyanide.Match(mol)
    for match in isocyanide.GetUMapList():
        # determine carbon and heteroatoms
        atom1 = mol.GetAtom(match[0])
        atom2 = mol.GetAtom(match[1])
        if atom1.GetAtomicNum() == 6:
            carbon = atom1
            hetero = atom2
        else:
            carbon = atom2
            hetero = atom1
        bond = mol.GetBond(match[0], match[1])
        # change atom and bond characteristics
        carbon.SetFormalCharge(0)
        hetero.SetFormalCharge(0)
        bond.SetBondOrder(2)

    # update mol for manual changes
    mol.DoTransformations(obconversion.GetOptions(obconversion.GENOPTIONS), obconversion)

    # get lists of charged atoms
    anysmarts.Match(mol)
    anychargeatoms = list(anysmarts.GetUMapList())

    # convert dative bonds
    originalmol = ob.OBMol(mol)
    if standardize or neutralize:
        molcopy = ob.OBMol(mol)
        obconversion.AddOption('b', obconversion.GENOPTIONS)
        molcopy.DoTransformations(obconversion.GetOptions(obconversion.GENOPTIONS), obconversion)
        obconversion.RemoveOption('b', obconversion.GENOPTIONS)
        anysmarts.Match(molcopy)
        newanychargeatoms = list(anysmarts.GetUMapList())
        newsmiles = obconversion.WriteString(molcopy)
        # make sure this didn't create any invalid bonds i.e. quadruple bond "$" then update original mol
        if len(newanychargeatoms) < len(anychargeatoms) and '$' not in newsmiles:
            changes.append('dative bonds converted')
            mol = molcopy
            originalmol = ob.OBMol(mol)
        # check for a mix of valid and invalid dative bond types
        elif '$' in newsmiles:
            count = newsmiles.count('$')
            converted = (len(anychargeatoms) - len(newanychargeatoms)) / 2
            if count < converted:
                return mol, [''], 'error standardizing SMILES: manual dative bond conversion required, structure contains both valid and invalid dative bonds'

    # neutralize the structure and see if anything changed
    molcopy = ob.OBMol(mol)
    obconversion.AddOption('neutralize', obconversion.GENOPTIONS)
    molcopy.DoTransformations(obconversion.GetOptions(obconversion.GENOPTIONS), obconversion)
    obconversion.RemoveOption('neutralize', obconversion.GENOPTIONS)
    anysmarts.Match(molcopy)
    newanychargeatoms = list(anysmarts.GetUMapList())
    if not neutralize and len(newanychargeatoms) < len(anychargeatoms):
        changes.append('structure contains neutralizable atoms')
    elif neutralize and len(newanychargeatoms) < len(anychargeatoms):
        changes.append('charged atoms neutralized')
        mol = molcopy

    if len(newanychargeatoms) > 0:
        changes.append('structure contains permanently charged atoms')

    # save standardized, canonical structure including counter-ions, and chiral marks
    canonicalchiralsmiles = obconversion.WriteString(mol).strip()

    # save standardized, canonical structure including counter-ions, but no chiral marks
    if obconversion.IsOption('i', obconversion.OUTOPTIONS) is None:
        obconversion.AddOption('i', obconversion.OUTOPTIONS)
    canonicalnochiralsmiles = obconversion.WriteString(mol).strip()
    obconversion.RemoveOption('i', obconversion.OUTOPTIONS)

    # remove all but the largest contiguous fragment
    if (standardize or neutralize) and '.' in smiles:
        changes.append('salts stripped')
        obconversion.AddOption('r', obconversion.GENOPTIONS)
        mol.DoTransformations(obconversion.GetOptions(obconversion.GENOPTIONS), obconversion)
        originalmol.DoTransformations(obconversion.GetOptions(obconversion.GENOPTIONS), obconversion)
        obconversion.RemoveOption('r', obconversion.GENOPTIONS)

    # check for non-organic atoms
    organicsmarts.Match(mol)
    carbonsmarts.Match(mol)
    if filtertype == 'organic' and (len(organicsmarts.GetUMapList()) > 0 or len(carbonsmarts.GetUMapList()) == 0):
        return mol, [''], 'input error: structure is not organic'

    # save salt-stripped canonical structure with chiral marks
    parentchiralsmiles = obconversion.WriteString(mol).strip()

    # save salt-stripped canonical structure with no chiral marks, both neutralized and charged forms
    if obconversion.IsOption('i', obconversion.OUTOPTIONS) is None:
        obconversion.AddOption('i', obconversion.OUTOPTIONS)
    parentnochiralsmiles = obconversion.WriteString(mol).strip()
    chargednochiralsmiles = obconversion.WriteString(originalmol).strip()
    obconversion.RemoveOption('i', obconversion.OUTOPTIONS)

    # count the aromatic atoms to check if aromaticity was broken by manipulations
    aromaticbeforecount = len(re.findall(aromatch, smiles))
    aromaticaftercount = len(re.findall(aromatch, parentnochiralsmiles))
    if aromaticaftercount < aromaticbeforecount:
        changes.append('error reading SMILES: aromaticity broken')
        return mol, [''], ', '.join(set(changes))

    # return results
    return mol, [parentnochiralsmiles, chargednochiralsmiles, parentchiralsmiles, canonicalnochiralsmiles, canonicalchiralsmiles], ', '.join(set(changes))


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
    # Phenol (N+) nitrobenzene (dative bond representation)
    'c1cc([O-])ccc1[N+]([O-])=O',
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
    # mix of dative bonds that can and cannot be converted
    'O=[N+]([O-])C[N+]#[C-]',
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
    'N(=O)(=O)c(ccc(N)c1)c1',
    'CC34CCC1C(CCc2cc(O)ccc12)C3CCC4O',
    '[Sn]c1ccccc1',
    'CCCCC[O-].[Na+]',
    'CCCCC[O-].[K+].[Br-]',
    'CC[C@H](C)CCC[NH3+].[Cl-]',
    'Cc1ccc(cc1)C#Cc2ncn3c2COc4ccccc34',
    'C1C[N+]2=C(C=CC=C2)C2=CC=CC=[N+]12',
    'CCC(C1OC2(C=CC1C)OC1CC=C(C)C(OC3CC(OC)C(C(O3)C)OC3CC(OC)C(C(O3)C)NC)C(C)C=CC=C3C4(C(C(=O)OC(C2)C1)C=C(C)C(C4OC3)O)O)C',
    '[SiH3]CCC',
    '[Si]CCC',
    '[C-]#[N+]c1ccccc1'
]

if __name__ == '__main__':
    for smiles in test_list:
        m, s, n = convertsmiles(smiles, standardize=True, neutralize=True, filtertype='organic')
        if len(s) > 1:
            print(smiles, s[0], s[1], s[2], s[3], n)
        else:
            print(smiles, s[0], n)
