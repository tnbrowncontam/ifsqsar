"""smiles_norm.py
Implements a convert function that takes a chemical structure as a SMILES and converts it to a standardized
representation.
"""

from openbabel import openbabel as ob
import re

# initialize smarts for handling silicon implicit hydrogens
silicon3 = ob.OBSmartsPattern()
silicon3.Init('[#14v3H0]')
silicon2 = ob.OBSmartsPattern()
silicon2.Init('[#14v2H0]')
silicon1 = ob.OBSmartsPattern()
silicon1.Init('[#14v1H0]')
silicon0 = ob.OBSmartsPattern()
silicon0.Init('[#14v0H0]')

# initialize smarts for handling isocyanides
isocyanide = ob.OBSmartsPattern()
isocyanide.Init('[CD1-]#[O,N;+]')

# initialize smarts for handling azides
azide = ob.OBSmartsPattern()
azide.Init('[N-]=[N+]=N')

# initialize smarts for finding atoms with any charge
chargedsmarts = ob.OBSmartsPattern()
chargedsmarts.Init('[*!+0]')

# initialize smarts for organic structure filters
anyatom = ob.OBSmartsPattern()
anyatom.Init('[*]')
organicatom = ob.OBSmartsPattern()
organicatom.Init('[#1,#5,#6,#7,#8,#9,#14,#15,#16,#17,#35,#53]')
carbonatom = ob.OBSmartsPattern()
carbonatom.Init('[#6]')
inorganicatom = ob.OBSmartsPattern()
inorganicatom.Init('[!#1H0]')
organometallicatom = ob.OBSmartsPattern()
organometallicatom.Init('[H0!$([#1,#5,#6,#7,#8,#9,#14,#15,#16,#17,#35,#53])$(*~[#6])!$(*~[!#6])]')


# initialize re for matching aromatic atoms
aromatch = re.compile('(?<!\[[A-Z])[cnosp]')

def convertsmiles(smiles, obconversion=None, neutralize=True, filtertype='not inorganic'):
    """Converts passed SMILES to a normalized form.
    Returns the new OBMol, a normalized SMILES string and notes on changes made.
    An existing OBConversion instance can be passed so that one does not need to instantiated."""

    # instantiate OBMol
    mol = ob.OBMol()
    changes = []

    # error checking of the smiles string
    smiles = smiles.lstrip().rstrip()
    if smiles[0] in '()#=-+]1234567890%':
        return mol, '', 'error reading SMILES: invalid first character in SMILES'
    elif smiles.count('(') != smiles.count(')') or smiles.count('[') != smiles.count(']'):
        return mol, '', 'error reading SMILES: missing brackets'

    # instantiate obconversion if necessary
    if obconversion is None:
        obconversion = ob.OBConversion()

    # handle silicon valence, fill up to valence 4 with implicit hydrogens
    if '[Si]' in smiles:
        obconversion.SetInAndOutFormats('smi', 'smi')
        obconversion.ReadString(mol, smiles)
        silicon0.Match(mol)
        for atom in silicon0.GetUMapList():
            thisatom = mol.GetAtom(atom[0])
            thisatom.SetImplicitHCount(4)
            mol.AddHydrogens(thisatom)
        silicon1.Match(mol)
        for atom in silicon1.GetUMapList():
            thisatom = mol.GetAtom(atom[0])
            thisatom.SetImplicitHCount(3)
            mol.AddHydrogens(thisatom)
        silicon2.Match(mol)
        for atom in silicon2.GetUMapList():
            thisatom = mol.GetAtom(atom[0])
            thisatom.SetImplicitHCount(2)
            mol.AddHydrogens(thisatom)
        silicon3.Match(mol)
        for atom in silicon3.GetUMapList():
            thisatom = mol.GetAtom(atom[0])
            thisatom.SetImplicitHCount(1)
            mol.AddHydrogens(thisatom)
        obconversion.AddOption('h', obconversion.OUTOPTIONS)
        vsmiles = obconversion.WriteString(mol).strip()
        obconversion.RemoveOption('h', obconversion.OUTOPTIONS)
    else:
        vsmiles = smiles

    # load structure
    obconversion.SetInAndOutFormats('smi', 'smi')
    obconversion.ReadString(mol, vsmiles)

    # check obmol to see if atoms were added, if not then there was a smiles error
    if mol.NumAtoms() == 0:
        return mol, '', 'SMILES error: no atoms loaded, check input'

    # remove all but the largest contiguous fragment
    if '.' in vsmiles:
        print('.')
        obconversion.AddOption('r', obconversion.GENOPTIONS)
        mol.DoTransformations(obconversion.GetOptions(obconversion.GENOPTIONS), obconversion)
        obconversion.RemoveOption('r', obconversion.GENOPTIONS)
        changes.append('salts stripped')

    # check for filters
    organometallicatom.Match(mol)
    if filtertype in ['organic', 'not inorganic']:
        anyatom.Match(mol)
        organicatom.Match(mol)
        inorganicatom.Match(mol)
        carbonatom.Match(mol)
        isinorganic = False
        if len(carbonatom.GetUMapList()) == 0 or \
           (len(inorganicatom.GetUMapList()) == len(anyatom.GetUMapList()) and len(inorganicatom.GetUMapList()) <= 3):
            isinorganic = True
        if filtertype == 'organic' and (isinorganic or
           len(organicatom.GetUMapList()) < len(anyatom.GetUMapList())):
            return ob.OBMol(), '', 'SMILES error: structure fails organic filter'
        elif filtertype == 'not inorganic' and (isinorganic or
             len(organicatom.GetUMapList()) + len(organometallicatom.GetUMapList()) < len(anyatom.GetUMapList())):
            return ob.OBMol(), '', 'SMILES error: structure fails not inorganic filter'

    # do not inchify organometallics, it disconnects the structures
    if len(organometallicatom.GetUMapList()) == 0:
        # inchify smiles and read into obmol
        obconversion.AddOption('I', obconversion.OUTOPTIONS)
        if obconversion.IsOption('i', obconversion.OUTOPTIONS) is None:
            obconversion.AddOption('i', obconversion.OUTOPTIONS)
        inchismiles = obconversion.WriteString(mol).strip()
        obconversion.RemoveOption('I', obconversion.OUTOPTIONS)
        obconversion.RemoveOption('i', obconversion.OUTOPTIONS)
        obconversion.SetInAndOutFormats('smi', 'can')
        obconversion.ReadString(mol, inchismiles)
    else:
        if obconversion.IsOption('i', obconversion.OUTOPTIONS) is None:
            obconversion.AddOption('i', obconversion.OUTOPTIONS)
        inchismiles = obconversion.WriteString(mol).strip()
        obconversion.RemoveOption('i', obconversion.OUTOPTIONS)
        obconversion.SetInAndOutFormats('smi', 'can')
        obconversion.ReadString(mol, inchismiles)
        # convert dative bonds because smiles was not inchified
        obconversion.AddOption('b', obconversion.GENOPTIONS)
        mol.DoTransformations(obconversion.GetOptions(obconversion.GENOPTIONS), obconversion)
        obconversion.RemoveOption('b', obconversion.GENOPTIONS)

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

    # handle azides
    azide.Match(mol)
    for match in azide.GetUMapList():
        # determine negative and positive charged atoms
        atom1 = mol.GetAtom(match[0])
        atom2 = mol.GetAtom(match[1])
        atom3 = mol.GetAtom(match[2])
        if atom1.GetFormalCharge() == -1:
            neg = atom1
        elif atom2.GetFormalCharge() == -1:
            neg = atom2
        else:
            neg = atom3
        if atom1.GetFormalCharge() == 1:
            pos = atom1
        elif atom2.GetFormalCharge() == 1:
            pos = atom2
        else:
            pos = atom3
        bond = mol.GetBond(neg.GetIdx(), pos.GetIdx())
        # change atom and bond characteristics
        pos.SetFormalCharge(0)
        neg.SetFormalCharge(0)
        bond.SetBondOrder(3)

    # update mol for manual changes
    mol.DoTransformations(obconversion.GetOptions(obconversion.GENOPTIONS), obconversion)

    # get lists of charged atoms
    chargedsmarts.Match(mol)
    anychargeatoms = list(chargedsmarts.GetUMapList())

    # neutralize the structure and see if anything changed
    molcopy = ob.OBMol(mol)
    obconversion.AddOption('neutralize', obconversion.GENOPTIONS)
    molcopy.DoTransformations(obconversion.GetOptions(obconversion.GENOPTIONS), obconversion)
    obconversion.RemoveOption('neutralize', obconversion.GENOPTIONS)
    chargedsmarts.Match(molcopy)
    newanychargeatoms = list(chargedsmarts.GetUMapList())
    if not neutralize and len(newanychargeatoms) < len(anychargeatoms):
        changes.append('structure contains neutralizable atoms')
    elif neutralize and len(newanychargeatoms) < len(anychargeatoms):
        changes.append('charged atoms neutralized')
        mol = molcopy

    if len(newanychargeatoms) > 0 and not neutralize:
        changes.append('structure contains permanently charged atoms')
    elif len(newanychargeatoms) > 0 and neutralize:
        return ob.OBMol(), '', 'SMILES error: structure fails neutralize filter, contains permanently charged atoms'

    # output normalized smiles
    normsmiles = obconversion.WriteString(mol).strip()

    # count the aromatic atoms to check if aromaticity was broken by manipulations
    aromaticbeforecount = len(re.findall(aromatch, smiles))
    aromaticaftercount = len(re.findall(aromatch, normsmiles))
    if aromaticaftercount < aromaticbeforecount:
        changes.append('SMILES error: aromaticity broken')
        return ob.OBMol(), '', ', '.join(changes)

    # return results
    return mol, normsmiles, ', '.join(changes)

