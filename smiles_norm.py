"""smiles_norm.py
Implements a convert function that takes a chemical structure as a SMILES and converts it to a standardized
representation.
"""

import pkg_resources
if pkg_resources.get_distribution('openbabel').version.split('.')[0] == '3':
    from openbabel import openbabel as ob
    obver = '3'
else:
    import openbabel as ob
    obver = '2'
del pkg_resources
import re

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

# initialize smarts for organic structure checking
organicsmarts = ob.OBSmartsPattern()
organicsmarts.Init('[*!$([#1,#5,#6,#7,#8,#9,#14,#15,#16,#17,#35,#53])]')
carbonsmarts = ob.OBSmartsPattern()
carbonsmarts.Init('[#6]')

# initialize smarts for neutralizable atoms
neutralizesmarts = ob.OBSmartsPattern()
neutralizesmarts.Init('[#6,#7,#8,#15,#16;!+0]')

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
    if smiles[0] in '()#=-+]1234567890':
        return mol, [''], 'error reading SMILES: invalid first character in SMILES'
    elif smiles.count('(') != smiles.count(')') or smiles.count('[') != smiles.count(']'):
        return mol, [''], 'error reading SMILES: missing brackets'

    # instantiate obconversion if necessary
    if obconversion is None:
        obconversion = ob.OBConversion()
        obconversion.SetInAndOutFormats('smi', 'can')

    # read smiles into openbabel
    obconversion.ReadString(mol, smiles)

    # check obmol to see if atoms were added, if not then there was a smiles error
    if mol.NumAtoms() == 0:
        return mol, [''], 'error reading SMILES: no atoms loaded, check input'

    # get lists of charged atoms
    permsmarts.Match(mol)
    permchargeatoms = [i[0] for i in permsmarts.GetUMapList()]
    anysmarts.Match(mol)
    anychargeatoms = list(set([i[0] for i in anysmarts.GetUMapList()]).difference(permchargeatoms))

    # check for permanently charged atoms
    if len(permchargeatoms) > 0:
        changes.append('structure contains permanently charged atoms')

    # convert dative bonds
    if standardize or neutralize:
        lastmol = ob.OBMol(mol)
        obconversion.AddOption('b', obconversion.GENOPTIONS)
        mol.DoTransformations(obconversion.GetOptions(obconversion.GENOPTIONS), obconversion)
        obconversion.RemoveOption('b', obconversion.GENOPTIONS)
        # get new lists of charged atoms after converting dative bonds
        lastpermchargeatoms = permchargeatoms.copy()
        permsmarts.Match(mol)
        permchargeatoms = [i[0] for i in permsmarts.GetUMapList()]
        lastanychargeatoms = anychargeatoms.copy()
        anysmarts.Match(mol)
        anychargeatoms = list(set([i[0] for i in anysmarts.GetUMapList()]).difference(permchargeatoms))
        # nonperm charges reduced, perm charges not reduced: add dative bonds note
        if len(anychargeatoms) < len(lastanychargeatoms) and len(permchargeatoms) == len(lastpermchargeatoms):
            changes.append('dative bonds converted')
        # nonperm charges not reduced, perm charges reduced: revert to previous obmol
        elif len(anychargeatoms) == len(lastanychargeatoms) and len(permchargeatoms) < len(lastpermchargeatoms):
            mol = lastmol
            permchargeatoms = lastpermchargeatoms
            anychargeatoms = lastanychargeatoms
        # both nonperm and perm charges reduced: flag for manual curation
        elif len(anychargeatoms) < len(lastanychargeatoms) and len(permchargeatoms) < len(lastpermchargeatoms):
            return mol, [''], 'error standardizing SMILES: manual dative bond conversion required, structure contains unresolvable mix of permanent and non-permanent charges'

    # if there charged atoms that are not permanently charged, then modify the molecule to neutralize them
    molcopy = ob.OBMol(mol)
    if not neutralize and len(anychargeatoms):
        changes.append('structure contains neutralizable atoms')
    elif neutralize and len(anychargeatoms):
        # get list of atom types to be neutralized
        neutralizesmarts.Match(mol)
        neutralizableatoms = [i[0] for i in neutralizesmarts.GetUMapList()]

        # setting begin modify prevents openbabel from re-evaluating
        # the structure as modifications are being made
        mol.BeginModify()

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
                formalcharge = atom.GetFormalCharge()
                atom.SetFormalCharge(0)
                if obver == '3' and atomid in neutralizableatoms:
                    newHcount = atom.GetImplicitHCount() - formalcharge
                    assert newHcount >= 0
                    atom.SetImplicitHCount(newHcount)
                changes.append('structure contains neutralizable atoms')

        # ending molecule modification causes openbabel to re-evaluate the
        # valence of the atoms with their charges now set to zero
        # (this recalculates implicit valence in openbabel 2X only, not in 3X!)
        mol.EndModify()

        # re-add hydrogen atoms with the new valence if they were deleted
        # note this does not necessarily cause them to be shown in the SMILES
        if len(hydrochargeatoms):
            obconversion.AddOption('h', obconversion.GENOPTIONS)
            mol.DoTransformations(obconversion.GetOptions(obconversion.GENOPTIONS), obconversion)
            obconversion.RemoveOption('h', obconversion.GENOPTIONS)

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
        molcopy.DoTransformations(obconversion.GetOptions(obconversion.GENOPTIONS), obconversion)
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
    chargednochiralsmiles = obconversion.WriteString(molcopy).strip()
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
    'CC[C@H](C)CCC[NH3+].[Cl-]',
    'Cc1ccc(cc1)C#Cc2ncn3c2COc4ccccc34',
    'C1C[N+]2=C(C=CC=C2)C2=CC=CC=[N+]12'
]

if __name__ == '__main__':
    for smiles in test_list:
        m, s, n = convertsmiles(smiles, standardize=True, neutralize=True, filtertype='organic')
        if len(s) > 1:
            print(smiles, s[0], s[1], s[2], s[3], n)
        else:
            print(smiles, s[0], n)
