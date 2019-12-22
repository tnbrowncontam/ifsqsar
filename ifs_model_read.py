"""ifs_model_read.py
Reads model files from the development code for group contribution QSARs (IFS) developed by Trevor N. Brown.
Applies the QSARs and checks domain of applicability for single structures as SMILES or for a file of structures.
"""

import openbabel as ob
import numpy as np
import base64
import pickle


def normcdfapprox(x):
    """Returns normal distribution CDF."""
    # approximation from:
    # Hector Vazquez-Leal, Roberto Castaneda-Sheissa, Uriel Filobello-Nino,
    # Arturo Sarmiento-Reyes, and Jesus Sanchez Orea, "High Accurate Simple
    # Approximation of Normal Distribution Integral," Mathematical Problems
    # in Engineering, vol. 2012, Article ID 124029, 22 pages, 2012.
    # doi:10.1155/2012/124029
    return 1./(np.exp(-x*358./23.+111*np.arctan(x*37./294.))+1.)


def calculate_fragment_similarity(counts_array_i, counts_array_j, stdev_array, simil_cut=None):
    """Calculate Tanimoto similarity coefficient between two arrays of fragment counts."""
    # all b-values = 1 meaning:
    #  1) a-values * b-values = a-values
    #  2) b-values**2 = b-values
    #  3) sum(b-values) = # non-zero elements
    #  so only the number of non-zero elements needs to be known
    b = float(np.logical_or(counts_array_i, counts_array_j).sum())
    if b == 0.:
        tanimoto = 0.
    else:
        # find non-zero a-values
        amask1 = np.logical_and(counts_array_i, counts_array_j)
        a1 = amask1.sum()
        # at this point sum(a)/sum(b) is the maximum possible CSS for this
        # pair of chemicals, so if it's less than the similarity cutoff
        # then skip the rest of the calculation
        if simil_cut is not None and a1/b < simil_cut:
            tanimoto = 0.
        else:
            # some a-values will be equal to one because the fragment
            # counts are the same (a1s). only the number needs to be
            # known for similar reasons as the b-values
            amask2 = np.logical_and(amask1, counts_array_i != counts_array_j)
            a1 -= amask2.sum()
            # calculate non-zero and non-one a-values if any.
            # for fragments with stdev = 0, the a value is = 0,
            # so do not include these in the array
            zmask = stdev_array != 0
            aa = np.abs(counts_array_i[amask2*zmask] - counts_array_j[amask2*zmask]) / stdev_array[amask2*zmask]
            a = 1 - (2 * normcdfapprox(aa) - 1)
            # calculate tanimoto => sum(ab) / [sum(a2) + sum(b2) - sum(ab)]
            asum = a1 + a.sum()
            tanimoto = asum / ((a1+(a**2).sum()) + b - asum)
    return tanimoto


class Model:
    """Class that loads an arbitrary Model from a Model file and applies it
    to any molecules passed to it as an openbabel mol"""
    
    def __init__(self, model_name):
        """Load Model from file, define variables."""
        
        # define Model namespace
        self.model_namespace = {}

        # read Model file
        modelfile = open(model_name, 'r')
        file_lines = []
        modelcommands = False
        chemical_lines = []
        trainset = False
        
        for line in modelfile:
            splitline = line.rstrip('\n').split('\t')
            # read in Model command sequence
            if splitline[0] == '<startmodel>':
                modelcommands = True
            elif splitline[0] == '<endmodel>':
                modelcommands = False
            elif splitline[0] == '<startdomain>':
                modelcommands = True
            elif splitline[0] == '<enddomain>':
                modelcommands = False
            elif splitline[0] == '<startdomaindata>':
                modelcommands = True
            elif splitline[0] == '<enddomaindata>':
                modelcommands = False
            elif modelcommands:
                file_lines.append(splitline)
            # read in chemical data lines
            elif splitline[0] == '<startdataset>':
                trainset = True
            elif splitline[0] == '<enddataset>':
                trainset = False
            elif trainset:
                chemical_lines.append(splitline)
        modelfile.close()

        # read settings
        for line in file_lines:
            if line[0] == '<setting>':
                self.model_namespace['s'+line[1]] = line[2]
        
        # read in training dataset count arrays and values
        self.model_namespace['train_arrays'] = {}
        self.model_namespace['train_values'] = {}
        self.model_namespace['train_array'] = None
        self.model_namespace['train_value'] = None
        self.model_namespace['train_value_similarity'] = None
        self.model_namespace['valid_array'] = None
        self.model_namespace['valid_value'] = None
        header = True
        value_index = 999999
        value_similarity_index = 999999
        value_similarity = 999999
        model_descriptors_index = 999999
        fit_test_index = 999999
        train_validate_index = 999999
        for line in chemical_lines:
            if header:
                fit_test_index = line.index('fit_test')
                train_validate_index = line.index('train_validate')
                for c in range(len(line)):
                    if line[c][:line[c].find('|')] == 'model_descriptors':
                        model_descriptors_index = c
                    elif line[c] == 'value_similarity':
                        value_similarity_index = c
                value_index = line.index(self.model_namespace['svalue_name'])
                header = False
            elif model_descriptors_index != 999999:
                countlist = line[model_descriptors_index].split(',')
                if countlist == ['']:
                    countlist = []
                for c in range(len(countlist)):
                    countlist[c] = float(countlist[c])
                value = float(line[value_index])
                fit_test = int(line[fit_test_index])
                train_validate = line[train_validate_index]
                if value_similarity_index != 999999:
                    value_similarity = line[value_similarity_index]
                # load train values and counts by fold for similarity readacross
                if fit_test != 0 and fit_test not in self.model_namespace['train_arrays']:
                    self.model_namespace['train_arrays'][fit_test] = np.array([countlist], dtype=float)
                    self.model_namespace['train_values'][fit_test] = np.array([[value]], dtype=float)
                elif fit_test != 0 and fit_test in self.model_namespace['train_arrays']:
                    self.model_namespace['train_arrays'][fit_test] = np.vstack((
                        self.model_namespace['train_arrays'][fit_test], np.array([countlist], dtype=float)))
                    self.model_namespace['train_values'][fit_test] = np.vstack((
                        self.model_namespace['train_values'][fit_test], np.array([[value]], dtype=float)))
                # load all values and counts each as one chunk for domain testing
                if train_validate == 'train' and self.model_namespace['train_array'] is None:
                    self.model_namespace['train_array'] = np.array([countlist], dtype=float)
                    self.model_namespace['train_value'] = np.array([[value]], dtype=float)
                    if value_similarity_index != 999999:
                        self.model_namespace['train_value_similarity'] = np.array([[value_similarity]], dtype=float)
                elif train_validate == 'train' and self.model_namespace['train_array'] is not None:
                    self.model_namespace['train_array'] = np.vstack((self.model_namespace['train_array'],
                                                                     np.array([countlist], dtype=float)))
                    self.model_namespace['train_value'] = np.vstack((self.model_namespace['train_value'],
                                                                     np.array([[value]], dtype=float)))
                    if value_similarity_index != 999999:
                        self.model_namespace['train_value_similarity'] = np.vstack((
                            self.model_namespace['train_value_similarity'],
                            np.array([[value_similarity]], dtype=float)))
                elif train_validate == 'validate' and type(self.model_namespace['valid_array']) != np.ndarray:
                    self.model_namespace['valid_array'] = np.array([countlist], dtype=float)
                    self.model_namespace['valid_value'] = np.array([[value]], dtype=float)
                elif train_validate == 'validate' and type(self.model_namespace['valid_array']) == np.ndarray:
                    self.model_namespace['valid_array'] = np.vstack((self.model_namespace['valid_array'],
                                                                     np.array([countlist], dtype=float)))
                    self.model_namespace['valid_value'] = np.vstack((self.model_namespace['valid_value'],
                                                                     np.array([[value]], dtype=float)))

        # read constants, fragments and holistic descriptors into the Model namespace
        self.model_namespace['fragment_counts'] = np.array((1, 0))
        for line in file_lines:
            if line[0] == '<cons>':
                self.model_namespace['fragment_counts'].resize((1, int(line[1])+1))
                self.model_namespace['fragment_counts'][0, int(line[1])] = float(line[2])
                self.model_namespace['f'+line[1]] = None
            elif line[0] == '<frag>':
                self.model_namespace['fragment_counts'].resize((1, int(line[1])+1))
                self.model_namespace['fragment_counts'][0, int(line[1])] = 0.
                self.model_namespace['f'+line[1]] = ob.OBSmartsPattern()
                self.model_namespace['f'+line[1]].Init(line[2])
            elif line[0] == '<holi>':
                self.model_namespace['fragment_counts'].resize((1, int(line[1])+1))
                self.model_namespace['fragment_counts'][0, int(line[1])] = 0.
                self.model_namespace['f'+line[1]] = line[2]
        
        # read coefficients into the Model namespace
        max_row = 0
        max_col = 0
        for line in file_lines:
            if line[0] == '<coeff>':
                if len(line)-2 > max_row:
                    max_row = len(line)-2
                if int(line[1])+1 > max_col:
                    max_col = int(line[1])+1
        self.model_namespace['coefficients'] = np.zeros((max_row, max_col), dtype=float)
        for line in file_lines:
            if line[0] == '<coeff>':
                for c in range(2, len(line)):
                    self.model_namespace['coefficients'][c-2, int(line[1])] = float(line[c])

        # read fragment standard deviations
        self.model_namespace['fragment_stdev'] = np.zeros((1, 0), dtype=float)
        for line in file_lines:
            if line[0] == '<fstdev>':
                self.model_namespace['fragment_stdev'] = np.hstack((self.model_namespace['fragment_stdev'],
                                                                    np.array([[float(line[2])]], dtype=float)))

        # read atom coverages
        self.model_namespace['atomcoverage'] = []
        for line in file_lines:
            if line[0] == '<atomcoverage>':
                self.model_namespace['atomcoverage'].append([ob.OBSmartsPattern(), ob.OBSmartsPattern(), line[3]])
                self.model_namespace['atomcoverage'][-1][0].Init(line[1])
                self.model_namespace['atomcoverage'][-1][1].Init(line[2])

        # load xtxi matrix for calculating leverages
        for line in file_lines:
            if line[0] == '<xtxi>':
                # try to decode as python 3 string
                try:
                    self.model_namespace['xtxi'] = pickle.loads(base64.b64decode(line[1].encode('utf-8')))
                # if that fails then decode as a python 2 string
                except UnicodeDecodeError:
                    self.model_namespace['xtxi'] = pickle.loads(base64.b64decode(line[1]), encoding='latin1')

        # read commands into the Model namespace converting types
        self.model_namespace['command_sequence'] = []
        for i in range(len(file_lines)):
            line = file_lines[i]
            # variable definitions
            if line[0] == '<define>':
                try:
                    line[2] = float(line[2])
                except ValueError:
                    pass
                self.model_namespace['command_sequence'].append(line)
            # apply coefficients
            elif line[0] == '<applycoeff>':
                line[1] = int(line[1])
                line[2] = int(line[2])
                self.model_namespace['command_sequence'].append(line)
            # apply similarity
            elif line[0] == '<applysimil>':
                line[1] = int(line[1])
                line[2] = int(line[2])
                self.model_namespace['command_sequence'].append(line)
            # sum counts
            elif line[0] == '<sumcounts>':
                line[1] = int(line[1])
                line[2] = int(line[2])
                self.model_namespace['command_sequence'].append(line)
            # calculated CSS
            elif line[0] == '<calculatecss>':
                line[1] = int(line[1])
                line[2] = int(line[2])
                self.model_namespace['command_sequence'].append(line)
            # calculate leverage
            elif line[0] == '<calculateleverage>':
                line[1] = int(line[1])
                line[2] = int(line[2])
                self.model_namespace['command_sequence'].append(line)
            elif line[0] == '<checkatomviolations>':
                self.model_namespace['command_sequence'].append(line)
            # if statement
            elif line[0] == '<if>':
                try:
                    line[1] = float(line[1])
                except ValueError:
                    pass
                try:
                    line[3] = float(line[3])
                except ValueError:
                    pass
                self.model_namespace['command_sequence'].append(line)
            # else statement
            elif line[0] == '<else>':
                self.model_namespace['command_sequence'].append(line)
            # end if statement
            elif line[0] == '<endif>':
                self.model_namespace['command_sequence'].append(line)
            elif line[0] in ['<add>', '<multiply>', '<power>', '<logarithm>', '<round>']:
                try:
                    line[2] = float(line[2])
                except ValueError:
                    pass
                try:
                    line[3] = float(line[3])
                except ValueError:
                    pass
                self.model_namespace['command_sequence'].append(line)
            elif line[0] in ['<concat>']:
                self.model_namespace['command_sequence'].append(line)
            
    def apply_model(self, molecule):
        """Take an openbabel molecule, apply the Model and return the result."""
        # add or delete hydrogens depending on model
        if self.model_namespace['smolecule_format'] == 'old_format':
            molecule.AddHydrogens()
        else:
            molecule.DeleteHydrogens()
        # parse through fragments applying SMARTS or calculating holistic descriptors
        if len(self.model_namespace['fragment_counts'].shape) >= 2:
            ci = self.model_namespace['fragment_counts'].shape[1]
        else:
            ci = 0
        for f in range(ci):
            key = 'f'+str(f)
            if self.model_namespace[key] is None:
                pass
            else:
                self.model_namespace[key].Match(molecule)
                self.model_namespace['fragment_counts'][0, f] = len(self.model_namespace[key].GetUMapList())
        # parse through command sequence
        blocks = [True]
        for command in self.model_namespace['command_sequence']:
            # variable definitions
            if blocks[0] and command[0] == '<define>':
                self.model_namespace[command[1]] = command[2]
            # apply coefficients
            elif blocks[0] and command[0] == '<applycoeff>':
                coeff_count = (self.model_namespace['fragment_counts'][0, command[1]:command[2]] *
                               self.model_namespace['coefficients'][0, command[1]:command[2]]).sum()
                for c in range(1, self.model_namespace['coefficients'].shape[0]):
                    coeff_count += (self.model_namespace['fragment_counts'][0, command[1]:command[2]] *
                                    self.model_namespace['coefficients'][c, command[1]:command[2]]).sum()
                coeff_count /= float(self.model_namespace['coefficients'].shape[0])
                if command[3] == '<addto>':
                    self.model_namespace[command[4]] += coeff_count
                elif command[3] == '<subfrom>':
                    self.model_namespace[command[4]] -= coeff_count
            # sum counts
            elif blocks[0] and command[0] == '<sumcounts>':
                count = (self.model_namespace['fragment_counts'][0, command[1]:command[2]]).sum()
                if command[3] == '<addto>':
                    self.model_namespace[command[4]] += count
                elif command[3] == '<subfrom>':
                    self.model_namespace[command[4]] -= count
            # calculate CSS
            elif blocks[0] and command[0] == '<calculatecss>':
                # compare this chemical to all other chemicals and find the top n most similar chemicals
                topn = 5
                topgroup = []
                mintop = 0.
                for t in range(self.model_namespace['train_array'].shape[0]):
                    fragsim = calculate_fragment_similarity(self.model_namespace['fragment_counts'][0, command[1]:command[2]],
                                                            self.model_namespace['train_array'][t, command[1]:command[2]],
                                                            self.model_namespace['fragment_stdev'][0, command[1]:command[2]],
                                                            mintop)
                    topgroup.append((fragsim, 1-self.model_namespace['train_value_similarity'][t, 0]))
                    topgroup.sort(reverse=True)
                    if len(topgroup) > topn:
                        for i in reversed(range(topn, len(topgroup))):
                            if topgroup[i][0] != topgroup[topn-1][0]:
                                topgroup.pop(i)
                    mintop = topgroup[-1][0]
                # for chemicals with the same similarity scores take the median value similarity,
                # otherwise calculate CSS as normal
                csstype = 'notmedian'
                if csstype == 'median':
                    css = 1
                    findmedian = []
                    n = 0
                    for i in range(len(topgroup)):
                        if topgroup[i][0] != topgroup[topn-1][0]:
                            css *= (topgroup[i][0]*(1-topgroup[i][1]))**0.5
                            n += 1
                        else:
                            findmedian.append(1-topgroup[i][1])
                    findmedian.sort()
                    median = findmedian[len(findmedian)//2]
                    for i in range(topn-n):
                        css *= (topgroup[topn-1][0]*median)**0.5
                # old method of calculating CSS; if there is a tie in fragment similarity use the
                # chemical with the lowest (worst) value similarity
                else:
                    css = 1
                    for i in range(topn):
                        css *= (topgroup[i][0]*(1-topgroup[i][1]))**0.5
                css = css**(1/float(topn))
                if command[3] == '<addto>':
                    self.model_namespace[command[4]] += css
                elif command[3] == '<subfrom>':
                    self.model_namespace[command[4]] -= css
            # calculate leverage
            elif blocks[0] and command[0] == '<calculateleverage>':
                x = np.mat(self.model_namespace['fragment_counts'][0, command[1]:command[2]])
                leverage = (x * self.model_namespace['xtxi'] * x.T)[0, 0]
                # print('leverage', leverage)
                if command[3] == '<addto>':
                    self.model_namespace[command[4]] += leverage
                elif command[3] == '<subfrom>':
                    self.model_namespace[command[4]] -= leverage
            # check atom coverage to see if prediction is valid
            elif blocks[0] and command[0] == '<checkatomviolations>':
                violations = []
                for atomcheck in self.model_namespace['atomcoverage']:
                    atomcheck[0].Match(molecule)
                    atomcheck[1].Match(molecule)
                    if len(atomcheck[0].GetUMapList()) != len(atomcheck[1].GetUMapList()):
                        violations.append(atomcheck[2])
                if len(violations):
                    outline = ''
                    for v in violations:
                        outline += v+', '
                    self.model_namespace[command[1]] = 1
                    self.model_namespace[command[2]] = outline[:-2]
                else:
                    self.model_namespace[command[1]] = 0
                    self.model_namespace[command[2]] = ''
            # if statement
            elif command[0] == '<if>':
                if command[1] in self.model_namespace:
                    left = self.model_namespace[command[1]]
                else:
                    left = command[1]
                if command[3] in self.model_namespace:
                    right = self.model_namespace[command[3]]
                else:
                    right = command[3]
                if command[2] == '<greaterthan>' and left > right:
                    blocks.insert(0, True and blocks[0])
                elif command[2] == '<greaterthan>' and left <= right:
                    blocks.insert(0, False and blocks[0])
                if command[2] == '<lessthan>' and left < right:
                    blocks.insert(0, True and blocks[0])
                elif command[2] == '<lessthan>' and left >= right:
                    blocks.insert(0, False and blocks[0])
                if command[2] == '<equalto>' and left == right:
                    blocks.insert(0, True and blocks[0])
                elif command[2] == '<equalto>' and left != right:
                    blocks.insert(0, False and blocks[0])
            # else statement
            elif command[0] == '<else>':
                blocks[0] = (not blocks[0]) and blocks[1]
            # end if statement
            elif command[0] == '<endif>':
                blocks.pop(0)
            # math functions with left and right side inputs
            elif blocks[0] and command[0] == '<add>':
                if command[2] in self.model_namespace:
                    left = self.model_namespace[command[2]]
                else:
                    left = command[2]
                if command[3] in self.model_namespace:
                    right = self.model_namespace[command[3]]
                else:
                    right = command[3]
                self.model_namespace[command[1]] = left+right
            elif blocks[0] and command[0] == '<multiply>':
                if command[2] in self.model_namespace:
                    left = self.model_namespace[command[2]]
                else:
                    left = command[2]
                if command[3] in self.model_namespace:
                    right = self.model_namespace[command[3]]
                else:
                    right = command[3]
                self.model_namespace[command[1]] = left*right
            elif blocks[0] and command[0] == '<power>':
                if command[2] in self.model_namespace:
                    left = self.model_namespace[command[2]]
                else:
                    left = command[2]
                if command[3] in self.model_namespace:
                    right = self.model_namespace[command[3]]
                else:
                    right = command[3]
                self.model_namespace[command[1]] = left**right
            elif blocks[0] and command[0] == '<logarithm>':
                if command[2] in self.model_namespace:
                    left = self.model_namespace[command[2]]
                else:
                    left = command[2]
                if command[3] in self.model_namespace:
                    right = self.model_namespace[command[3]]
                else:
                    right = command[3]
                # left is base
                self.model_namespace[command[1]] = np.log(right)/np.log(left)
            elif blocks[0] and command[0] == '<round>':
                if command[2] in self.model_namespace:
                    left = self.model_namespace[command[2]]
                else:
                    left = command[2]
                if command[3] in self.model_namespace:
                    right = self.model_namespace[command[3]]
                else:
                    right = command[3]
                self.model_namespace[command[1]] = round(left, int(right))
            elif blocks[0] and command[0] == '<concat>':
                if command[2] in self.model_namespace:
                    left = str(self.model_namespace[command[2]])
                else:
                    left = command[2]
                if command[3] in self.model_namespace:
                    right = str(self.model_namespace[command[3]])
                else:
                    right = command[3]
                self.model_namespace[command[1]] = left+right
            # math functions with only one input
            elif blocks[0] and command[0] == '<ln>':
                if command[2] in self.model_namespace:
                    right = self.model_namespace[command[2]]
                else:
                    right = command[2]
                self.model_namespace[command[1]] = np.log(right)

        if 'WARN' in self.model_namespace and 'ERROR' in self.model_namespace:
            return self.model_namespace['RETURN'], int(self.model_namespace['WARN']), \
                   self.model_namespace['ERROR'], self.model_namespace['NOTE']
        else:
            return self.model_namespace['RETURN'], None, None, None


def apply_model_to_file(model, filename, outfilename=False):
    """apply_model_to_file(model,filename)
    -take a model object and apply it to smiles in a file
    -output a new file with the results
    v0.0.3 - original coding"""

    ##open file
    try:
        infile = open(filename, 'r')
    except IOError:
        print('File not found:', filename)
        return

    ##read fields from header file then parse file contents
    smiles = []
    data = []
    header = False
    for line in infile:
        columns = line.rstrip('\n').split('\t')

        ##read column header
        if not header:
            try:
                assert 'smiles' in columns
            except:
                print('"smiles" missing from column header of file!')
                infile.close()
                return
            header = columns
            continue

        ##add new chemical
        smiles.append(columns[header.index('smiles')])
        data.append(line.rstrip('\n'))

    ##close infile
    infile.close()

    ##apply model to the smiles
    obconversion = ob.OBConversion()
    obconversion.SetInAndOutFormats('smi', 'can')

    results = []
    for s in smiles:
        # instantiate OBMol
        mol = ob.OBMol()
        # read smiles into openbabel
        obconversion.ReadString(mol, s)
        results.append(model.apply_model(mol))

    ##output results to file
    if outfilename:
        outfile = open(outfilename, 'w')
    else:
        outfile = open(filename.replace('.txt', '') + '_' + model.model_namespace['svalue_name'] + '.txt', 'w')
    headerline = ''
    for h in header:
        headerline += h + '\t'
    outfile.write(headerline + 'prediction\tUL\terror\tnote\n')
    for i in range(len(data)):
        outfile.write(data[i] + '\t' + str(results[i][0]) + '\t' + str(results[i][1]) + '\t' + str(results[i][2]) + '\t' + str(results[i][3]) + '\n')
    outfile.close()

if __name__ == "__main__":
    m = Model('./ifs_qsar_dsm_linr.txt')
    apply_model_to_file(m, './dataset_ws.txt')
