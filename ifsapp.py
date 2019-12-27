"""ifsapp.py
Implements a simple GUI for applying group contribution QSARS (IFS) developed by Trevor N. Brown.
"""

import openbabel as ob
import ifs_model_read
import smiles_norm
import tkinter as tk
from tkinter import filedialog
from idlelib.redirector import WidgetRedirector
import sys
import os


class _ReadOnlyText(tk.Text):
    """Subclass of tk.Text that is read-only."""

    def __init__(self, *args, **kwargs):
        """Replace insert and delete bindings."""
        # subclass to tk.Text
        tk.Text.__init__(self, *args, **kwargs)
        self.redirector = WidgetRedirector(self)
        # freeze user changes
        self.insert = self.redirector.register('insert', lambda *args, **kw: 'break')
        self.delete = self.redirector.register('delete', lambda *args, **kw: 'break')
        # bind ctrl-a as select all
        self.bind("<Control-Key-a>", self.select_all)

    def select_all(self, event):
        """Select all event bound to ctrl-a."""
        # select all text
        self.tag_add(tk.SEL, '1.0', tk.END)
        self.mark_set(tk.INSERT, '1.0')
        self.see(tk.INSERT)
        return 'break'


class AppManagerClass:
    """A GUI interface for reading in structures as SMILES and applying
    QSARs to the structures."""

    def __init__(self):
        """GUI enters single mode by default. Available QSARs are loaded."""
        # create frame
        self.frame = tk.Frame(root)
        self.frame.pack_propagate(0)
        self.frame.pack()
        # load models
        if hasattr(sys, '_MEIPASS'):
            self.models = [ifs_model_read.Model(os.path.join(sys._MEIPASS, 'ifs_qsar_fhlb_linr.txt')),
                           ifs_model_read.Model(os.path.join(sys._MEIPASS, 'ifs_qsar_hhlb_linr.txt')),
                           ifs_model_read.Model(os.path.join(sys._MEIPASS, 'ifs_qsar_hhlt_linr.txt')),
                           ifs_model_read.Model(os.path.join(sys._MEIPASS, 'ifs_qsar_dsm_linr.txt')),
                           ifs_model_read.Model(os.path.join(sys._MEIPASS, 'ifs_qsar_tm_linr.txt')),
                           ifs_model_read.Model(os.path.join(sys._MEIPASS, 'ifs_qsar_ADB_UFZ__E_linr.txt')),
                           ifs_model_read.Model(os.path.join(sys._MEIPASS, 'ifs_qsar_ADB_UFZ__S_linr.txt')),
                           ifs_model_read.Model(os.path.join(sys._MEIPASS, 'ifs_qsar_ADB_UFZ__A_linr.txt')),
                           ifs_model_read.Model(os.path.join(sys._MEIPASS, 'ifs_qsar_ADB_UFZ__B_linr.txt')),
                           ifs_model_read.Model(os.path.join(sys._MEIPASS, 'ifs_qsar_ADB_UFZ__L_linr.txt')),
                           ifs_model_read.Model(os.path.join(sys._MEIPASS, 'ifs_qsar_V.txt')),
                           ]
        else:
            self.models = [ifs_model_read.Model('ifs_qsar_fhlb_linr.txt'),
                           ifs_model_read.Model('ifs_qsar_hhlb_linr.txt'),
                           ifs_model_read.Model('ifs_qsar_hhlt_linr.txt'),
                           ifs_model_read.Model('ifs_qsar_dsm_linr.txt'),
                           ifs_model_read.Model('ifs_qsar_tm_linr.txt'),
                           ifs_model_read.Model('ifs_qsar_ADB_UFZ__E_linr.txt'),
                           ifs_model_read.Model('ifs_qsar_ADB_UFZ__S_linr.txt'),
                           ifs_model_read.Model('ifs_qsar_ADB_UFZ__A_linr.txt'),
                           ifs_model_read.Model('ifs_qsar_ADB_UFZ__B_linr.txt'),
                           ifs_model_read.Model('ifs_qsar_ADB_UFZ__L_linr.txt'),
                           ifs_model_read.Model('ifs_qsar_V.txt'),
                           ]
        # setup openbabel converter
        self.obcon = ob.OBConversion()
        self.obcon.SetInAndOutFormats('smi', 'can')
        # set to single mode
        self.setup_single_mode()

    def setup_single_mode(self):
        """Delete any batch mode widgets present and load single mode widgets."""
        # delete batch mode widgets
        if hasattr(self, 'buttoninput'):
            self.buttoninput.destroy()
            delattr(self, 'buttoninput')
        if hasattr(self, 'textinput'):
            self.textinput.destroy()
            delattr(self, 'textinput')
        if hasattr(self, 'buttonoutput'):
            self.buttonoutput.destroy()
            delattr(self, 'buttonoutput')
        if hasattr(self, 'textoutput'):
            self.textoutput.destroy()
            delattr(self, 'textoutput')
        if hasattr(self, 'buttongotosingle'):
            self.buttongotosingle.destroy()
            delattr(self, 'buttongotosingle')
        if hasattr(self, 'buttoninfo'):
            self.buttoninfo.destroy()
            delattr(self, 'buttoninfo')
        if hasattr(self, 'framebatch'):
            self.framebatch.destroy()
            delattr(self, 'framebatch')
        if hasattr(self, 'buttoncalcbatch'):
            self.buttoncalcbatch.destroy()
            delattr(self, 'buttoncalcbatch')
        # single mode - text box to enter smiles
        self.labelsmiles = tk.Label(self.frame, text='Enter a SMILES', font='TkDefaultFont 10')
        self.labelsmiles.grid(row=0, column=0)
        self.entrysmiles = tk.Entry(self.frame, width=50, font='TkTextFont 10')
        self.entrysmiles.grid(row=0, column=1)
        # single mode - text box to show results
        self.labelresult = tk.Label(self.frame, text='Model Results', font='TkDefaultFont 10')
        self.labelresult.grid(row=1, column=0)
        self.frametext = tk.Frame(self.frame)
        self.frametext.grid(row=1, column=1)
        self.scrbar = tk.Scrollbar(self.frametext)
        self.scrbar.pack(side=tk.RIGHT, fill=tk.Y)
        self.textresult = _ReadOnlyText(self.frametext, height=5, width=48, font='TkTextFont 10')
        self.textresult.pack(side=tk.LEFT)
        self.scrbar.config(command=self.textresult.yview)
        self.textresult.config(yscrollcommand=self.scrbar.set)
        # single mode - buttons to switch to batch mode and display info
        self.framesingle = tk.Frame(self.frame)
        self.framesingle.grid(row=2, column=0)
        self.buttongotobatch = tk.Button(self.framesingle, text='Batch Mode', font='TkDefaultFont 10',
                                         height=1, width=15, command=self.setup_batch_mode)
        self.buttongotobatch.grid(row=0, column=0)
        self.buttoninfo = tk.Button(self.framesingle, text='Info', font='TkDefaultFont 10',
                                    height=1, width=15, command=self.info)
        self.buttoninfo.grid(row=1, column=0)
        # single mode - button to calculate results
        self.buttoncalcsingle = tk.Button(self.frame, text='Apply IFS QSARs', font='TkDefaultFont 10',
                                          height=2, width=25, command=self.calculate_single)
        self.buttoncalcsingle.grid(row=2, column=1)
        
    def setup_batch_mode(self):
        """Delete any single mode widgets present and load batch mode widgets."""
        # delete single mode widgets
        if hasattr(self, 'labelsmiles'):
            self.labelsmiles.destroy()
            delattr(self, 'labelsmiles')
        if hasattr(self, 'entrysmiles'):
            self.entrysmiles.destroy()
            delattr(self, 'entrysmiles')
        if hasattr(self, 'labelresult'):
            self.labelresult.destroy()
            delattr(self, 'labelresult')
        if hasattr(self, 'textresult'):
            self.textresult.destroy()
            delattr(self, 'textresult')
        if hasattr(self, 'scrbar'):
            self.scrbar.destroy()
            delattr(self, 'scrbar')
        if hasattr(self, 'frametext'):
            self.frametext.destroy()
            delattr(self, 'frametext')
        if hasattr(self, 'buttongotobatch'):
            self.buttongotobatch.destroy()
            delattr(self, 'buttongotobatch')
        if hasattr(self, 'buttoninfo'):
            self.buttoninfo.destroy()
            delattr(self, 'buttoninfo')
        if hasattr(self, 'framesingle'):
            self.framesingle.destroy()
            delattr(self, 'framesingle')
        if hasattr(self, 'buttoncalcsingle'):
            self.buttoncalcsingle.destroy()
            delattr(self, 'buttoncalcsingle')
        # clear filenames if they exist
        self.inputfilename = ''
        self.outputfilename = ''
        # batch mode - button to select input file
        self.buttoninput = tk.Button(self.frame, text='Select Input File', font='TkDefaultFont 10',
                                     height=1, width=15, command=self.select_input_file)
        self.buttoninput.grid(row=0, column=0)
        self.textinput = _ReadOnlyText(self.frame, height=1, width=50, font='TkTextFont 10')
        self.textinput.grid(row=0, column=1)
        # batch mode - button to select output file
        self.buttonoutput = tk.Button(self.frame, text='Select Ouput File', font='TkDefaultFont 10',
                                      height=1, width=15, command=self.select_output_file)
        self.buttonoutput.grid(row=1, column=0)
        self.textoutput = _ReadOnlyText(self.frame, height=1, width=50, font='TkTextFont 10')
        self.textoutput.grid(row=1, column=1)
        # batch mode - buttons to switch to single mode and display info
        self.framebatch = tk.Frame(self.frame)
        self.framebatch.grid(row=2, column=0)
        self.buttongotosingle = tk.Button(self.framebatch, text='Single Mode', font='TkDefaultFont 10',
                                          height=1, width=15, command=self.setup_single_mode)
        self.buttongotosingle.grid(row=0, column=0)
        self.buttoninfo = tk.Button(self.framebatch, text='Info', font='TkDefaultFont 10',
                                    height=1, width=15, command=self.info)
        self.buttoninfo.grid(row=1, column=0)
        # batch mode - button to calculate results
        self.buttoncalcbatch = tk.Button(self.frame, text='Apply IFS QSARs', font='TkDefaultFont 10',
                                         height=2, width=25, command=self.calculate_batch)
        self.buttoncalcbatch.grid(row=2, column=1)

    def toggle_disabled(self, setstate):
        """Disable or enable all widgets."""
        if hasattr(self, 'buttoninput'):
            self.buttoninput.config(state=setstate)
            self.buttoninput.update()
        if hasattr(self, 'textinput'):
            self.textinput.config(state=setstate)
            self.textinput.update()
        if hasattr(self, 'buttonoutput'):
            self.buttonoutput.config(state=setstate)
            self.buttonoutput.update()
        if hasattr(self, 'textoutput'):
            self.textoutput.config(state=setstate)
            self.textoutput.update()
        if hasattr(self, 'buttongotosingle'):
            self.buttongotosingle.config(state=setstate)
            self.buttongotosingle.update()
        if hasattr(self, 'buttoninfo'):
            self.buttoninfo.config(state=setstate)
            self.buttoninfo.update()
        if hasattr(self, 'buttoncalcbatch'):
            self.buttoncalcbatch.config(state=setstate)
            self.buttoncalcbatch.update()
        if hasattr(self, 'labelsmiles'):
            self.labelsmiles.config(state=setstate)
            self.labelsmiles.update()
        if hasattr(self, 'entrysmiles'):
            self.entrysmiles.config(state=setstate)
            self.entrysmiles.update()
        if hasattr(self, 'labelresult'):
            self.labelresult.config(state=setstate)
            self.labelresult.update()
        if hasattr(self, 'buttongotobatch'):
            self.buttongotobatch.config(state=setstate)
            self.buttongotobatch.update()
        if hasattr(self, 'buttoninfo'):
            self.buttoninfo.config(state=setstate)
            self.buttoninfo.update()
        if hasattr(self, 'buttoncalcsingle'):
            self.buttoncalcsingle.config(state=setstate)
            self.buttoncalcsingle.update()
        
    def info(self):
        """Display the module README.md in a popup window."""
        # create a popup to show information about the program
        popup = tk.Tk()
        popup.wm_title('Info')
        scrbar = tk.Scrollbar(popup)
        scrbar.pack(side=tk.RIGHT, fill=tk.Y)
        text = _ReadOnlyText(popup, font='Consolas 10')
        with open('README.md', 'r') as file:
            infostring = file.read()
        infostring = infostring.replace('<pre>', '').replace('</pre>', '').replace('\[', '[')
        text.insert('1.0', infostring)
        text.pack(side=tk.LEFT)
        scrbar.config(command=text.yview)
        text.config(yscrollcommand=scrbar.set)
        popup.mainloop()

    def calculate_single(self):
        """Get the SMILES from the GUI and apply QSARs."""
        # get smiles from the gui, apply models and write results to gui
        smiles = self.entrysmiles.get()
        self.toggle_disabled(tk.DISABLED)
        # convert smiles to standardized obmol
        molecule, newsmiles, conversionnote = smiles_norm.convert(smiles, obconversion=self.obcon)
        # check conversion results and output smiles and conversion note
        result = 'input SMILES: ' + smiles + '\n'
        conversionsuccess = True
        if 'error reading SMILES' in conversionnote:
            conversionsuccess = False
            result += 'normalized SMILES: \nSMILES conversion note: ' + conversionnote + '; check input\n'
        elif 'aromaticity broken' in conversionnote:
            conversionsuccess = False
            result += 'normalized SMILES: \nSMILES conversion note: ' + conversionnote + \
                      '; structure or conversion error\n'
        elif 'structure contains permanently charged atoms' in conversionnote:
            conversionsuccess = False
            result += 'normalized SMILES: \nSMILES conversion note: ' + conversionnote + \
                      '; QSARs only handle neutrals\n'
        elif 'charged atom(s) neutralized' in conversionnote:
            result += 'normalized SMILES: ' + newsmiles + '\nSMILES conversion note: ' + conversionnote + \
                      '; QSARs only handle neutrals\n'
        else:
            result += 'normalized SMILES: ' + newsmiles + '\nSMILES conversion note: ' + conversionnote + '\n'
        # prepare model outputs
        if conversionsuccess:
            result += '\nproperty (units)\tprediction\terror\tUL\tnote\n'
            for m in self.models:
                pred, warn, error, note = m.apply_model(molecule)
                result += m.model_namespace['svalue_name'] + ' ('+m.model_namespace['sunits'] + ')' + '\t' + \
                    str(pred) + '\t' + str(error) + '\t' + str(warn) + '\t' + str(note) + '\n'
        # display results
        self.textresult.delete('1.0', tk.END)
        self.textresult.insert('1.0', str(result))
        self.toggle_disabled(tk.NORMAL)
        
    def calculate_batch(self):
        """Load SMILES from input file, apply QSARs, and write to output file."""
        # read in a file and apply models to each smiles, then write results to the output file
        if self.inputfilename == '' or self.outputfilename == '':
            return
        self.toggle_disabled(tk.DISABLED)
        infile = open(self.inputfilename, 'r')
        outfile = open(self.outputfilename, 'w')
        header = 0
        smiles_index = 0
        sep = '\t'
        if self.inputfilename[-4:] == '.csv':
            sep = ','
        for line in infile:
            strip = line.rstrip('\n').split(sep)
            if header == 0:
                header = 1
                for s in range(len(strip)):
                    if strip[s].lower().strip().rstrip() == 'smiles':
                        smiles_index = s
                        header = 2
                        break
                if header == 2:
                    outline = ''
                    for field in strip:
                        outline += field + sep
                    result = 'normalized SMILES' + sep + 'SMILES conversion note' + sep
                    for m in self.models:
                        result += m.model_namespace['svalue_name'] + ' (' + m.model_namespace['sunits'] + ')' + \
                                  sep + 'error' + sep + 'UL' + sep + 'note' + sep
                    outfile.write(outline+result[:-1] + '\n')
                    continue
                else:
                    result = 'smiles' + sep + 'normalized SMILES' + sep + 'SMILES conversion note' + sep
                    for m in self.models:
                        result += m.model_namespace['svalue_name'] + ' (' + m.model_namespace['sunits'] + ')' + \
                                  sep + 'error' + sep + 'UL' + sep + 'note' + sep
                    outfile.write(result[:-1] + '\n')
            outline = ''
            for field in strip:
                outline += field + sep
            smiles = strip[smiles_index]
            # convert smiles to standardized obmol
            molecule, newsmiles, conversionnote = smiles_norm.convert(smiles, obconversion=self.obcon)
            # check conversion results and output smiles and conversion note
            result = newsmiles + sep + conversionnote
            conversionsuccess = True
            if 'error reading SMILES' in conversionnote:
                conversionsuccess = False
                result += '; check input' + sep
            elif 'aromaticity broken' in conversionnote:
                conversionsuccess = False
                result += '; structure or conversion error' + sep
            elif 'structure contains permanently charged atoms' in conversionnote:
                conversionsuccess = False
                result += '; QSARs only handle neutrals' + sep
            elif 'charged atom(s) neutralized' in conversionnote:
                result += '; QSARs only handle neutrals' + sep
            else:
                result += sep
            if conversionsuccess:
                for m in self.models:
                    pred, warn, error, note = m.apply_model(molecule)
                    result += str(pred) + sep + str(error) + sep + str(warn) + sep + str(note) + sep
            else:
                for m in range(len(self.models)):
                    result += 'Prediction not possible' + sep + 'None' + sep + 'None' + sep + 'None' + sep
            outfile.write(outline+result[:-1] + '\n')
        self.toggle_disabled(tk.NORMAL)
        self.textinput.delete('1.0', tk.END)
        self.textoutput.delete('1.0', tk.END)

    def select_input_file(self):
        """Spawn a load file popup."""
        # create popup to select output file
        self.inputfilename = filedialog.askopenfilename(title='Select Input File',
                                                        filetypes=[('tab-delimited txt files', '*.txt'),
                                                                   ('csv files', '*.csv')])
        self.textinput.delete('1.0', tk.END)
        self.textinput.insert('1.0', self.inputfilename)

    def select_output_file(self):
        """Spawn a save file popup."""
        # create popup to select input file
        self.outputfilename = filedialog.asksaveasfilename(title='Select Output File',
                                                           defaultextension='.txt',
                                                           filetypes=[('tab-delimited txt files', '*.txt')])
        self.textoutput.delete('1.0', tk.END)
        self.textoutput.insert('1.0', self.outputfilename)

    def update(self):
        """Check for events."""
        # update
        root.after(50, app_manager.update)


# main tk gui loop
if __name__ == "__main__":
    root = tk.Tk()
    root.wm_title('IFSApp')
    app_manager = AppManagerClass()
    root.after(50, app_manager.update)
    root.mainloop()
