'''
Created on 14.02.2015

@author: David
'''

from tkinter import Tk, Frame, LabelFrame, Label, Button, Listbox, Scrollbar,\
    Toplevel, Entry
from tkinter.scrolledtext import ScrolledText
from tkinter.ttk import Combobox
from tkinter import W, E, N, S

import tkinter.filedialog

from tkinter.messagebox import showinfo
import sequenceIO
import colorsys

from Bio.Restriction.Restriction_Dictionary import rest_dict
import tkinter
from test.test_importlib import source

from SeqUtils import stripCharsNotInList
from fileinput import filename


class SpeciesSearchDialog():
    
    def __init__(self, parent, caller):
        # main window
        self.parent=parent
        # dialog that called this second dialog
        self.caller=caller
        self.gui = Toplevel(parent.guiRoot)
        self.gui.grab_set()
        self.gui.focus()
        
        self.gui.columnconfigure(0, weight=1)
        self.gui.rowconfigure(1, weight=1)
        
        self.entrySearch = Entry(self.gui)
        
        self.buttonSearch = Button(self.gui, text=" Search ")
        self.buttonAdd = Button(self.gui, text=" Add Species ")
        self.buttonClose = Button(self.gui, text=" Close Window ")
        
        self.frameResults = Frame(self.gui)
        self.frameResults.columnconfigure(0, weight=1)
        self.frameResults.rowconfigure(0, weight=1)
        self.scrollResults = Scrollbar(self.frameResults, orient="vertical")
        self.scrollResults.grid(row=0, column=1, sticky="ns")
        self.listResults = Listbox(self.frameResults, width=70, height=20)
        self.listResults.grid(row=0, column=0, sticky="nswe")
        
        self.listResults.config(yscrollcommand=self.scrollResults.set)
        self.scrollResults.config(command=self.listResults.yview)
        
        
        self.entrySearch.grid(row=0, column=0, columnspan=2, sticky="we", padx=5, pady=5)
        self.frameResults.grid(row=1, column=0, columnspan=3, sticky="nswe", padx=5, pady=5)
        self.buttonSearch.grid(row=0, column=2, padx=5, pady=5, sticky="e")
        self.buttonAdd.grid(row=2, column=1, padx=5, pady=5, sticky="e")
        self.buttonClose.grid(row=2, column=2, padx=5, pady=5, sticky="e")
        
        self.gui.protocol("WM_DELETE_WINDOW", self.actionClose)        
        self.buttonClose.bind("<ButtonRelease>", self.actionClose)
        self.buttonAdd.bind("<ButtonRelease>", self.actionAdd)
        self.buttonSearch.bind("<ButtonRelease>", self.actionSearch)
        self.entrySearch.bind("<Return>", self.actionSearch)
        
        self.gui.update()
        self.gui.minsize(self.gui.winfo_width(), self.gui.winfo_height())
        
        self.gui.mainloop()
        
    def actionAdd(self, event):
        try:
            selection = self.listResults.selection_get()
            selectionSplit = selection.split(":\t", 1)
            selectionSplit2 = selectionSplit[1].split("\t")
            if not (selectionSplit[0], selectionSplit2[0]) in self.parent.optimizer.speciesList:   
                self.parent.optimizer.speciesList.append((selectionSplit[0], selectionSplit2[0].strip()))
            self.caller.gui.event_generate("<<Update>>")
        except tkinter.TclError :
            # no selection
            pass
    
    def actionSearch(self, event):
        query = self.entrySearch.get()
        if(query):
            
            self.listResults.delete(0, "end")
            
            results = self.parent.optimizer.SPSUMHandler.search(query)
            # sort results by nr of CDS
            results = sorted(results.items(), reverse=True, key=lambda x: (int(x[1][1])))
            for name, (taxid, ncs) in results:
                self.listResults.insert("end", taxid + ":\t " + name + " \t(" + ncs + " CDS)")
        
    def actionClose(self, event=None):
        self.caller.gui.event_generate("<<Update>>", when="tail")
        self.gui.destroy()

class SpeciesListDialog():
    
    def __init__(self, parent):
        self.parent = parent
        self.gui = Toplevel(parent.guiRoot)
        self.gui.grab_set()
        self.gui.focus()
        
        self.gui.columnconfigure(0, weight=1)
        self.gui.rowconfigure(1, weight=1)
        
        Label(self.gui, text="Registered Species:").grid(row=0, column=0, pady=5, padx=5, sticky="w")
        self.listRegisteredSpecies = Listbox(self.gui, width=70)
        self.buttonAdd = Button(self.gui, text=" + ")
        self.buttonDel = Button(self.gui, text=" - ")
        self.listRegisteredSpecies.grid(row=1, column=0, columnspan=3, sticky="nswe", pady=5, padx=5)
        self.buttonAdd.grid(row=2, column=1, pady=5, padx=5)
        self.buttonDel.grid(row=2, column=2, pady=5, padx=5)
        
        # Set (minimum + max) Window size 
        self.gui.update()
        self.gui.minsize(self.gui.winfo_width(), self.gui.winfo_height())
#         self.gui.maxsize(self.gui.winfo_width(), self.gui.winfo_height())
        
        self.actionUpdate(None)
        self.gui.bind("<<Update>>", self.actionUpdate)
        self.gui.protocol("WM_DELETE_WINDOW", self.actionClose)
        
        self.buttonDel.bind("<ButtonRelease>", self.actionDel)
        self.buttonAdd.bind("<ButtonRelease>", self.actionAdd)
        
        
        self.gui.mainloop()
        
    def actionClose(self):
        self.parent.guiRoot.event_generate("<<Update>>", when="tail")
        self.gui.destroy()
        
    def actionUpdate(self, event):
        self.listRegisteredSpecies.delete(0, "end")
        for (taxid, name) in self.parent.optimizer.speciesList:
            self.listRegisteredSpecies.insert("end", taxid + ": " + name)
            
    def actionDel(self, event):
        try:
            selection = self.listRegisteredSpecies.selection_get()
            selectionSplit = selection.split(": ")            
            self.parent.optimizer.speciesList.remove((selectionSplit[0], selectionSplit[1]))
            self.gui.event_generate("<<Update>>")
        except tkinter.TclError :
            # no selection
            pass
        
    def actionAdd(self, Event):
        SpeciesSearchDialog(self.parent, self)
        

class AddRestrictionDialog():
    
    def __init__(self, parent):
        self.parent = parent
        self.gui = Toplevel(parent.guiRoot)
        self.gui.grab_set()
        self.gui.focus()
        
        self.gui.columnconfigure(0, weight=1)
        self.gui.columnconfigure(1, weight=1)
        
        Label(self.gui, text="Enzyme:").grid(row=0, column=0, sticky="e", padx=5)
        self.entryEnzyme = Entry(self.gui)
        self.entryEnzyme.grid(row=0, column=1, sticky="w", padx=5, pady=10)
        self.frameButtonGroup = Frame(self.gui)
        self.frameButtonGroup.grid(row=1, column=0, columnspan=2, pady=10)
        self.buttonOK = Button(self.frameButtonGroup, text=" OK ")
        self.buttonCancel = Button(self.frameButtonGroup, text=" Cancel ")
        self.buttonOK.grid(row=0, column=1, sticky="w", padx=5)
        self.buttonCancel.grid(row=0, column=0, sticky="e", padx=5)
        
        # Set (minimum + max) Window size 
        self.gui.update()
        self.gui.minsize(self.gui.winfo_width(), self.gui.winfo_height())
        self.gui.maxsize(self.gui.winfo_width(), self.gui.winfo_height())
        
        self.entryEnzyme.focus()
        
        self.buttonOK.bind("<ButtonRelease>", self.actionOK)
        self.buttonCancel.bind("<ButtonRelease>", self.actionCancel)
        self.entryEnzyme.bind("<Return>", self.actionOK)
        
        self.gui.mainloop()
#         self.gui.grab_release()
#         self.gui.destroy()

    def actionOK(self, event):
        enzymeString = self.entryEnzyme.get()
        if enzymeString in rest_dict.keys():
            
            if enzymeString in self.parent.optimizer.restrictionEnzymeList:
                showinfo("", (enzymeString + " was already added to the list"))
                return
                
            self.parent.optimizer.restrictionEnzymeList.append(enzymeString)
            self.parent.guiRoot.event_generate("<<Update>>", when="tail")
            self.gui.destroy()
        else:
            showinfo("", (enzymeString + " is not a valid restriction enzyme"))
    
    def actionCancel(self, event):
        self.gui.destroy()
        

class OptimizerMainWindow():
    '''
    classdocs
    '''
    
    #TODO: change that name
    def reactToClick(self, event):
        a = AddRestrictionDialog(self)

    def __init__(self, optimizer):
        
        # always have a reference to model/controller
        self.optimizer = optimizer
        
        # setup main GUI and make stretchable
        self.guiRoot = Tk()
        self.guiRoot.title("OPTIMIZR")
        self.guiRoot.columnconfigure(1, weight=1)
        self.guiRoot.rowconfigure(0, weight=1)
        
        # left (settings) and right (sequences) part
        self.frameLeft = Frame(self.guiRoot)
        self.frameLeft.grid(row=0, column=0, sticky=W+E+N+S)
        self.frameLeft.columnconfigure(0, weight=1)
        self.frameRight = Frame(self.guiRoot)
        self.frameRight.grid(row=0, column=1, sticky=W+E+N+S)
        self.frameRight.columnconfigure(0, weight=1)
        self.frameRight.rowconfigure(0, weight=1)
        self.frameRight.rowconfigure(1, weight=1)
        
        
        self.frameSpeciesControll = LabelFrame(self.frameLeft, text="Species", pady=10, padx=10)
        self.frameSpeciesControll.columnconfigure(1, weight=1)
        self.frameOptimizationControll = LabelFrame(self.frameLeft, text="Optimization", pady=10, padx=10)
        self.frameRestrictionControll = LabelFrame(self.frameLeft, text="Restriction Enzymes", pady=10, padx=10)
        self.frameSpeciesControll.grid(row=0, column=0, sticky=W+E, padx=10, pady=10)
        self.frameOptimizationControll.grid(row=1, column=0, sticky=W+E, padx=10, pady=10)
        self.frameRestrictionControll.grid(row=2, column=0, sticky=W+E, padx=10, pady=10)
        
        # Species Controll
        Label(self.frameSpeciesControll, text="Source:").grid(row=0, column=0)
        Label(self.frameSpeciesControll, text="Target:").grid(row=1, column=0)
        
        self.comboSourceSpecies = Combobox(self.frameSpeciesControll, state="readonly")
        self.comboSourceSpecies.grid(row=0, column=1, pady=5, sticky="ew")
        self.comboTargetSpecies = Combobox(self.frameSpeciesControll, state="readonly")
        self.comboTargetSpecies.grid(row=1, column=1, pady=5, sticky="we")        
        self.buttonSpeciesList = Button(self.frameSpeciesControll, text="Edit Species List")
        self.buttonSpeciesList.grid(row=2, column=1, pady=5, sticky="e")
        
        self.comboSourceSpecies.bind("<<ComboboxSelected>>", self.actionOptimizerSettingsChanged)
        self.comboTargetSpecies.bind("<<ComboboxSelected>>", self.actionOptimizerSettingsChanged)
        
        # Optimization Controll
        Label(self.frameOptimizationControll, text="Optimization Strategy:").grid(row=0, column=0)
        self.comboOptimizationStrategy = Combobox(self.frameOptimizationControll, state="readonly")
        self.comboOptimizationStrategy.grid(row=0, column=1)
        self.comboOptimizationStrategy["values"] = self.optimizer.possibleOptimizationStrategies
        self.comboOptimizationStrategy.bind("<<ComboboxSelected>>", self.actionOptimizerSettingsChanged)
        
        # Restriction Enzymes
        self.listRestriction = Listbox(self.frameRestrictionControll)
        self.listRestriction.grid(row=0, column=0, columnspan=3, pady=5, sticky=W+E)
        self.frameRestrictionControll.columnconfigure(0, weight=1)
        self.buttonRestricionAdd = Button(self.frameRestrictionControll, text=" + ")
        self.buttonRestricionDel = Button(self.frameRestrictionControll, text=" - ")
        self.buttonRestricionAdd.grid(row=1, column=1, padx=5)
        self.buttonRestricionDel.grid(row=1, column=2, padx=5)
        
        # Source Sequence Frame
        self.frameSourceSequence = LabelFrame(self.frameRight, text="Source Sequence", padx=10, pady=10)
        self.frameResultSequence = LabelFrame(self.frameRight, text="Result Sequence", padx=10, pady=10)
        self.frameSourceSequence.grid(row=0, column=0, sticky="wens", padx=10, pady=10)
        self.frameResultSequence.grid(row=1, column=0, sticky="wens", padx=10, pady=10)
        
        
        self.buttonSourceLoad = Button(self.frameSourceSequence, text=" Load ")
        self.textSourceSeq = ScrolledText(self.frameSourceSequence, height=10)
        self.buttonSourceLoad.grid(row=0, column=1, sticky="e", pady=5)
        self.textSourceSeq.grid(row=1, column=0, columnspan=2, sticky="wens")
        self.frameSourceSequence.columnconfigure(0, weight=1)
        self.frameSourceSequence.rowconfigure(1, weight=1)
        self.textSourceSeq.frame.columnconfigure(1, weight=1)
        self.textSourceSeq.frame.rowconfigure(0, weight=1)
        
        self.buttonOptimize = Button(self.frameResultSequence, text=" OPTIMIZE! ")
        self.buttonOptimize.bind("<ButtonRelease>", self.actionOptimize)
        
        self.buttonRemoveRestriction = Button(self.frameResultSequence, text=" RESTRICTION-B-GONE! ")
        self.buttonRemoveRestriction.bind("<ButtonRelease>", self.actionRemoveRestricion)
        
        self.buttonSaveResult = Button(self.frameResultSequence, text=" Save ")
        self.textResultSequence = ScrolledText(self.frameResultSequence, height=10)
        self.buttonOptimize.grid(column=0, row=0, pady=5, sticky="w")
        self.buttonRemoveRestriction.grid(column=1, row=0, pady=5, padx=10, sticky="w")
        self.textResultSequence.grid(row=1, column=0, columnspan=4, sticky="wens")
        self.buttonSaveResult.grid(row=2, column=3, pady=5, sticky="e")
        self.frameResultSequence.columnconfigure(2, weight=1)
        self.frameResultSequence.rowconfigure(1, weight=1)
        self.textResultSequence.frame.columnconfigure(1, weight=1)
        self.textResultSequence.frame.rowconfigure(0, weight=1)
        
        self.textSourceSeq.bind("<<Modified>>", self.actionSequenceModified)
        self.textResultSequence.bind("<<Modified>>", self.actionSequenceModified)
        
        #generate color tags for textboxes
        for i in range(101):
            
            #green for normal codons
            (r,g,b) = colorsys.hsv_to_rgb(210/360, i/100, 1.0)            
            colorHex = "#%02x%02x%02x" % (int(r*255), int(g*255), int(b*255))
            
            self.textSourceSeq.tag_config("normal"+str(i), background=colorHex)
            self.textResultSequence.tag_config("normal"+str(i), background=colorHex)
            
            #red for codons with restriction sites
            (r,g,b) = colorsys.hsv_to_rgb(5/360, i/100, 1.0)            
            colorHex = "#%02x%02x%02x" % (int(r*255), int(g*255), int(b*255))
            
            self.textSourceSeq.tag_config("restrict"+str(i), background=colorHex)
            self.textResultSequence.tag_config("restrict"+str(i), background=colorHex)
        
        
        # Set (minimum + max) Window size 
        self.guiRoot.update()
        self.guiRoot.minsize(self.guiRoot.winfo_width(), self.guiRoot.winfo_height())
        
        self.buttonRestricionAdd.bind("<ButtonRelease>", self.reactToClick)
        self.buttonRestricionDel.bind("<ButtonRelease>", self.actionRestrictionEnzymeDelete)
        self.buttonSpeciesList.bind("<ButtonRelease>", self.actionEditSpeciesButton)
        
        self.buttonSourceLoad.bind("<ButtonRelease>", self.actionLoadSequence)
        self.buttonSaveResult.bind("<ButtonRelease>", self.actionSaveSequence)
        
        # TEST
#         self.listRestriction.insert("end", "EcoRI")
#         self.listRestriction.insert("end", "BamHI")
#         
        
        # dummy event to manually trigger update
        self.guiRoot.bind("<<Update>>", self.actionUpdate)
        
        self.actionUpdate(None)
        
        self.guiRoot.mainloop()
    
    def actionRestrictionEnzymeDelete(self, event):
        try:
            selectedEnzyme = self.listRestriction.selection_get()
            self.optimizer.restrictionEnzymeList.remove(selectedEnzyme)
            self.guiRoot.event_generate("<<Update>>")
        except tkinter.TclError :
            # no selection
            pass    
        
    def actionUpdate(self, event):
#         print("update called")

        # clear list of restriction enzymes
        self.listRestriction.delete(0, "end")        
        for r in self.optimizer.restrictionEnzymeList:
            self.listRestriction.insert("end", r)
            
        self.comboSourceSpecies.delete(0, "end")
        self.comboTargetSpecies.delete(0, "end")
        
        speciesValues = list()        
        for (taxid, name) in self.optimizer.speciesList:
            speciesValues.append(taxid + ": " + name)
            
        self.comboSourceSpecies["values"] = speciesValues
        self.comboTargetSpecies["values"] = speciesValues
        
        if self.comboSourceSpecies.get() not in speciesValues:
            self.comboSourceSpecies.set("")
        if self.comboTargetSpecies.get() not in speciesValues:
            self.comboTargetSpecies.set("")
            
        self.textSourceSeq.edit_modified(True)
        self.textResultSequence.edit_modified(True)
        
        self.optimizer.saveConfig("config.ini")
            
    def actionEditSpeciesButton(self, event):
        speciesListDialog = SpeciesListDialog(self)
        
    def actionOptimizerSettingsChanged(self, event=None):
#         print("Something happened")
        strategy = self.comboOptimizationStrategy.get()
        sourceString = self.comboSourceSpecies.get()
        targetString = self.comboTargetSpecies.get()
        
        if not (strategy and sourceString and targetString):
            return
        
        sourceTaxid = sourceString.split(":")[0]
        targetTaxid = targetString.split(":")[0]
        
        self.optimizer.setOptimizer(sourceTaxid, targetTaxid, strategy)
        
        self.textSourceSeq.edit_modified(True)
        self.textResultSequence.edit_modified(True)
#         self.optimizer.testPrint()

    def actionOptimize(self, event=None):
        self.optimizer.runOptimization()
        self.textSourceSeq.edit_modified(True)
        self.textResultSequence.edit_modified(True)
        
    def actionRemoveRestricion(self, event=None):
        self.optimizer.runRestricionRemoval()
        self.textSourceSeq.edit_modified(True)
        self.textResultSequence.edit_modified(True)
        
    def actionSequenceModified(self, event=None):
        # necessary if, otherwise -> infinite loop
        if self.textSourceSeq.edit_modified():
            
            seq = self.textSourceSeq.get("1.0", "end").strip()
            seq = stripCharsNotInList(seq.upper(), ['A', 'C', 'G', 'T'])
            self.optimizer.setSourceSeq(seq)
            
            oldInsert = self.textSourceSeq.index("insert")
            self.textSourceSeq.delete("1.0", "end")
            
            sourceCodons = self.optimizer.getCodonsForPrint(True)
            if not sourceCodons:
                self.textSourceSeq.insert("end", self.optimizer.sourceSequence)
            else:
                for (co, sc, r) in sourceCodons:
                    if sc:
                        if not r:
                            self.textSourceSeq.insert("end", co, "normal"+str(int(sc*100)))
                            #print("normal"+str(int(sc*100)))
                        else:
                            self.textSourceSeq.insert("end", co, "restrict"+str(int(sc*100)))
                    else:
                        # remainder without color
                        self.textSourceSeq.insert("end", co)
                    
                
            self.textSourceSeq.mark_set("insert", oldInsert)

            # reset the modified status at the very end
            self.textSourceSeq.edit_modified(False)

        if self.textResultSequence.edit_modified():
            
            seq = self.textResultSequence.get("1.0", "end").strip()
#             self.optimizer.setOptimizedSeq(seq)
            
            oldInsert = self.textResultSequence.index("insert")
            self.textResultSequence.delete("1.0", "end")
            
            targetCodons = self.optimizer.getCodonsForPrint(False)
            if not targetCodons:
                self.textSourceSeq.insert("end", self.optimizer.optimizedSequence)
            else:
                for (co, sc, r) in targetCodons:
                    if sc:
                        if not r:
                            self.textResultSequence.insert("end", co, "normal"+str(int(sc*100)))
                            #print("normal"+str(int(sc*100)))
                        else:
                            self.textResultSequence.insert("end", co, "restrict"+str(int(sc*100)))
                    else:
                        # remainder without color
                        self.textResultSequence.insert("end", co)
                    
                
            self.textSourceSeq.mark_set("insert", oldInsert)
            
            self.textResultSequence.edit_modified(False)

    def actionLoadSequence(self, event=None):
        filename = tkinter.filedialog.askopenfilename()
        if filename:
            seq = sequenceIO.readFile(filename)
            self.textSourceSeq.delete("1.0", "end")
            self.textSourceSeq.insert("end", seq)            
            self.textSourceSeq.edit_modified(True)
            
    def actionSaveSequence(self, event=None):
        filename = tkinter.filedialog.asksaveasfilename()
        if filename:
#             print("file is " + filename)
            with open(filename, mode='w') as fd:
                fd.write(self.optimizer.optimizedSequence)
            
            