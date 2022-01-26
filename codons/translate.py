# -*- coding: utf-8 -*-
"""
Created on Wed Jan 26 00:33:29 2022

@author: Andrew Freiburger
"""
from pprint import pprint
import json, os, re

# allows case insensitive dictionary searches
class CaseInsensitiveDict(dict):   # sourced from https://stackoverflow.com/questions/2082152/case-insensitive-dictionary
    @classmethod
    def _k(cls, key):
        return key.lower() if isinstance(key, str) else key

    def __init__(self, *args, **kwargs):
        super(CaseInsensitiveDict, self).__init__(*args, **kwargs)
        self._convert_keys()
        
    def __getitem__(self, key):
        return super(CaseInsensitiveDict, self).__getitem__(self.__class__._k(key))
    
    def __setitem__(self, key, value):
        super(CaseInsensitiveDict, self).__setitem__(self.__class__._k(key), value)
        
    def __delitem__(self, key):
        return super(CaseInsensitiveDict, self).__delitem__(self.__class__._k(key))
    
    def __contains__(self, key):
        return super(CaseInsensitiveDict, self).__contains__(self.__class__._k(key))
    
    def has_key(self, key):
        return super(CaseInsensitiveDict, self).has_key(self.__class__._k(key))
    
    def pop(self, key, *args, **kwargs):
        return super(CaseInsensitiveDict, self).pop(self.__class__._k(key), *args, **kwargs)
    
    def get(self, key, *args, **kwargs):
        return super(CaseInsensitiveDict, self).get(self.__class__._k(key), *args, **kwargs)
    
    def setdefault(self, key, *args, **kwargs):
        return super(CaseInsensitiveDict, self).setdefault(self.__class__._k(key), *args, **kwargs)
    
    def update(self, E=None, **F):
        super(CaseInsensitiveDict, self).update(self.__class__(E))
        super(CaseInsensitiveDict, self).update(self.__class__(**F))
        
    def _convert_keys(self):
        for k in list(self.keys()):
            v = super(CaseInsensitiveDict, self).pop(k)
            self.__setitem__(k, v)
            
            

class Codons():
    def __init__(self,
                 codons_table: str = 'standard', # the translation table for codons to amino acids
                 amino_acids_form: str = 'full_name', # selects the scale of amino acid nomenclature
                 hyphenated: bool = None # selects whether the printed protein will be hyphenated between the protein residues
                 ):
        self.proteins = []
        
        # define the simulation paths
        self.paths = {}
        self.paths['changed_codons'] = os.path.join(os.path.dirname(__file__), 'rosetta_stone', 'changed_codons.json')
        self.paths['standard_table'] = os.path.join(os.path.dirname(__file__), 'rosetta_stone', 'standard_table.json')
        self.paths['amino_acid_synonyms'] = os.path.join(os.path.dirname(__file__), 'rosetta_stone', 'amino_acid_synonyms.json')
        
        self.parameters = {}
        self.parameters['residue_delimiter'] = '-' 
        
        # define the proper codons table
        self.codons_table = json.load(open(self.paths['standard_table']))
        if codons_table != 'standard':
            self._convert_codon_tables(codons_table)    
        self.codons_table = CaseInsensitiveDict(self.codons_table)
        
        # define the amino acid nomenclature
        if amino_acids_form != 'full_name':
            self.amino_acid_synonyms = json.load(open(self.paths['amino_acid_synonyms']))
            for codon in self.codons_table:
                amino_acid = self.codons_table[codon] 
                if amino_acid != 'stop':
                    self.codons_table[codon] = self.amino_acid_synonyms[amino_acid][amino_acids_form]
            
            if amino_acids_form == 'one_letter' and not hyphenated:
                self.parameters['residue_delimiter'] = ''
                
    
    def _convert_codon_tables(self,codons_table):
        # convert the standard table into the desired table
        self.changed_codons = json.load(open(self.paths['changed_codons']))
        if codons_table not in self.changed_codons:
            raise IndexError(f'The {codons_table} parameter is not presently supported by the options: {list(self.changed_codons.keys())}. Edit the < changed_codons.json > file to offer the exchanges that you desire for your simulation.')
            
        self.changed_codons = self.changed_codons[codons_table]
        for cd in self.changed_codons:
            self.codons_table[cd] = self.changed_codons[cd]
        
        
    def translation(self,
                       sequence: str # the genomic code as a string
                       ):
        codon = ''
        amino_acids = []
        for nuc in sequence:
            if re.search('[atucg]', nuc, flags = re.IGNORECASE):
                codon += nuc
                if len(codon) == 3:
                    amino_acid = self.codons_table[codon]
                    codon = ''
                    if amino_acid is None:
                        continue
                    elif amino_acid == 'stop':
                        
                        protein = self.parameters['residue_delimiter'].join(amino_acids)
                        self.proteins.append(protein)
                        continue
                    amino_acids.append(amino_acid)
                