# -*- coding: utf-8 -*-
"""
Created on Wed Jan 26 00:33:29 2022

@author: Andrew Freiburger
"""
from math import ceil
from Bio.Blast import NCBIWWW
from pprint import pprint
from chemw import Proteins
import datetime
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
                 sequence: str = None,  # the genetic sequence can be optionally provided, for easy use in the other functions.
                 codons_table: str = 'standard', # the translation table for codons to amino acids
                 amino_acids_form: str = 'full_name', # selects the scale of amino acid nomenclature
                 hyphenated: bool = None, # selects whether the printed protein will be hyphenated between the protein residues
                 verbose: bool = False,
                 printing: bool = True
                 ):
        self.verbose = verbose
        self.printing = printing
        self.proteins = {}
        self.fasta = []
        self.transcribed_sequence = None
        self.protein_blast_results = None
        self.nucleotide_blast_results = None
        self.protein_mass = Proteins()
        
        # define the simulation paths
        self.paths = {}
        self.paths['changed_codons'] = os.path.join(os.path.dirname(__file__), 'rosetta_stone', 'changed_codons.json')
        self.paths['standard_table'] = os.path.join(os.path.dirname(__file__), 'rosetta_stone', 'standard_table.json')
        self.paths['amino_acid_synonyms'] = os.path.join(os.path.dirname(__file__), 'rosetta_stone', 'amino_acid_synonyms.json')
        
        self.parameters = {}
        self.parameters['residue_delimiter'] = '-' 
        
        # refine the sequence into the FASTA format
        self.sequence = sequence
        
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
                if amino_acid not in ['stop', 'start']:
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
            
    def make_fasta(self,
                   sequence: str,  # the genetic nucleotide sequence
                   description: str = 'sequence' # the description of the genetic information
                   ):
        if sequence is None:
            return None
        
        return '\n'.join([f'>{description}',sequence])
            
    def transcribe(self,
                   sequence: str = None, # the genomic code as a string
                   description: str = '',  # a description of the sequence
                   ):
        if sequence:
            self.sequence = sequence
            
        # determine the capitalization of the sequence
        for ch in sequence:
            if re.search('[a-zA-Z]',ch): 
                upper_case = ch.isupper()
                break
            
        # substitute the nucleotides with the appropriate capitalization
        transcription = 'Transcribed sequence '
        if re.search('u|U', self.sequence):
            transcription += 'RNA to DNA'
            if upper_case:
                self.transcribed_sequence = re.sub('U', 'T', self.sequence)
            else:
                self.transcribed_sequence = re.sub('u', 't', self.sequence)
        if re.search('t|T', self.sequence):
            transcription += 'DNA to RNA'
            if upper_case:
                self.transcribed_sequence = re.sub('T', 'U', self.sequence)
            else:
                self.transcribed_sequence = re.sub('t', 'u', self.sequence)
                
        print('The sequence is transcribed.')
        self.transcribed_sequence = f'>{transcription}\n{self.transcribed_sequence}*\n'
        return self.transcribed_sequence
        
        
    def translate(self,
                       sequence: str = None, # the genomic code as a string
                       ):
        if sequence:
            self.sequence = sequence
            
        self.multi_fasta = []
        self.missed_codons = []
        codon = ''
        amino_acids = None
        description = 'polypeptide'
        for nuc in self.sequence:
            if re.search('[atucg]', nuc, flags = re.IGNORECASE):
                codon += nuc
                if len(codon) == 3:
                    if codon not in self.codons_table:
                        self.missed_codons.append(codon)
                    else:
                        amino_acid = self.codons_table[codon]
                        if amino_acid == 'start':
                            amino_acids = []
                        elif amino_acid == 'stop' and type(amino_acids) is list:
                            if len(amino_acids) >= 1:
                                protein = self.parameters['residue_delimiter'].join(amino_acids)
                                mass = self.protein_mass.mass(protein)
                                
                                self.proteins[protein] = mass
                                description = ' - '.join(['>Protein', description, f'{len(protein)}residues', f'{mass}amu'])
                                fasta_file = self.make_fasta(description, protein)
                                self.multi_fasta.append(fasta_file)
                                amino_acids = None
                        else:
                            if type(amino_acids) is list and re.search('[a-z]+',amino_acid, flags = re.IGNORECASE):
                                amino_acids.append(amino_acid)
                                if self.verbose:
                                    print(codon, '\t', amino_acid)
                    
                    codon = ''
                    
        self.multi_fasta = '\n'.join(self.multi_fasta)
        if self.printing:
           print(self.multi_fasta)
           if self.missed_codons != []:
               print(f'The {self.missed_codons} codons were not captured by the employed codons table.')
            
        return self.proteins
    
    def blast_protein(self,  # https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome
                      sequence: str = None,
                      database: str = 'nr', # the blastp database that will be searched with the collected FASTA sequences
                      description: str = '',  # a description of the sequence
                      ):
        if sequence:
            self.multi_fasta = self.make_fasta(sequence) 
            
        # estimate the completion time
        estimated_time = datetime.datetime.now()+datetime.timedelta(seconds = len(self.multi_fasta)*3)    # approximately one second per amino acid   
        print(f'The database search for the parameterized protein(s) will complete circa {estimated_time}.')
        
        # acquire the BLAST results
        self.protein_blast_results = NCBIWWW.qblast('blastp', database, self.multi_fasta)
        
    def blast_nucleotide(self,
                         sequence: str = None,
                         database: str = 'nt',
                         ):
        self.nucleotide_blast_results = []
        if sequence:
            self.sequence = sequence 
            
        # estimate the completion time
        estimated_time = datetime.datetime.now()+datetime.timedelta(seconds = len(self.sequence)/10)    # approximately one second per nucleic acid
        print(f'The database search for the parameterized protein(s) will complete circa {estimated_time}.')
        
        # acquire the BLAST results
        sections = ceil(len(self.sequence)/5000)
        section_size = int(len(self.sequence)/sections)
        sequence_sections = [sequence[i:i+section_size] for i in range(0, sections, section_size)]
        for sequence in sequence_sections:
            nucleotide_blast_result = NCBIWWW.qblast('blastn', database, self.sequence)
            self.nucleotide_blast_results.append(nucleotide_blast_result)
            
        self.nucleotide_blast_results = '\n\n'.join(self.nucleotide_blast_results)
                
    def export(self, export_name = None, export_directory = None):
        # define the simulation_path
        if export_directory is None:
            export_directory = os.getcwd()
        elif not os.path.exists(export_directory):
            error = '--> ERROR: The provided directory does not exist'
            print(error)
            self.messages.append(error)

        translation = ''
        transcription = ''
        ending = ''
        if self.proteins != []:
            translation = 'translation'
            ending = f'-{len(self.proteins)}_proteins'
        if self.transcribed_sequence:
            transcription = 'transcription'
        if export_name is None:
            export_name = '-'.join([re.sub(' ', '_', str(x)) for x in [datetime.date.today(), 'codons', translation, transcription]])
        export_name += ending
            
        count = 0
        export_path = os.path.join(export_directory, export_name)
        while os.path.exists(export_path):
            if not re.search('(-[0-9]+$)', export_path):
                export_path += f'-{count}'   
            else:
                export_path = re.sub('([0-9]+)$', str(count), export_name)
            count += 1
            
        os.mkdir(export_path)
        
        # export the genetic and protein sequences
        self.paths['genetic_sequence'] = os.path.join(export_path, 'genetic_sequence.fasta')
        with open(self.paths['genetic_sequence'], 'w') as genes:
            genes.write(self.make_fasta(self.sequence))
            
        if self.proteins != []:
            self.paths['protein_sequence'] = os.path.join(export_path, 'protein_sequence.fasta')
            with open(self.paths['protein_sequence'], 'w') as proteins:
                proteins.write(self.multi_fasta)
            
        if self.transcribed_sequence:
            self.paths['transcribed_sequence'] = os.path.join(export_path, 'transcribed_sequence.fasta')
            with open(self.paths['transcribed_sequence'], 'w') as genes:
                genes.write(self.transcribed_sequence)
                
        if self.protein_blast_results:
            self.paths['protein_blast_results'] = os.path.join(export_path, 'protein_blast_results.xml')
            with open(self.paths['protein_blast_results'], 'w') as protein_data:
                protein_data.write(self.protein_blast_results.read())
                
        if self.nucleotide_blast_results:
            self.paths['nucleotide_blast_results'] = os.path.join(export_path, 'nucleotide_blast_results.xml')
            with open(self.paths['nucleotide_blast_results'], 'w') as protein_data:
                protein_data.write(self.nucleotide_blast_results.read())
                
                 