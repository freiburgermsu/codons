from pandas import DataFrame
from shutil import rmtree
from math import isclose
from glob import glob
import codons
import re, os


def test_init():
    # import the class modules
    cd = codons.Codons()
    for TF in [cd.verbose, cd.final, cd.printing, cd.verbose]:
        assert type(TF) is bool
    for dic in [cd.paths, cd.proteins, cd.parameters, cd.codons_table, cd.]:
        assert type(dic) is dict
    for path in ['changed_codons', 'standard_table', 'amino_acid_synonyms']:
        assert type(cd.paths[path]) is bool
    for string in [cd.sequence, cd.parameters['residue_delimiter']:
        assert type(string) is str
           
def test_transcribe():
    # calculate the MW for chemicals of known MW 
    dna_sequence = 'cactaagaaa gatgctgctg ctgctaaaaa taagatgcgc cacaagcgca cttccaccaa'
    
    # calculate the MW for the dictionary of chemicals    
    chem_mw = chemw.ChemMW()
    for chemical in test_chemicals:
        chem_mw.mass(chemical)
        tolerance = chem_mw.mw*0.001 # 99.9% accuracy
        if not isclose(chem_mw.raw_mw, test_chemicals[chemical], rel_tol = tolerance):
            assert False
        else:
            assert True
            
    # affirm that iterated entities are zero
    for zero in [chem_mw.groups, chem_mw.layer, chem_mw.skip_characters]:
        assert zero == 0
            