Translating, Transcribing, and Investigating Genetic Sequences and their corresponding Proteins 
--------------------------------------------------------------------------------------------------------

|PyPI version| |Actions Status| |Downloads| |License|

.. |PyPI version| image:: https://img.shields.io/pypi/v/codons.svg?logo=PyPI&logoColor=brightgreen
   :target: https://pypi.org/project/codons/
   :alt: PyPI version

.. |Actions Status| image:: https://github.com/freiburgermsu/codons/workflows/Test%20Codons/badge.svg
   :target: https://github.com/freiburgermsu/codons/actions
   :alt: Actions Status

.. |License| image:: https://img.shields.io/badge/License-MIT-blue.svg
   :target: https://opensource.org/licenses/MIT
   :alt: License

.. |Downloads| image:: https://pepy.tech/badge/Codons
   :target: https://pepy.tech/project/Codons
   :alt: Downloads

The Codons module is a lightweight tool for 1) conducting transcription and translation of genetic sequences, either in a FASTA format or a string; 2) making and reading FASTA or multi-FASTA files; and 3) conducting BLAST searches of protein and nucleotide sequences. Example Notebooks of these features are offered in the "examples" directory of the Codons `GitHub repository <https://github.com/freiburgermsu/codons/tree/main/examples>`_.



Installation
+++++++++++++

The following command installs ``Codons`` in a command prompt/terminal environment::
 
 pip install codons

_________________


+++++++++++
__init__
+++++++++++

The data environment, in a `Python IDE <https://www.simplilearn.com/tutorials/python-tutorial/python-ide>`_, is defined: 

.. code-block:: python

  import codons
  cd = codons.Codons(sequence = None, codons_table = 'standard', amino_acids_form = 'full_name', hyphenated = None, verbose = False, printing = True)

- *sequence* ``str``: specifies the genetic sequence that will be processed through subsequent functions, which can alternatively be provided in each function ad hoc.
- *codons_table* ``str``: specifies the framework for translating codons into amino acids, where the `standard translation table <https://en.wikipedia.org/wiki/DNA_and_RNA_codon_tables>`_ is used by default.
- *amino_acids_form* ``str``: specifies whether the amino acid ``full_name``, ``three_letter``, or ``one_letter`` nomenclature will be used in the protein sequence. 
- *hyphenated* ``bool``: specifies whether amino acid residues of the protein sequence are delimited by hyphens, where ``None`` defaults to ``True`` for ``amino_acids_for = full_name`` and ``amino_acids_for = three_letter`` and ``False`` for ``amino_acids_for = one_letter``.
- *verbose* & *printing* ``bool``: specifies whether troubleshooting information or MW results will be printed, respectively.

++++++++++++++++
read_fasta()
++++++++++++++++

A genetic sequence is converted from DNA -> RNA, or RNA -> DNA, where the directionality of the conversion is automatically listed in the FASTA description:

.. code-block:: python

 sequences, descriptions, fasta_file = transcribed_sequence = cd.read_fasta(fasta_path = None, fasta_link = None):

- *fasta_path* ``str``: The path to a FASTA file that will be loaded, parsed, and returned.
- *fasta_link* ``str``: The URL link to a FASTA file that will be imported, parsed, and returned.

**Returns**: 

- *sequences* & *descriptions* ``list``: The sequences and descriptions that are contained within the FASTA file.
- *fasta_file* ``str``: The FASTA file as a string that is specified by the path or URL link argument.

++++++++++++++++
make_fasta()
++++++++++++++++

A simple function that returns, and optionally exports, a FASTA-formatted file from the parameterized description and sequence:

.. code-block:: python

 fasta_file = cd.make_fasta(sequence, description = 'sequence', export_path = None):

- *sequence* ``str``: The genetic or protein sequence that will constitute the FASTA file. 
- *description* ``str``: A description of the sequence that will be the first line of the FASTA file. 
- *export_path* ``str``: The path to which the FASTA file will be exported, where ``None`` specifies that the file will not be exported.
 
**Returns**: 

- *fasta_file* ``str``: The FASTA-formatted file as a string, based upon the parameterized sequence and description.

++++++++++++++++
transcribe()
++++++++++++++++

A genetic sequence is converted from DNA -> RNA, or RNA -> DNA, where the directionality of the conversion is automatically listed in the FASTA description:

.. code-block:: python

 transcribed_sequence = cd.transcribe(sequence = None, description = '', fasta_path = None, fasta_link = None)

- *sequence* ``str``: The genetic seqeuence that will be transcribed. The sequence is case-insensitive, and can even possess line numbers or column-spaces, which the code ignores. The absence of a passed sequence executes the sequence that is loaded into the ``Codons`` object.
- *description* ``str``: A description of the genetic sequence that will be added to the FASTA-formatted output of the function. 
- *fasta_path* & *fasta_link* ``str``: The path or URL link to a FASTA file that will be transcribed.

**Returns**: 

- *transcribed_sequence* ``str``: The translated sequence as a single string.

++++++++++++++++
translate()
++++++++++++++++

A genetic sequence is translated into a FASTA-formatted sequence of amino acids for each protein that is coded by the genetic code:

.. code-block:: python

 proteins = cd.translate(sequence = None, fasta_path = None, fasta_link = None)

- *sequence* ``str``: The genetic sequence , of either DNA or RNA, that will be translated into a protein sequence. The sequence is case-insensitive, and can even possess line numbers or column-spaces, which the code ignores. The absence of a passed sequence executes the sequence that is loaded into the ``Codons`` object.
- *fasta_path* & *fasta_link* ``str``: The path or URL link to a FASTA file that will be translated.


++++++++++++++++
blast_protein()
++++++++++++++++

A protein sequence or a FASTA-formatted file of protein sequences is searched in through the `BLAST database <https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp&PAGE_TYPE=BlastSearch&BLAST_SPEC=&LINK_LOC=blasttab&LAST_PAGE=blastn>`_ of the NIH for information about the protein(s):

.. code-block:: python

 blast_results = cd.blast_protein(sequence = None, database = 'nr', description = 'Protein sequence description',  fasta_path = None, fasta_link = None,  export_name = None, export_directory = None)

- *sequence* ``str``: The genetic seqeuence, of either DNA or RNA, that will be parsed and translated into a protein sequence. The sequence is case-insensitive, and can even possess line numbers or column-spaces, which the code ignores. The absence of a passed sequence executes the sequence that is loaded into the ``Codons`` object.
- *database* ``str``: The BLAST database that will be searched for the protein sequence. Permissible options include: ``nr``, ``refseq_select``, ``refseq_protein``, ``landmark``, ``swissprot``, ``pataa``, ``pdb``, ``env_nr``, ``tsa_nr``.
- *description* ``str``: A description of the genetic sequence that will be added to the FASTA-formatted output of the function. 
- *fasta_path* & *fasta_link* ``str``: The path or URL link to a protein FASTA or multi-FASTA file that will be systematically searched.
- *export_name* & *export_directory* ``str``: The name of the folder and directory to which the scraped BLAST data will be saved in a file: ``nucleotide_blast_results.xml``. The ``None`` values enable the code to construct a unique folder name that describes the contents and saves it to the current working directory.

**Returns**

- *blast_results* `Bio.Blast.NCBIXML`: An API accessible format that facilitates investigating the acquired BLAST from the search content.


++++++++++++++++++++++++++++
blast_nucleotide()
++++++++++++++++++++++++++++

A genetic sequence is translated into a FASTA-formatted sequence of amino acids for each protein that is coded by the genetic code:

.. code-block:: python

 cd.blast_nucleotide(sequence = None, database= 'nt', description = 'Genetic sequence description', export_name = None, export_directory = None)

- *sequence* ``str``: The genetic sequence, of either DNA or RNA, that will be parsed and translated into a protein sequence. The sequence is case-insensitive, and can even possess line numbers or column-spaces, which the code ignores. The absence of a passed sequence executes the sequence that is loaded into the ``Codons`` object.
- *database* ``str``: The BLAST database that will be searched for the nucleotide sequence. Permissible options include: ``nr``, ``nt``, ``refseq_select``, ``refseq_rna``, ``refseq_representative_genomes``, ``wgs``, ``refseq_genomes``, ``est``, ``SRA``, ``TSA``, ``HTGS``, ``pat``, ``pdb``, ``RefSeq_Gene``, ``gss``, ``dbsts``.
- *description* ``str``: A description of the genetic sequence that will be added to the FASTA-formatted output of the function. 
- *fasta_path* & *fasta_link* ``str``: The path or URL link to a protein FASTA or multi-FASTA file that will be systematically searched.
- *export_name* & *export_directory* ``str``: The name of the folder and directory to which the scraped BLAST data will be saved in a file: ``protein_blast_results.xml``. The ``None`` values enable the code to construct a unique folder name that describes the contents and saves it to the current working directory.

**Returns**

- *blast_results* `Bio.Blast.NCBIXML`: An API accessible format that facilitates investigating the acquired BLAST from the search content.

++++++++++++++++
export()
++++++++++++++++

The genetic sequence and any corresponding protein or nucleotide content from the aforementioned functions, which reside in the ``Codons`` object, are exported:

.. code-block:: python

 cd.export(export_name = None, export_directory = None)

- *export_name* ``str``: optionally specifies a name for the folder of exported content, where `None` enables the code to design a unique folder name for simulation and descriptive tags of the contents.
- *export_directory* ``str``: optionally specifies a path to where the content will be exported, where `None` selects the current working directory.

++++++++++++++++++++++++++
Accessible content
++++++++++++++++++++++++++
The ``Codons`` object retains numerous components that are accessible to the user: 

- *genes* ``dict``: A dictionary of all genes in the genetic sequence, with sub-content of a list of all coding Codons in the gene and the corresponding protein sequence and mass.
- *protein_fasta* & *gene_fasta* ``str``: Assembled FASTA-formatted files for the translated proteins from a parameterized genetic sequence and for a genetic sequence, respectively.
- *transcribed_sequence* & *sequence* ``str``: The transcribed genetic sequence from the ``transcription()`` function, and the genetical sequence that is used in any of the ``Codons`` functions.
- *amino_acid_synonyms* ``dict``: The synonyms for each amino acid, with keys of the full amino acid name.
- *codons_table* & *changed_codons* ``dict``: The translation table between genetic codons and amino acid residues, which is accessed with case-insensitivity, and the translation conversions that were changed based upon the user's specification.
- *missed_codons* ``dict``: A collections of the codons that were parsed yet never matched with a ``codons_table`` key.
- *paths* & *parameters* ``dict``: Collections of the paths and parameters are are defined for the simulation.
- *export_path* ``str``: The complete export path for the ``Codons`` contents.
- *protein_blast_results* & *nucleotide_blast_results* ``str``: The BLAST search results for the passed proteins and nucleotides, respectively.

Execution
+++++++++++

Codons is executed through the following sequence of the aforementioned functions, which is exemplified in the `example Notebook of our GitHub repository <./examples>`_:

.. code-block:: python

   import codons
   cd = codons.Codons(sequence = None, codons_table = 'standard', amino_acids_form = 'full_name', hyphenated = None, verbose = False, printing = True)
   # < Codons function(s) > 
   cd.export(export_name = None, export_directory = None)