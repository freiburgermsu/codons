Execution
---------------

Codons is executed through the following sequence of the aforementioned functions, which is exemplified in the `example Notebook of our GitHub repository <./examples>`_:

.. code-block:: python

   import codons
   cd = codons.Codons(sequence = None, codons_table = 'standard', amino_acids_form = 'full_name', hyphenated = None, verbose = False, printing = True)
   # < Codons function(s) > 
   cd.export(export_name = None, export_directory = None)