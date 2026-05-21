PyDecoys
========

.. automodule:: pydecoys

.. note::
   Some functionalities of PyDecoys require Biopython. Those functionalities
   are documented as such.

IO functions
------------
Functions that perform IO operations.

.. autofunction:: from_fasta
.. autofunction:: to_fasta

Iterator functions
------------------
Functions that return an iterator of decoy sequences from an interator of
target sequences.

.. autofunction:: from_SeqRecords
.. autofunction:: from_seqs
.. autofunction:: from_tuples

.. note::
   Those functions try to be lazy whenver possible, but since
   :class:`strategies.ContextulGenerator` needs prior access to the target
   database, lazily loading targets isn't possible with context-dependant
   strategies. Thus, target sequences are collected into a list beforehand when
   using a context-dependant strategy.

   Since IO functions all use :func:`from_SeqRecords` internally, this is also
   valid for them.

Single record functions
-----------------------
Functions that return a single decoy sequence from a single target sequence.

.. autofunction:: SeqRecord_as_decoy
.. autofunction:: seq_as_decoy
.. autofunction:: tuple_as_decoy

.. note::
   Those functions do not give context to
   :class:`strategies.ContextfulGenerator` objects. Context should be loaded
   beforehad.

Utils
-----
.. autofunction:: register
