pydecoys
========

.. currentmodule:: pydecoys

Provides:

   #. A CLI app to easily generate a decoy fasta file from a target fasta file
   #. A comprehensive Python API with IO, Generator and single data functions
   #. Easy implementation of custom decoy strategies
   #. Full integration with Biopython

.. note::
   Some functionalities of PyDecoys require Biopython. Those functionalities
   are documented as such.

**Functions that perform IO operations:**

.. autofunction:: from_fasta
.. autofunction:: to_fasta

**Functions that return an iterator of decoy sequences from an interator of
target sequences:**

.. autofunction:: from_SeqRecords
.. autofunction:: from_tuples
.. autofunction:: from_seqs

.. note::
   Those functions try to be lazy whenver possible, but since
   :class:`strategies.ContextulGenerator` needs prior access to the target
   database, lazily loading targets isn't possible with context-dependant
   strategies. Thus, target sequences are collected into a list beforehand when
   using a context-dependant strategy.

   Since IO functions all use :func:`from_SeqRecords` internally, this is also
   valid for them.

**Functions that return a single decoy sequence from a single target
sequence:**

.. autofunction:: SeqRecord_as_decoy
.. autofunction:: tuple_as_decoy
.. autofunction:: seq_as_decoy

**Utils:**

.. autofunction:: register
.. autofunction:: get_contextualized_strategy
