pydecoys.strategies
===================

.. currentmodule:: pydecoys.strategies

.. automodule:: pydecoys.strategies

.. type:: SeqLike

   ``(str | Seq | MutableSeq)``

   `SeqLike` objects can be indexed and spliced; `str` at runtime.

.. type:: Seq_

   ``Seq``

   `Seq` type that doesn't require Biopython; `str` at runtime.

.. type:: MutableSeq_

   ``MutableSeq``

   `MutableSeq` type that doesn't require Biopython; `str` at runtime.

.. note:: None of those types need `Biopython` installed.

.. data:: RAND
   :type: typing.Final[random.Random]

   Random number generator for stochastic decoy strategies.

   This RNG has a fixed seed to guarantee reproducibility and repeatability
   of decoy databases and thus any experiments using them. Custom stochastic
   decoy strategies should use this RNG or other RNG with a fixed seed.

   .. warning::
      This doesn't mean that a stochastic decoy generator is deterministic.
      Since using the RNG changes the RNG state itself, the same target protein
      will generate different decoys accross different calls within the same
      Python runtime. This also means that if a protein sequence appears more
      than once in a target dataset, a stochastic method will inflate the decoy
      dataset size.

.. data:: AMINOACIDS
   :type: typing.Final[str]
   :value: 'QWERTYIPASDFGHKLCVNM'

   Standard 20 aminoacids single-letter codes, majuscule.

.. autoclass:: DecoyGenerator
   :show-inheritance:
   :special-members: __call__
   :members:

.. autoclass:: ContextfulGenerator
   :show-inheritance:
   :members:

   .. note::
      Functions that don't pass the targets to a ContextfulGenerator are
      documented as such.

   .. note::
      PyDecoys doesn't discriminate between ContextfulGenerators that have been
      already given context or not. You should implement this yourself if you
      want this behavior.

.. autoclass:: EnzymeSpecificGenerator
   :show-inheritance:
   :special-members: __call__
   :members:

.. autoclass:: ReversePep
   :show-inheritance:
   :special-members: __call__
   :members:

.. autoclass:: ShufflePep
   :show-inheritance:
   :special-members: __call__
   :members:

.. note::
   `EnzymeSpecificGenerator` classes don't have setters and should be treated as
   immutable.

   Internally, they use a regex pattern to identify aminoacids that shouldn't be
   altered (cleavage sites and maybe terminal aminoacids). Since the regex pattern
   is constructed and compiled at initialization, changing either property would
   just cause a mismatch between the visible specifications and the actual
   behavior.
