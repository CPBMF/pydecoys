strategies
==========

.. currentmodule:: pydecoys.strategies

.. automodule:: pydecoys.strategies

.. autotype:: SeqLike
.. autotype:: Seq_
.. autotype:: MutableSeq_

.. autodata:: RAND

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

.. autodata:: AMINOACIDS

Neither of those types need `Biopython` installed.

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
