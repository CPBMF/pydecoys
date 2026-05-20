strategies
==========

.. currentmodule:: pydecoys.strategies

.. automodule:: pydecoys.strategies

Interfaces
----------
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
      already given context or not. To avoid passing context

Reverse funtions
----------------
.. autofunction:: reverse
.. autofunction:: reverse_keep_n
.. autofunction:: reverse_keep_c
.. autofunction:: reverse_keep_term

Shuffle functions
-----------------
.. autofunction:: shuffle
.. autofunction:: shuffle_keep_n
.. autofunction:: shuffle_keep_c
.. autofunction:: shuffle_keep_term

Pseudo-reverse and -shuffle
---------------------------
Pre-initialized pseudo-reversers and pseudo-shufflers covering most proteases
are available:

.. autodata:: pseudoreverse_trypsin
.. autodata:: pseudoreverse_stricttrypsin
.. autodata:: pseudoreverse_argc
.. autodata:: pseudoreverse_aspn
.. autodata:: pseudoreverse_chymo
.. autodata:: pseudoreverse_gluc
.. autodata:: pseudoreverse_lysc
.. autodata:: pseudoreverse_lysn
.. autodata:: pseudoreverse_stricttrypsin_keepn
.. autodata:: pseudoshuffle_trypsin
.. autodata:: pseudoshuffle_stricttrypsin
.. autodata:: pseudoshuffle_argc
.. autodata:: pseudoshuffle_aspn
.. autodata:: pseudoshuffle_chymo
.. autodata:: pseudoshuffle_gluc
.. autodata:: pseudoshuffle_lysc
.. autodata:: pseudoshuffle_lysn
.. autodata:: pseudoshuffle_stricttrypsin_keepn

.. autoclass:: PseudoReverseRule
   :members:
   :special-members: __call__

.. autoclass:: PseudoShuffleRule
   :members:
   :special-members: __call__

.. note::
   Both classes don't have setters and should be treated as immutable.
   Internally, they use a regex pattern to split given sequences into their
   enzymatic peptides and cleavage sites. Since the regex pattern is constructed
   and compiled at initialization, changing either property would just cause a
   mismatch between the visible specifications and the actual behavior.

Utils
-----
.. autodata:: rand

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

.. autotype:: SeqLike
.. autotype:: Seq_
.. autotype:: MutableSeq_
