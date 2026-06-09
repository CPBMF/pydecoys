pydecoys.strategies
===================

.. currentmodule:: pydecoys.strategies

.. automodule:: pydecoys.strategies

**Types:**

.. autotype:: SeqLike

   `SeqLike` objects can be indexed and spliced; `str` at runtime.

   .. note:: This type doesn't depend on `Biopython`.


.. autotype:: DecoyGenerator

   TypeAlias specifying the signature for decoy strategies.

   A decoy strategy should be a ``Callable[[T], T]`` where ``T`` is a
   :type:`SeqLike`.

**Global CONSTS:**

.. data:: RAND
   :type: typing.Final[random.Random]
   :value: random.Random(_SEED)

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

.. data:: STD_AMINOACIDS
   :type: typing.Final[str]
   :value: 'QWERTYIPASDFGHKLCVNM'

   Standard 20 aminoacids single-letter codes, majuscule.

.. data:: EXT_AMINOACIDS
   :type: typing.Final[str]
   :value: STD_AMINOACIDS + 'OU' + 'BJZX'

   Extended aminoacids single-letter codes, majuscule.

   Non-standard aminoacids
   -----------------------

   - O: Pyrrolysine
   - U: Selenocysteine

   Special codes
   -------------

   - B: Aspartic acid (D) or Asparagine (N)
   - J: Leucine (L) or Isoleucine (I)
   - Z: Glutamic acid (E) or Glutamine (Q)
   - X: Any aminoacid

**Protocols and ABCs:**

.. autoclass:: ContextfulGenerator
   :show-inheritance:
   :special-members: __call__
   :members:

   .. note::
      Functions that don't pass the targets to :meth:`learn_context` are
      documented as such.

.. autoclass:: EnzymeSpecificGenerator
   :show-inheritance:
   :special-members: __call__
   :members:

**Utils:**

.. autofunction:: seq_cast

**Factory functions:**

.. autofunction:: keepsn
.. autofunction:: keepsc
.. autofunction:: keepsterm

**Decorators:**

.. autofunction:: register_function
.. autofunction:: register_class
.. autofunction:: register_cleavage_aware

**Register functions:**

.. autofunction:: register_cleavage_agent
.. autofunction:: register_callable

**Registry utils:**

.. autofunction:: view_strategies
.. autofunction:: view_cleavage_agents
