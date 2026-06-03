Creating custom strategies
==========================

.. currentmodule:: pydecoys

The :py:mod:`pydecoys.strategies` module exposes the API used for decoy
strategies.

The basic type signature all strategies must follow is
:py:type:`strategies.DecoyGenerator`. This is as simple type alias:

.. code-block:: python
    :linenos:

    from typing import Callable

    type SeqLike = 'str | Seq | MutableSeq'
    type DecoyGenerator[T: SeqLike] = Callable[[T], T]

Note that :type:`SeqLike` is a Union type of ``str``, ``Seq`` and
``MutableSeq``, but it doesn't actually require `Biopython`.

Creating a simple strategy
--------------------------

To illustrate how to create custom strategies, let's create a naïve randomizer.
It takes a target protein and returns a fully random decoy protein of same
length:

.. code-block:: python
    :linenos:

    # Use PyDecoys own RNG or another RNG with a fixed seed to ensure
    # reproducibility
    from pydecoys.strategies import RAND

    # PyDecoys also has a str containing all 20 standard aminoacid
    # single-letter codes
    from pydecoys.strategies import AMINOACIDS

    def random(sequence: str) -> str:
        length = len(sequence)
        new = RAND.choices(AMINOACIDS, k=length)
        return "".join(new)

This will work perfectly already, and for most usecases it'll be enough.
However, it might fire static type checkers: ``random`` doesn't yet implement
the correct type signature we saw earlier.

Also, our implementation isn't interfaceable with Biopython, which means it
might break (or at least return a ``str``) when it should return ``Seq`` or
``MutableSeq``.

Fixing these problems is easy:

.. code-block:: python
    :linenos:

    from pydecoys.strategies import RAND, AMINOACIDS, SeqLike, seq_cast

    def random[T: SeqLike](sequence: T) -> T:
        length = len(sequence)
        new = RAND.choices(AMINOACIDS, k=length)
        # `seq_cast` handles casting the return value to `sequence`'s type
        return seq_cast(sequence, "".join(new))

Note that :type:`strategies.SeqLike` and :func:`strategies.seq_cast` work
without `Biopython` installed. Although this function correctly interfaces with
`Biopython`, it doesn't need `Biopython` to work!

Strategies with state
---------------------

In the previous example, we implemented a naïve randomizer. In practice, a much
more useful strategy is to give each aminoacid a weighted likelihood based on
the proportions of aminoacis from the target database. For this, the generator
would need some way of gathering context from the database prior to decoy
generation.

PyDecoys has a specific protocol for that:
:class:`strategies.ContextfulGenerator`. When a function that takes a set of
targets receives a :class:`strategies.ContextfulGenerator` as strategy, it
can pass all of the set database to it beforehand.

We need a :type:`strategies.DecoyGenerator`, so we'll override its
``__call__`` method. We also need to implement the
:meth:`strategies.ContextfulGenerator.learn_context` and
:meth:`strategies.ContextfulGenerator.reset` methods, and add a
:attr:`strategies.ContextfulGenerator.is_set` attribute.

.. code-block:: python
    :linenos:

    from typing import Iterable

    from pydecoys.strategies import (
        RAND,
        AMINOACIDS,
        SeqLike,
        seq_cast,
    )

    # Since `ContextfulGenerator` is a protocol, no need to inherit it
    class SmartRandomizer:
        def __init__(self):
            # We'll use a None placeholder until contextualization
            self._weights = None

        def learn_context(self, sequences: Iterable[SeqLike]) -> None:
            # We init weights
            self._weights = [0] * 20

            for seq in sequences:
                for aa in seq:
                    # The position of each weight must overlap with its aminoacid
                    # pos
                    pos = AMINOACIDS.find(aa)
                    self._weights[pos] += 1

        def reset(self):
            self._weights = None

        @property
        def is_set(self) -> bool:
            # We can use a property to ensure is_set is always updated
            return self._weights is not None

        def __call__[T: SeqLike](self, sequence: T) -> T:
            length = len(sequence)
            new = RAND.choices(AMINOACIDS, weights=self._weights, k=length)
            return seq_cast(sequence, "".join(new))

This gets us a working decoy strategy that correctly implements weights before
running! When calling :mod:`pydecoys` IO or Iterable functions, it'll
automatically pass the database to the instance if it's unset, and reset it
again at the end.

.. warning::
    When using IO or Iterable functions from PyDecoys with an unset context-dependant
    strategy, the API will handle setting it, generating the decoys and then resetting
    it at the end. For single-data functions, using an unset strategy will raise a
    ``ValueError``, since they never receive a full dataset. See
    :func:`pydecoys.get_contextualized_strategy` to circumvent this.

Adding new enzymes
------------------

The :class:`strategies.ReversePep`, :class:`strategies.ShufflePep` and
:class:`strategies.RandomizePep` classes allow you to set new enzyme
specifications for `reversepep`, `shufflepep` and `randomizepep`
strategies. Setting new enzymes is an easy instantiation. Let's set
high-specificity chymotrypsin for `reversepep`:

.. code-block:: python
    :linenos:

    import pydecoys
    from pydecoys.strategies import ReversePep, ShufflePep, RandomizePep

    reversepep_chymohs = ReversePep('FWY', nocut='P', sense='C')
    pydecoys.register('reversepep-chymohs', reversepep_chymohs)

That's it. Now you can use a `reversepep` strategy with high specificity
chymotrypsin by providing ``'reversepep-chymohs'`` or ``reversepep_chymohs``
as a strategy.

Enzyme-specific strategies
--------------------------

You might want to set your own enzyme-specific strategy. Luckily, there's an
ABC for that: :class:`strategies.EnzymeSpecificGenerator`. This class sets
up the :meth:`strategies.EnzymeSpecificGenerator.__init__` method we used
earlier to add the new enzyme.

Most importantly, it sets the
:meth:`strategies.EnzymeSpecificGenerator.split_sequence` method. This method
splits a given sequence into its enzymatic fragments (minus cleavage sites) and
the cleavage sites themselves. It yields tuples of an enzymatic fragments and
`False` or a cleavage site and `True`, in the order they appear in the
sequence.

Let's redo the naïve randomizer, but this time let's randomize peptides:

.. code-block:: python
    :linenos:

    from collections.abc import Sequence

    from pydecoys.strategies import (
        RAND,
        AMINOACIDS,
        SeqLike,
        seq_cast,
        EnzymeSpecificGenerator
    )

    class RandomizePep(EnzymeSpecificGenerator):
        def __call__[T: SeqLike](self, sequence: T) -> T:
            decoy_list = []
            for frag, cleavage in self.split_sequence(sequence):
                if not cleavage:
                    length = len(pep)
                    new = RAND.choices(AMINOACIDS, k=length)
                    frag = "".join(new)
                decoy_list.append(frag)

            decoy = "".join(decoy_list)
            return seq_cast(sequence, decoy)

    # We can now create pep randomizers:
    randompep_trypsin = RandomPep('KR', nocut='P', sense='C')

Keeping N- and C-termini
------------------------

You might also want to keep the N-, C- or both termini from the target protein
in the decoy. To accomplish that, use one of the following factories:
:func:`strategies.keepns`, :func:`strategies.keepsc` and
:func:`strategies.keepsterm`. Each factory returns a new callable
that preserves the terminal aminoacids from the argument callable.

:py:class:`strategies.ContextfulGenerator` objects will preserve their
functionality, but will discard the aminoacid that shouldn't be altered from
each sequence when learning context.

.. code-block:: python
    :linenos:

    from pydecoys.strategies import keepsn

    # The `sequence` value is passed directly without the aminoacid that
    # should be preserved. The aa is reinserted after the wrapped function
    # returns.
    random_keepn = keepsn(random)
