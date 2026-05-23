Creating custom strategies
==========================

.. currentmodule:: pydecoys

The :py:mod:`pydecoys.strategies` module exposes the API used for decoy strategies.

The basic protocol all strategies must follow is
:py:class:`strategies.DecoyGenerator`. This protocol is just a callable class
with some overloads:

.. code-block:: python
    :linenos:

    class DecoyGenerator(Protocol):
        @overload
        def __call__(self, sequence: Seq) -> Seq: ...

        @overload
        def __call__(self, sequence:  MutableSeq) -> MutableSeq: ...

        @overload
        def __call__(self, sequence: str) -> str: ...

        def __call__(self, sequence: SeqLike) -> SeqLike:
            pass

Creating a simple strategy
--------------------------

To illustrate how to create custom strategies, let's create a naïve randomizer
that takes a target protein and returns a fully random decoy protein of same
size:

.. code-block:: python
    :linenos:

    # Use PyDecoys own RNG or another RNG with a fixed seed to ensure
    # reproducibility
    from pydecoys.strategies import RAND

    # PyDecoys also has a str containing all 20 aminoacid letter-codes
    from pydecoys.strategies import AMINOACIDS

    def random(sequence: str) -> str:
        length = len(sequence)
        new = RAND.choices(AMINOACIDS, k=length)
        return "".join(new)

This will work perfectly already, and for most usecases it'll be enough.
However, it might fire static type checkers: ``random`` doesn't yet implement the
correct overloads. PyDecoys :py:class:`strategies.DecoyGenerator` protocol has
three overloads for ``__call__``.

Our implementation also isn't interfaceable with Biopython, which means it
might break or at least return a ``str`` when it should return ``Seq`` or
``MutableSeq``.

Fixing these problems is easy:

.. code-block:: python
    :linenos:

    from typing import overload

    from Bio.Seq import Seq, MutableSeq
    from pydecoys.strategies import RAND, AMINOACIDS, SeqLike, seq_cast

    @overload
    def random(sequence: Seq) -> Seq: ...

    @overload
    def random(sequence: MutableSeq) -> MutableSeq: ...

    @overload
    def random (sequence: str) -> str: ...

    def random(sequence: SeqLike) -> SeqLike:
        length = len(sequence)
        new = RAND.choices(AMINOACIDS, k=length)
        # `seq_cast` handles casting the return value to `sequence`'s type
        return seq_cast(sequence, "".join(new))

But now we introduced `Biopython` as a dependency!

While :py:type:`strategies.SeqLike` and :py:type:`strategies.seq_cast` work
without `Biopython` installed, we're still using ``Seq`` and ``MutableSeq``. To
solve that, PyDecoys introduces two alias that are ``str`` at runtime, and
don't require `Biopython`:

.. code-block:: python
    :linenos:

    from typing import overload

    from pydecoys.strategies import (
        RAND,
        AMINOACIDS,
        SeqLike,
        seq_cast,
        Seq_,
        MutableSeq_
    )

    @overload
    def random(sequence: Seq_) -> Seq_: ...

    @overload
    def random(sequence: MutableSeq_) -> MutableSeq_: ...

    @overload
    def random (sequence: str) -> str: ...

    def random(sequence: SeqLike) -> SeqLike:
        length = len(sequence)
        new = RAND.choices(AMINOACIDS, k=length)
        return seq_cast(sequence, "".join(new))

Now the function is perfectly typed without ever needing `Biopython` as a
dependency! Better yet, if it does receive `Biopython` objects, it'll handle
them correctly.

Strategies with state
---------------------

In the previous example, we implemented a naïve randomizer. In practice, a much
more useful strategy is to give each aminoacid a weighted likelihood based on
the proportions of aminoacis from the target database. For this, the generator
would need some way of gathering context from the database prior to decoy
generation, so it could save it as state.

PyDecoys has a specific protocol for that:
:py:class:`strategies.ContextfulGenerator`. When a function that takes a set of
targets receives a :py:class:`strategies.ContextfulGenerator` as strategy, it
passes all the set database to it beforehand.

It's a :py:class:`strategies.DecoyGenerator`, so we'll need to override its
``__call__`` method. We also need to override the
:py:meth:`strategies.ContextfulGenerator.learn_context` method. This is
the method our implementation will use to gather its state.

.. code-block:: python
    :linenos:

    # Boilerplate imports
    from collections.abc import Sequence
    from typing import overload

    from pydecoys.strategies import (
        RAND,
        AMINOACIDS,
        SeqLike,
        seq_cast,
        Seq_,
        MutableSeq_
    )

    # Since `ContextfulGenerator` is a protocol, no need to inherit it
    class SmartRandomizer:
        def __init__(self):
            # We just define a list of 0's
            self._weights = [0] * 20

        def learn_context(self, sequences: Sequence[SeqLike]) -> None:
            # We'll use the whole dataset to populate weights
            for seq in sequences:
                for aa in seq:
                    # The position of each weight must overlap with its aminoacid
                    # pos
                    pos = AMINOACIDS.find(aa)
                    self._weights[pos] += 1

        # Boilerplate overloads
        @overload
        def __call__(self, sequence: Seq_) -> Seq_:

        @overload
        def __call__(self, sequence: MutableSeq_) -> MutableSeq_:

        @overload
        def __call__(self, sequence: str) -> str:

        def __call__(self, sequence: SeqLike) -> SeqLike:
            length = len(sequence)
            new = RAND.choices(AMINOACIDS, weights=self._weights, k=length)
            return seq_cast(sequence, "".join(new))

This gets us a working decoy strategy that correctly implements weights before
running! When calling :py:mod:`pydecoys` IO or Iterable functions, it'll
automatically pass the database to the instance. One problem only: if we use
the same instance again with one of those functions, it'll learn context again.

To deal with that, we can use a simple boolean tag:

.. code-block:: python
    :linenos:

    # imports...

    class SmartRandomizer:
        def __init__(self):
            self._weights = [0] * 20
            # Add a flag
            self._context = False

        def learn_context(self, sequences: Sequence[SeqLike]) -> None:
            # Just a simple guard
            if self._context:
                return

            for seq in sequences:
                for aa in seq:
                    pos = AMINOACIDS.find(aa)
                    self._weights[pos] += 1

            # We need to set it afterwards
            self._context = True

        # __call__...

All set! Now it's ready to be used or registered through :py:func:`register`!

New enzymes for ReversePep and ShufflePep
-----------------------------------------
The :py:class:`strategies.ReversePep` and :py:class:`strategies.ShufflePep`
allow you to set new enzyme specifications for `reversepep` and `shufflepep`
strategies. Setting new enzymes is an easy instantiation. Let's set
high-specificity chymotrypsin for each:

.. code-block:: python
    :linenos:

    import pydecoys
    from pydecoys.strategies import ReversePep, ShufflePep

    reversepep_chymohs = ReversePep('FWY', nocut='P', sense='C')
    shufflepep_chymohs = ShufflePep('FWY', nocut='P', sense='C')

That's it. The ``reversepep_chymohs`` and ``shufflepep_chymohs`` are decoy
generators that reverse or shuffle high-specificity chymotrypsin fragments!

If you wish to keep some terminal aminoacids as well, that's simple:

.. code-block:: python
    :linenos:

    # imports...

    reversepep_chymohs_keepn = \
        ReversePep('FWY', nocut='P', sense='C', keep_term='N')

    reversepep_chymohs_keepc = \
        ReversePep('FWY', nocut='P', sense='C', keep_term='C')

    reversepep_chymohs_keepterm = \
        ReversePep('FWY', nocut='P', sense='C', keep_term='both')

Note that you still need to register those instances to use them via a ``str``
key.

Enzyme-specific strategies
--------------------------
You might want to set your own enzyme-specific strategy. Luckily, there's an
ABC for that: :py:class:`strategies.EnzymeSpecificGenerator`. This class sets
up the :py:meth:`strategies.EnzymeSpecificGenerator.__init__` method we used
earlier to add the new enzyme. It also sets
:py:attr:`strategies.EnzymeSpecificGenerator.cut`, 
:py:attr:`strategies.EnzymeSpecificGenerator.nocut`,
:py:attr:`strategies.EnzymeSpecificGenerator.sense` and
:py:attr:`strategies.EnzymeSpecificGenerator.keep_term` get-only properties!

Most importantly, it sets a ``_pattern`` attribute that matches the cleavage
sites and maybe the terminal aminoacids (based on the specifications at
instantiation). This pattern is already compiled and wrapped into a capture
group.

Let's redo the naïve randomizer, but this time let's randomize peptides:

.. code-block:: python
    :linenos:

    # Boilerplate imports
    from collections.abc import Sequence
    from typing import overload

    from pydecoys.strategies import (
        RAND,
        AMINOACIDS,
        SeqLike,
        seq_cast,
        Seq_,
        MutableSeq_,
        EnzymeSpecificGenerator
    )

    class RandomizePep(EnzymeSpecificGenerator):
        @overload
        def __call__(self, sequence: Seq_) -> Seq_:

        @overload
        def __call__(self, sequence: MutableSeq_) -> MutableSeq_:

        @overload
        def __call__(self, sequence: str) -> str:

        def __call__(self, sequence: SeqLike) -> SeqLike:
            # re module requires `str`
            sequence = str(sequence)

            # List of tuples with position and aa of sites to keep
            keep = []
            for pep in re.finditer(self._pattern, sequence):
                keep.append((pep.start(), pep.group()))

            # We save each portion of the decoy to a list
            decoy_list = []

            split = re.split(self._pattern, sequence)
            for pep in split:
                # split can sometimes return an empty str at the start of the
                # list
                if not pep:
                    continue
                length = len(pep)
                new = RAND.choices(AMINOACIDS, k=length)
                decoy_list.append("".join(new))

            decoy = "".join(decoy_list)

            # Reinsert aminoacids that should be kept
            for i, pep in keep:
                decoy = decoy[:i] + keep + decoy[i+1:]

            return seq_cast("".join(decoy))

    # We can now create pep randomizers:
    randompep_trypsin = RandomPep('KR', nocut='P', sense='C')
