Using the PyDecoys API
======================

.. currentmodule:: pydecoys

PyDecoys has a comprehensive API for integration in Python proteomics
workflows.

This section provides some basic examples. Refer to the API Reference for more
detailed explanations.

Basic funcionality
------------------

At the top level, PyDecoys exposes a bunch of useful functions that take a set
of target proteins and returns decoys. Each function uses a different
representation for targets and decoys.

IO functions
^^^^^^^^^^^^

To read a fasta file:


>>> from Bio import SeqIO
>>> import pydecoys
>>> decoys = pydecoys.from_fasta('path/to/targets.fasta', strategy='reverse')

You can also write to a fasta file:

>>> targets = SeqIO.parse('path/to/targets.fasta', format='fasta')]
>>> pydecoys.to_fasta(targets, 'path/to/decoys.fasta', strategy='reverse')

For convenience, you can directly pass a path or file handle to
:py:func:`to_fasta`:

>>> pydecoys.to_fasta(
...     'path/to/targets.fasta',
...     'path/to/decoys.fasta',
...     strategy='reverse'
... )

Finally, you might want to concatenate the targets and decoys files into a
single fasta file:

>>> count = pydecoys.to_fasta(
...     'path/to/targets.fasta',
...     'path/to/concat.fasta',
...     strategy='reverse',
...     concat=True
... )

Iterator functions
^^^^^^^^^^^^^^^^^^

Sometimes you might want to hold your data a bit longer before saving it to a
file. You can convert iterators of targets to decoys:

.. code-block:: python
    :linenos:

    from Bio import SeqIO
    import pydecoys

    targets = SeqIO.parse('path/to/targets.fasta', format='fasta')
    decoys = pydecoys.from_SeqRecords(targets, strategy='reverse')

The following functions don't require `Biopython`:

.. code-block:: python
    :linenos:

    import pydecoys

    tuple_targets: list[tuple[str, str]] = [
        ('seq1', 'DNIDYKAVYR'),
        ('seq2', 'QSYMCTVTHP'),
        ('seq3', 'CQWSLTEELL'),
    ]
    # Decoys are a `Generator` of `tuple`
    decoys = pydecoys.from_tuples(tuple_targets, strategy='reverse')

    str_targets: list[str] = [
        'DNIDYKAVYR',
        'QSYMCTVTHP',
        'CQWSLTEELL',
    ]
    # Decoys are a `Generator` of `str`
    decoys = pydecoys.from_seqs(str_targets, strategy='reverse')

Single data functions
^^^^^^^^^^^^^^^^^^^^^

PyDecoys has some functions to convert a single target to a single decoy:

.. code-block:: python
    :linenos:

    from types import GeneratorType
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord

    target = SeqIO.read('path/to/target.fasta', format='fasta')
    decoy = pydecoys.SeqRecord_as_decoy(target, strategy='reverse')
    assert instance(decoys, SeqRecord)

    tuple_target = ('seq1', 'DNIDYKAVYR')
    decoy = pydecoys.tuple_as_decoy(target, strategy='reverse')
    assert instance(decoy, tuple)

    str_target = 'DNIDYKAVYR'
    decoy = pydecoys.seq_as_decoy(str_target, strategy='reverse')
    assert instance(decoy, str)

.. note::
    Although the examples above use a ``str`` to represent an aminoacid
    sequence, :py:func:`from_tuples`, :py:func:`from_seqs`,
    :py:func:`tuple_as_decoy` and :py:func:`seq_as_decoy` can take ``Seq`` or
    ``MutableSeq`` objects as aa sequences, and will return the correct type.

Adding new strategies
---------------------

You can register a new decoy strategy using the register function:

.. code-block:: python
    :linenos:

    import pydecoys
    from pydecoys.strategies import SeqLike


    def custom_decoys(sequence: SeqLike) -> SeqLike:
        ...

    pydecoys.register('custom-decoys', custom_decoys)
    decoys = pydecoys.from_fasta('path/to/targets.fasta', 'custom-decoys')

Custom strategies can be used directly as functions:

.. code-block:: python
    :linenos:

    import pydecoys
    from pydecoys.strategies import SeqLike


    def custom_decoys(sequence: SeqLike) -> SeqLike:
        ...

    decoys = pydecoys.from_fasta('path/to/targets.fasta', custom_decoys)
