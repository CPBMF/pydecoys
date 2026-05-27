Using the PyDecoys API
======================

.. currentmodule:: pydecoys

PyDecoys has a comprehensive API for integration in Python proteomics
workflows.

This section provides some basic examples. Refer to the API Reference for more
detailed explanations.

Basic funcionality
------------------

PyDecoys' main API contains useful functions that take a set of target proteins
and return decoys. Each function uses a different Python `type` as
representation for targets and decoys.

IO functions
^^^^^^^^^^^^

To read a target fasta file and get decoys from it, use :func:`from_fasta`:


>>> from Bio import SeqIO
>>> import pydecoys
>>> decoys = pydecoys.from_fasta('path/to/targets.fasta', strategy='reverse')

With :func:`to_fasta`, you can write decoys to a fasta file using a target set:

>>> targets = SeqIO.parse('path/to/targets.fasta', format='fasta')]
>>> pydecoys.to_fasta(targets, 'path/to/decoys.fasta', strategy='reverse')

For convenience, you can directly pass a path or file handle to
:func:`to_fasta`:

>>> pydecoys.to_fasta(
...     'path/to/targets.fasta',
...     'path/to/decoys.fasta',
...     strategy='reverse'
... )

You might want to concatenate the targets and decoys files into a single fasta
file. For that, simply set the ``concat`` parameter to ``True``:

>>> count = pydecoys.to_fasta(
...     'path/to/targets.fasta',
...     'path/to/concat.fasta',
...     strategy='reverse',
...     # the resulting fasta will contain both targets and decoys
...     concat=True
... )

Iterable functions
^^^^^^^^^^^^^^^^^^

Sometimes you might want to hold your data a bit longer before saving it to a
file. PyDecoys has dedicated functions to read from and get iterable protein
sets. For example, to convert a target ``SeqRecord`` iterable into a decoy
``SeqRecord`` iterable:

.. code-block:: python
    :linenos:

    from Bio import SeqIO
    import pydecoys

    targets = SeqIO.parse('path/to/targets.fasta', format='fasta')
    decoys = pydecoys.from_SeqRecords(targets, strategy='reverse')

Similar functions are also available for ``tuple`` and ``str`` objects
representing a protein (not requiring `Biopython`):

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

PyDecoys also has some functions to convert a single target to a single decoy.
As with the previous example, there's a dedicated function for `SeqRecord`
objects:

.. code-block:: python
    :linenos:

    from types import GeneratorType
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord

    # `Bio.SeqIO.read` returns only one `SeqRecord`
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
    sequence, :func:`from_tuples`, :func:`from_seqs`, :func:`tuple_as_decoy`
    and :func:`seq_as_decoy` can take ``Seq`` or ``MutableSeq`` objects as aa
    sequences, and will return the correct type.

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

Custom strategy functions can also be passed directly:

.. code-block:: python
    :linenos:

    import pydecoys
    from pydecoys.strategies import SeqLike


    def custom_decoys(sequence: SeqLike) -> SeqLike:
        ...

    decoys = pydecoys.from_fasta('path/to/targets.fasta', custom_decoys)
