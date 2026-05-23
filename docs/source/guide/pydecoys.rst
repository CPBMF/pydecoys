Using the PyDecoys API
======================
PyDecoys has a comprehensive API for integration in Python proteomics
workflows.

This section provides some basic examples. Refer to the API Reference for more
detailed explanations.

Basic funcionality
------------------
At the top level, PyDecoys exposes a bunch of useful functions that take a set
of target proteins and returns decoys. Each function uses a different
representation for targets and decoys.

To read and write to fasta files:

.. code-block:: python

    from Bio import SeqIO
    import pydecoys

    # Get decoys from a set of targets in a fasta file
    decoys = pydecoys.from_fasta('path/to/targets.fasta', strategy='reverse')

    # Write decoys to a file from a set of targets
    targets = SeqIO.parse('path/to/targets.fasta', format='fasta')]
    pydecoys.to_fasta(targets, 'path/to/decoys.fasta', strategy='reverse')

    # For convenience, you can directly pass a path or file handle to `to_fasta`
    pydecoys.to_fasta('path/to/targets.fasta', 'path/to/decoys.fasta', strategy='reverse')


You can also convert iterators of targets to decoys:

.. code-block:: python
    :linenos:

    from Bio import SeqIO
    import pydecoys

    targets = SeqIO.parse('path/to/targets.fasta', format='fasta')
    decoys = pydecoys.from_SeqRecords(targets, strategy='reverse')

    # The following functions don't require Biopython:

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

The following functions don't need to take an iterator:

.. code-block:: python

    from types import GeneratorType
    from Bio import SeqIO

    target = SeqIO.read('path/to/target.fasta', format='fasta')
    decoys = pydecoys.from_SeqRecords(target, strategy='reverse')
    assert instance(decoys, GeneratorType)  # Results are still iterable

    str_target = 'DNIDYKAVYR'
    decoys = pydecoys.from_seqs(str_target, strategy='reverse')
    assert instance(decoys, GeneratorType)  # Results are still iterable

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
    # Or pass the function directly as strategy
    decoys = pydecoys.from_fasta('path/to/targets.fasta', custom_decoys)
