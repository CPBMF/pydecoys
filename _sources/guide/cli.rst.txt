Running the CLI
===============

Running `pydecoys` is simple. This section provides an overview of the CLI's
functionalities.

.. note::
    See the full list of available strategies at :ref:`available-strategies`.
    The full list of CLI options is at :ref:`cli`.

Decoy fasta from a target fasta
-------------------------------

The first argument `pydecoys` takes is an input fasta of target proteins. You
can also specify the output filename of the decoy proteins fasta with the
``-o/--output`` flag.

.. code-block:: sh

    # To generate a reversed decoy set from a given `targets.fasta`
    pydecoys targets.fasta -o rev.fasta

Different decoy strategies
--------------------------

By default, `pydecoys` will use the ``'reverse'`` decoy strategy. To choose
another strategy, use the ``-s/--strategy`` flag followed by the strategy key.

.. code-block:: sh

    # Change the strategy with the `-s/--strategy` flag
    pydecoys targets.fasta -o shuf.fasta -s shuffle

    # Shuffle tryptic enzymatic fragments, except cleavage sites
    pydecoys targets.fasta -o shuf.fasta -s shufflepep-trypsin

    # Except N-terminal aas and cleavage sites
    pydecoys targets.fasta -o shuf.fasta -s shufflepep-trypsin-keepn

    # Randomize (except N-terminal aas) then prepend each target to its decoy
    pydecoys targets.fasta -o shuf.fasta -s randomize-keepn-fuse

`stdin` and `stdout`
--------------------

Both ``input`` and ``output`` are optional, and by default `pydecoys` will
read from ``stdin`` and write to ``stdout``. As a result, you can pipe from
and to `pydecoys`:

.. code-block:: sh

    # The input defaults to stdin
    cat fasta_1.fasta fasta_2.fasta | pydecoys -o rev_1_2.fasta

    # And the output, to stdout
    pydecoys targets.fasta | grep some_seq

    # So you can chain commands
    cat fasta_1.fasta fasta_2.fasta | pydecoys | grep some_seq

Concatenating targets and decoys
--------------------------------

Sometimes, you might want one fasta file containing the targets and decoys.
You can use the ``--concat`` flag so that the output is a concatenated database
containing both target and decoy proteins:

.. code-block:: sh

    # If you want to output a concatenated fasta of targets and decoys:
    pydecoys targets.fasta -o targets_and_decoys.fasta --concat
