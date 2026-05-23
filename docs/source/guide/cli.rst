Running the CLI
===============

Running `pydecoys` is simple:

.. code-block:: sh

    # To generate a reversed decoy set from a given `targets.fasta`
    pydecoys targets.fasta -o rev.fasta

    # Change the strategy with the `-s/--strategy` flag
    pydecoys targets.fasta -o shuf.fasta -s shuffle

    # The input defaults to stdin
    cat database_1.fasta database_2.fasta | pydecoys -o rev_1_2.fasta

    # And the output, to stdout
    pydecoys targets.fasta | grep some_seq

    # If you want to output a concatenated fasta of targets and decoys:
    pydecoys --concat targets.fasta -o targets_and_decoys.fasta

.. note::
    See the full list of available strategies at :ref:`available-strategies`.
    The full list of CLI options is at :ref:`cli`
