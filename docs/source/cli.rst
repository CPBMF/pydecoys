CLI Reference
=============

Getting started
---------------
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

.. note::
    See the full list of available strategies at :ref:`available-strategies`.

List of options
---------------
The full list of available options:

`input`
^^^^^^^
input file (default: `stdin`)

`-o/--output`
^^^^^^^^^^^^^
output filename (default: `stdout`)

`-s/--strategy`
^^^^^^^^^^^^^^^
decoy generation strategy (default: `reverse`)

`-h/--help`
^^^^^^^^^^^
show help message and exit

`-v/--version`
^^^^^^^^^^^^^^
show version and exit
