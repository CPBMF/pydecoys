.. _cli:

CLI Reference
=============

Dependencies
------------
- Requires Python 3.12 or superior.
- Requires Biopython. This should be automatically installed. If not,
  run `pip install biopython`.


CLI options
-----------
`input`
^^^^^^^
input file (default: `stdin`)

`-o/--output`
^^^^^^^^^^^^^
output filename (default: `stdout`)

`-s/--strategy`
^^^^^^^^^^^^^^^
decoy generation strategy (default: `'reverse'`)

`--decoy-tag`
^^^^^^^^^^^^^
decoy tag identifier (default: `'decoy_'`)

`--prefix | --suffix`
^^^^^^^^^^^^^^^^^^^^^
add `decoy_tag` as prefix (default) or suffix

`--concat`
^^^^^^^^^^
if specified, output fasta will have both target and decoy entries

`-h/--help`
^^^^^^^^^^^
show help message and exit

`-v/--version`
^^^^^^^^^^^^^^
show version and exit
