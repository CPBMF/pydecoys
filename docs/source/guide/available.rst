.. _available-strategies:

Available strategies
====================
Available strategies for protein generation, that can be used from the CLI or
API:

=================== =================================================================================
str key             Description
=================== =================================================================================
reverse             Reverse the aminoacid order of the target protein
shuffle             Shuffle the target protein aminoacids
reversepep-<enzyme> Reverse the aminoacid order of the enzymatic fragments, except the cleavage sites
shufflepep-<enzyme> Shuffle the aminoacids of the enzymatic fragments, except the cleavage sites
=================== =================================================================================

Terminal-aa options
-------------------
If you want to preserve terminal aminoacids from the target protein, you can
choose one of the following:

========== =================================
str key    Description
========== =================================
\-keepn    Preserve the N-terminal aminoacid
\-keepc    Preserve the C-terminal aminoacid
\-keepterm Preserve both terminal aminoacids
========== =================================

Enzymes
-------
For enzyme-specific strategies, the followin enzymes are available:

============== ============= ==== ===== =====
Enzyme         str key       cut  nocut sense
============== ============= ==== ===== =====
Trypsin        trypsin       KR   P     C
Strict trypsin stricttrypsin KR   \-    C
Arg-C          argc          R    P     C
Asp-N          aspn          D    \-    N
Chymotrypsin   chymo         FLWY P     C
Glu-C          gluc          DE   P     C
Lys-C          lysc          K    P     C
Lys-N          lysn          K    \-    N
Pepsin-A       pepsina       FL   \-    C
CNBr           cnbr          M    \-    C
============== ============= ==== ===== =====

Examples
--------
For example, to reverse tryptic peptides:

>>> seq = 'QSYKPTRTHQ'
>>> seq_as_decoy(seq, 'reversepep-trypsin')
'TPKYSQRQHT'

If you want to reverse the whole protein, except the N-terminal:

>>> seq = 'QSYKPTRTHQ'
>>> seq_as_decoy(seq, 'reverse-keepn')
'QQHTRTPKYS'

If you want so simply shuffle it:

>>> seq = 'QSYKPTRTHQ'
>>> seq_as_decoy(seq, 'suffle')
'RSTTHKQQYP'
