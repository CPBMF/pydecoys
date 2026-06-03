.. _available-strategies:

Available strategies
====================
Available strategies for protein generation, that can be used from the CLI or
API:

=================== ==================================================================================
str key             Description
=================== ==================================================================================
reverse             Reverse the aminoacid order of the target protein
shuffle             Shuffle the target protein aminoacids
randomize           Randomize the target protein aminoacids using the proportions of the target tabase
reversepep-<enzyme> Reverse the aminoacid order of the enzymatic fragments, except the cleavage sites
shufflepep-<enzyme> Shuffle the aminoacids of the enzymatic fragments, except the cleavage sites
randomize-<enzyme>  Randomize the aminoacids of the enzymatic fragments, except the cleavage sites
=================== ==================================================================================

Enzymes
-------
For enzyme-specific strategies, the followin enzymes are available:

============== =============== ==== ===== ======= ===== ================
Enzyme         str key         cut  nocut nocut_n sense pattern
============== =============== ==== ===== ======= ===== ================
Trypsin        trypsin         KR   P     \-      C     r'([KR])(?!P)'
Trypsin/P      trypsinp        KR   \-    \-      C     r'([KR])'
Arg-C          argc            R    P     \-      C     r'(R)(?!P)'
Asp-N          aspn            BD   \-    \-      N     r'([BD])'
Chymotrypsin   chymo           FYWL P     \-      C     r'([FYWL])(?!P)'
V8-DE          v8de            BDEZ P     \-      C     r'([BDEZ])(?!P)'
Lys-C          lysc            K    P     \-      C     r'(K)(?!P)'
Lys-N          lysn            K    \-    \-      N     r'(K)'
PepsinA        pepsina         FL   \-    \-      C     r'([FL])'
CNBr           cnbr            M    \-    \-      C     r'(M)'
============== =============== ==== ===== ======= ===== ================

Terminal-aa options
-------------------
If you want to preserve terminal aminoacids from the target protein, you can
choose one of the following:

========== =========================================
str key    Description
========== =========================================
\-keepn    Preserve the N-terminal aminoacid
\-keepc    Preserve the C-terminal aminoacid
\-keepterm Preserve the N- and C-terminal aminoacids
========== =========================================

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

You can combine both:

>>> seq = 'QSYKPTRTHQ'
>>> seq_as_decoy(seq, 'reversepep-trypsin-keepn')
'QTPKYSRQHT'

If you want so simply shuffle it:

>>> seq = 'QSYKPTRTHQ'
>>> seq_as_decoy(seq, 'suffle')
'RSTTHKQQYP'
