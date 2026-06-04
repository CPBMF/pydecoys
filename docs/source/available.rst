.. _available-strategies:

Available strategies
====================

Available strategies for protein generation, that can be used from the CLI or
API:

.. list-table::
    :widths: 30 70
    :width: 100%
    :header-rows: 1

    * - str key
      - Description
    * - reverse
      - Reverse the aminoacid order of the target protein.
    * - shuffle
      - Shuffle the target protein's aminoacids.
    * - randomize
      - Randomize the target protein's aminoanicds. Context-dependat. Each
        aminoacid has a weighted likelihood based on the aminoacid proportions
        of a target dataset. Any character in the target dataset not in
        :data:`pydecoys.strategies.EXT_AMINOACIDS` is silently ignored.
    * - reversepep-<enzyme>
      - Reverse the aminoacid order of the enzymatic fragments, except the
        cleavage sites. The cleavage sites are preserved in-place.
    * - shufflepep-<enzyme>
      - Shufflle the enzymatic fragments' aminoacids in place, except the
        cleavage sites. The cleavage sites are preserved in-place.
    * - randomizepep-<enzyme>
      - Randomize the enzymatic fragments' aminoacids, except the cleavage
        sites. The cleavage sites are preserved in-place. Context-dependat.
        Each aminoacid has a weighted likelihood based on the aminoacid
        proportions of a target dataset. Any character in the target dataset
        not in :data:`pydecoys.strategies.EXT_AMINOACIDS` is silently ignored.
        The cleavage sites are also ignored in the proportion.

Enzymes
-------
For enzyme-specific strategies, the following enzymes are available:

.. list-table::
    :width: 100%
    :header-rows: 1

    * - Enzyme
      - str key
      - cut
      - nocut
      - nocut_n
      - sense
      - pattern
    * - AlphaLP
      - alphalp
      - TASV
      - \-
      - \-
      - C
      - ``r'([TASV])'``
    * - Arg-C
      - argc
      - R
      - P
      - \-
      - C
      - ``r'(R)(?!P)'``
    * - Asp-N
      - aspn
      - BD
      - \-
      - \-
      - N
      - ``r'([BD])'``
    * - Asp-N_ambic
      - aspnambic
      - DE
      - \-
      - \-
      - N
      - ``r'([DE])'``
    * - Chymotrypsin
      - chymo
      - FYWL
      - P
      - \-
      - C
      - ``r'([FYWL])(?!P)'``
    * - CNBr
      - cnbr
      - M
      - \-
      - \-
      - C
      - ``r'(M)'``
    * - Formic_acid
      - formicacid
      - D
      - \-
      - \-
      - both
      - ``r'(D)'``
    * - glutamyl endopeptidase
      - gluc
      - E
      - \-
      - E
      - C
      - ``r'(?<!E)(E)'``
    * - leukocyte elastase
      - elastase
      - ALIV
      - P
      - \-
      - C
      - ``r'([ALIV])(?!P)'``
    * - Lys-C
      - lysc
      - K
      - P
      - \-
      - C
      - ``r'(K)(?!P)'``
    * - Lys-C/P
      - lyscp
      - K
      - \-
      - \-
      - C
      - ``r'(K)'``
    * - Lys-N
      - lysn
      - K
      - \-
      - \-
      - N
      - ``r'(K)'``
    * - PepsinA
      - pepsina
      - FL
      - \-
      - \-
      - C
      - ``r'([FL])'``
    * - proline endopeptidase
      - proc
      - P
      - P
      - \-
      - C
      - ``r'([HKR]P)(?!P)'``
    * - Trypsin
      - trpysin
      - KR
      - P
      - \-
      - C
      - ``r'([KR])(?!P)'``
    * - Trypsin/P
      - trypsinp
      - KR
      - \-
      - \-
      - C
      - ``r'([KR])'``
    * - TrypChymo
      - trypchymo
      - FYWLKR
      - P
      - \-
      - C
      - ``r'([FYWLKR])(?!P)'``
    * - Tryp-N
      - trypn
      - KR
      - \-
      - \-
      - N
      - ``r'([KR])'``
    * - 2-iodobenzoate
      - 2iodobenzoate
      - W
      - \-
      - \-
      - C
      - ``r'(W)'``
    * - V8-DE
      - v8de
      - BDEZ
      - P
      - \-
      - C
      - ``r'([BDEZ])(?!P)'``
    * - V8-E
      - v8e
      - EZ
      - P
      - \-
      - C
      - ``r'([EZ])(?!P)'``

The enzyme specifications are taken from the `PSI MS Ontology`_.

Whenever PyDecoys refers to a 'cleavage site', what's meant is the capture
group of these regex patterns.

Since we need to actually capture the aminoacid(s) that induce cleavage to
preserve the cleavage sites in enzymatic-aware strategies, the patterns
presented here aren't the ones in `PSI MS Ontology`_. For paterns specifically
meant for splitting peptides, see the `PSI MS Ontology`_.

The patterns presented here should generate the same results, with the
difference that they capture the cleavage site aminoacids instead of
capturing the cleavage point itself. Another key difference is that, while
`PSI MS Ontology`_ regexes might not capture a cleavage point if the 'cut'
aminoacid is at either end of the protein, the paterns presented here will.

.. note::
    Some strategies presented here are redundant (they share the same regex
    pattern). They are separetely included in the case a strategy that needs to
    differentiate them is added (for example, a strategy that actually needs to
    consider cleavage sense), as well as for completeness sake.

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

.. _PSI MS Ontology: https://github.com/HUPO-PSI/psi-ms-CV
