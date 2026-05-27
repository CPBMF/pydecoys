# Copyright (C) 2026 CPBMF and INCT-TB, PUCRS, Porto Alegre, Brazil
# Copyright (C) 2026 Bruno Maestri A Becker
#
# This file is part of PyDecoys.
#
# PyDecoys is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# PyDecoys is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# PyDecoys. If not, see <https://www.gnu.org/licenses/>.

"""
Provides:

    #. A CLI app to easily generate a decoy fasta file from a target fasta file
    #. A comprehensive Python API with IO, Generator and single data functions
    #. Easy implementation of custom decoy strategies
    #. Full integration with Biopython

Check the documentation for guides and the full API reference:
<https://cpbmf.github.io/pydecoys/>.

Available functionality
-----------------------
:func:`from_fasta`
    Return a decoy generator from a target fasta (requires `Biopython`)
:func:`to_fasta`
    Write a decoy fasta from a target protein set or fasta (requires `Biopython`)
:func:`from_SeqRecords`
    Return a decoy generator from a target set of `SeqRecord` (requires `Biopython`)
:func:`from_seqs`
    Return a decoy generator from a target set of :type:`SeqLike`
:func:`from_tuples`
    Return a decoy generator from a target set of `tuple`
:func:`SeqRecord_as_decoy`
    Return a decoy `SeqRecord` from a given `SeqRecord` (requires `Biopython`)
:func:`seq_as_decoy`
    Return a decoy :type:`SeqLike` from a given :type:`SeqLike`
:func:`tuple_as_decoy`
    Return a decoy `tuple` from a given `tuple`
:func:`register`
    Register a custom decoy strategy on the API
:mod:`pydecoys.strategies`
    API for custom decoy generation strategies

Avaliable decoy strategies
--------------------------
Each decoy strategy is specified by a lowercase string. New strategies can be
created by following the :type:`strategies.DecoyGenerator` function
signature or the rest of the :mod:`pydecoys.strategies` API.

- reverse:                         Reverse protein
- reverse-keepn:                   Reverse protein, except N-terminal aa
- reverse-keepc:                   Reverse protein, except C-terminal aa
- reverse-keepterm:                Reserse protein, except terminal aas
- shuffle:                         Shuffle protein
- shuffle-keepn:                   Shuffle protein, except N-terminal aa
- shuffle-keepc:                   Shuffle protein, except C-terminal aa
- shuffle-keepterm:                Shuffle protein, except terminal aas
- randomize:                       Randomize protein
- randomize-keepn:                 Randomize protein, except N-terminal aa
- randomize-keepc:                 Randomize protein, except C-terminal aa
- randomize-keepterm:              Randomize protein, except terminal aas
- reversepep-[enzyme]:             Pseudo-reverse enzymatic fragments
- reversepep-[enzyme]-keepn:       Pseudo-reverse enzymatic fragments, except N-terminal aa
- reversepep-[enzyme]-keepc:       Pseudo-reverse enzymatic fragments, except C-terminal aa
- reversepep-[enzyme]-keepterm:    Pseudo-reverse enzymatic fragments, except terminal aas
- shufflepep-[enzyme]:             Pseudo-shuffle enzymatic fragments
- shufflepep-[enzyme]-keepn:       Pseudo-shuffle enzymatic fragments, except N-terminal aa
- shufflepep-[enzyme]-keepc:       Pseudo-shuffle enzymatic fragments, except C-terminal aa
- shufflepep-[enzyme]-keepterm:    Pseudo-shuffle enzymatic fragments, except terminal aas
- randomizepep-[enzyme]:           Pseudo-randomize enzymatic fragments
- randomizepep-[enzyme]-keepn:     Pseudo-randomize enzymatic fragments, except N-terminal aa
- randomizepep-[enzyme]-keepc:     Pseudo-randomize enzymatic fragments, except C-terminal aa
- randomizepep-[enzyme]-keepterm:  Pseudo-randomize enzymatic fragments, except terminal aas

You can check the full list of proteases and how to add more at the
documentation.
"""  # noqa: W505, E501

from importlib.metadata import version


__version__ = version('pydecoys')

__all__ = [
    'from_fasta',
    'to_fasta',
    'from_SeqRecords',
    'from_seqs',
    'from_tuples',
    'SeqRecord_as_decoy',
    'seq_as_decoy',
    'tuple_as_decoy',
    'register',
    'strategies',
]


from pydecoys import strategies
from pydecoys._pydecoys import (
    from_fasta,
    to_fasta,
    from_SeqRecords,
    from_seqs,
    from_tuples,
    SeqRecord_as_decoy,
    seq_as_decoy,
    tuple_as_decoy,
    register,
)
