# Copyright © 2026 Bruno Maestri A Becker
#
# This file is part of Decoys.
#
# Decoys is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# Decoys is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# Decoys. If not, see <https://www.gnu.org/licenses/>.


import typing as _t

from Bio.SeqRecord import SeqRecord

from .DecoyStrategy import (
    reverse,
    reverse_keep_n,
    reverse_keep_c,
    reverse_keep_term,
    shuffle,
    PseudoReverseRule,
    PseudoShuffleRule,
    DecoyGenerator,
)


pseudoreverse_trypsin = PseudoReverseRule("KR", nocut="P")
pseudoreverse_stricttrypsin = PseudoReverseRule("KR")
pseudoreverse_argc = PseudoReverseRule("R", nocut="P")
pseudoreverse_aspn = PseudoReverseRule("D", sense="N")
pseudoreverse_chymo = PseudoReverseRule("FLWY", nocut="P")
pseudoreverse_gluc = PseudoReverseRule("DE", nocut="P")
pseudoreverse_lysc = PseudoReverseRule("K", nocut="P")
pseudoreverse_lysn = PseudoReverseRule("K", sense="N")

pseudoshuffle_trypsin = PseudoShuffleRule("KR", nocut="P")
pseudoshuffle_stricttrypsin = PseudoShuffleRule("KR")
pseudoshuffle_argc = PseudoShuffleRule("R", nocut="P")
pseudoshuffle_aspn = PseudoShuffleRule("D", sense="N")
pseudoshuffle_chymo = PseudoShuffleRule("FLWY", nocut="P")
pseudoshuffle_gluc = PseudoShuffleRule("DE", nocut="P")
pseudoshuffle_lysc = PseudoShuffleRule("K", nocut="P")
pseudoshuffle_lysn = PseudoShuffleRule("K", sense="N")


_decoy_strategy: dict[str, DecoyGenerator] = {
    "reverse": reverse,
    "reverse-keepn": reverse_keep_n,
    "reverse-keepc": reverse_keep_c,
    "reverse-keepterm": reverse_keep_term,
    "shuffle": shuffle,
    "pseudoreverse-trypsin": pseudoreverse_trypsin,
    "pseudoreverse-stricttrypsin": pseudoreverse_stricttrypsin,
    "pseudoreverse-argc": pseudoreverse_argc,
    "pseudoreverse-aspn": pseudoreverse_aspn,
    "pseudoreverse-chymo": pseudoreverse_chymo,
    "pseudoreverse-gluc": pseudoreverse_gluc,
    "pseudoreverse-lysc": pseudoreverse_lysc,
    "pseudoreverse-lysn": pseudoreverse_lysn,
    "pseudoshuffle-trypsin": pseudoshuffle_trypsin,
    "pseudoshuffle-stricttrypsin": pseudoshuffle_stricttrypsin,
    "pseudoshuffle-argc": pseudoreverse_argc,
    "pseudoshuffle-aspn": pseudoshuffle_aspn,
    "pseudoshuffle-chymo": pseudoshuffle_chymo,
    "pseudoshuffle-gluc": pseudoshuffle_gluc,
    "pseudoshuffle-lysc": pseudoshuffle_lysc,
    "pseudoshuffle-lysn": pseudoshuffle_lysn,
}


def generate(
    sequences: _t.Iterable[SeqRecord] | SeqRecord,
    strategy: str,
    decoy_tag: str = 'decoy_',
    prefix: bool = True,
) -> _t.Generator[SeqRecord, None, None]:
    if not isinstance(strategy, str):
        raise TypeError("Need a string for the decoy strategy (lower case)")
    if not strategy:
        raise ValueError("Strategy required (lower case string)")
    if not strategy.islower():
        raise ValueError(f"Strategy string '{strategy}' should be lower case")

    if not isinstance(decoy_tag, str):
        raise TypeError("Need a string for the decoy tag")

    if isinstance(sequences, SeqRecord):
        sequences = [sequences]

    decoy_generator = _decoy_strategy.get(strategy)

    if decoy_generator is None:
        raise ValueError(f"Unknown strategy: '{strategy}'")

    for sequence in sequences:
        if sequence.seq is None:
            raise ValueError(f"Seq not present for SeqRecord '{sequence.id}'")

        id = sequence.id if sequence.id else ""
        id = decoy_tag + id if prefix else id + decoy_tag
        seq = decoy_generator(sequence.seq)

        yield SeqRecord(seq, id, description="")


def register(strategy: str, decoy_generator_fn: DecoyGenerator) -> None:
    if not isinstance(strategy, str):
        raise TypeError("Need a string for the decoy strategy (lower case)")
    if not strategy:
        raise ValueError("Strategy required (lower case string)")
    if not strategy.islower():
        raise ValueError(f"Strategy string '{strategy}' should be lower case")

    if strategy in _decoy_strategy:
        raise ValueError(f"Strategy '{strategy}' already exists")

    _decoy_strategy[strategy] = decoy_generator_fn
