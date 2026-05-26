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

"""Internal module responsible for lazily importing Biopython, registering
the relevant names and namespacing utility functions related to the Biopython
interface.
"""

from typing import Iterable

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def SeqRecord_to_tuple(record: SeqRecord) -> tuple[str, str]:
    return (record.id if record.id else "", str(record.seq))


def tuple_to_SeqRecord(record: tuple[str, str]) -> SeqRecord:
    return SeqRecord(Seq(record[1]), record[0], name="", description="")


def iter_SeqRecord(seq: Iterable[SeqRecord] | SeqRecord) -> Iterable[SeqRecord]:
    if isinstance(seq, SeqRecord):
        return [seq]
    return seq


def str_to_Seq(sequence: str) -> Seq:
    return Seq(sequence)
