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

"""CLI entrypoint for PyDecoys."""

from Bio import SeqIO

from . import from_SeqRecords


def main() -> None:
    import argparse

    from . import __version__
    from ._pydecoys import _decoy_strategy

    parser = argparse.ArgumentParser(
        prog="pydecoys",
        description="Proteomics decoy utilities for Python. Generate decoy sequences from a fasta file."  # noqa: E501
    )

    parser.add_argument(
        '-v', '--version',
        action='version',
        version=f'%(prog)s {__version__} (GPL-3.0-or-later)'
    )

    parser.add_argument(
        'input',
        type=argparse.FileType(),
        default='-',
        nargs='?',
        help='input file (default: `stdin`)'
    )
    parser.add_argument(
        'output',
        type=argparse.FileType('w'),
        default='-',
        nargs='?',
        help='output file (default: `stdout`)'
    )
    parser.add_argument(
        '-s', '--strategy',
        choices=_decoy_strategy.keys(),
        default='reverse',
        help='decoy generation strategy (default: \'%(default)s\')'
    )

    args = parser.parse_args()

    targets = SeqIO.parse(args.input, 'fasta')
    decoys = from_SeqRecords(targets, args.strategy)
    SeqIO.write(decoys, args.output, 'fasta')


if __name__ == '__main__':
    main()
