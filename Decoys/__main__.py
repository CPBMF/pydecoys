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


from os import PathLike

from Bio import SeqIO

from . import from_SeqRecords, __version__


def main(
    input: str | PathLike[str],
    output: str | PathLike[str],
    strategy: str
) -> None:
    targets = SeqIO.parse(input, 'fasta')
    decoys = from_SeqRecords(targets, strategy)
    SeqIO.write(decoys, output, 'fasta')


if __name__ == '__main__':
    import argparse

    from . import _decoy_strategy

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('input')
    parser.add_argument('output')
    parser.add_argument(
        '--version', '-v',
        action='version',
        version=f'%(prog)s {__version__} (GPL-3.0-or-later)'
    )
    parser.add_argument(
        '--strategy', '-s',
        choices=_decoy_strategy.keys(),
        default='reverse'
    )

    args = parser.parse_args()

    main(args.input, args.output, args.strategy)
