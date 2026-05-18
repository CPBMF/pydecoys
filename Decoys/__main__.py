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


from pathlib import Path

from Bio import SeqIO

from . import generate


def main(input: Path, output: Path, type: str):
    targets = SeqIO.parse(input, 'fasta')
    decoys = generate(targets, type)
    SeqIO.write(decoys, output, 'fasta')


if __name__ == '__main__':
    import argparse

    from . import _decoy_strategy

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('input', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument(
        '--type', '-t',
        choices=_decoy_strategy.keys(),
        default='reverse'
    )

    args = parser.parse_args()

    main(args.input, args.output, args.type)
