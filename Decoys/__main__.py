# Copyright © 2026 Bruno Maestri A Becker
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


import sys

from Bio import SeqIO

from . import generate


def main(input: str, output: str):
    targets = SeqIO.parse(input, 'fasta')
    decoys = generate(targets, 'pseudoreverse-stricttrypsin')
    SeqIO.write(decoys, output, 'fasta')


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
