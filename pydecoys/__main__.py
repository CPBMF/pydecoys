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

"""CLI entrypoint for PyDecoys."""

from pydecoys import to_fasta


def main(arg_list: list[str] | None = None) -> None:
    import argparse

    from pydecoys import __version__
    from pydecoys._builtins import decoy_strategy

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
        '-o', '--output',
        type=argparse.FileType('w'),
        default='-',
        help='output filename (default: `stdout`)'
    )
    parser.add_argument(
        '-s', '--strategy',
        choices=decoy_strategy.keys(),
        default='reverse',
        help='decoy generation strategy (default: \'%(default)s\')'
    )
    parser.add_argument(
        '-t', '--decoy-tag',
        default='decoy_',
        help='decoy tag identifier (default: \'%(default)s\')'
    )

    prefix = parser.add_mutually_exclusive_group()
    prefix.add_argument(
        '--prefix',
        action='store_true',
        dest='prefix',
        help='add `decoy_tag` as prefix (default)'
    )
    prefix.add_argument(
        '--suffix',
        action='store_false',
        dest='prefix',
        help='add `decoy_tag` as suffix'
    )
    parser.set_defaults(prefix=True)

    parser.add_argument(
        '--concat',
        action='store_true',
        help='if specified, output fasta will have both target and decoy entries'
    )

    args = parser.parse_args(arg_list)

    to_fasta(
        args.input,
        args.output,
        args.strategy,
        args.decoy_tag,
        args.prefix,
        args.concat
    )


if __name__ == '__main__':
    main()
