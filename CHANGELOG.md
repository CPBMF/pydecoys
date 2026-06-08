# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- A `fuses` factory function that creates DecoyGenerators that fuse targets and
  decoys.
- A `-fuse` option to available strategies to prepend targets to decoys.
- Markov-chain models for generating decoy proteins from target proteins and
  decoy enzymatic fragments from target enzymatic fragments.

### Changed

- Refactor strategy registration logic at module initialization so it's cleaner
  and easier to mantain.

### Fixed

- Stop `_builtins_test` from spilling state to other tests.

## [0.3.0] - 2026-06-04

### Added

- Support to non-standard aminoacids (O, U) and special codes (B, J, Z, X):
  - O: Pyrrolysine (Pyl)
  - U: Selenocysteine (Sec)
  - B: Aspartic acid (D) or Asparagine (N)
  - J: Leucine (L) or Isoleucine (I)
  - Z: Glutamic acid (E) or Glutamine (Q)
  - X: Unspecified aminoacid
- New cleavage agents from [PSI MS Ontology](https://github.com/HUPO-PSI/psi-ms-CV):
  - AlphaLP: `alphalp`
  - Asp-N_ambic: `aspnambic`
  - Formic_acid: `formicacid`
  - leukocyte elastase: `elastase`
  - Lys-C/P: `lyscp`
  - proline endopeptidase: `proc`
  - TrypChymo: `trypchymo`
  - Tryp-N: `trypn`
  - 2-iodobenzoate: `2iodobenzoate`
  - V8-E: `v8e`
- Regex-based instantiation of `EnzymeSpecificGenerator`. The regex pattern
  must match only the cleavage sites, and must capture them.
- A `pattern` attribute to `EnzymeSpecificGenerator` that returns the regex
  pattern used by the instance.
- Package type-hinting.
- Case-insensitivity to all decoy strategies.

### Changed

- Follow the [PSI MS Ontology](https://github.com/HUPO-PSI/psi-ms-CV) cleavage
  agent definitions.
- Add 'B' residue to `cut` specifications of Asp-N.
- Add 'BZ' residues to `cut` specifications of Glu-C (now V8-DE).
- Rename old Glu-C specification to V8-DE (the key `gluc` was changed to
  `v8de`). A new Glu-C (key `gluc`) specification was added with a slightly
  different regex pattern.
- Rename old Strict Trypsin specification to Trypsin/P (the key `stricttrypsin`
  was changed to `trypsinp`).
- Rename old Pepsin-A specification to PepsinA (the key is unchanged).
- Move old `__init__` method of `EnzymeSpecificGenerator` to a new class method
  `from_enzyme`. The new `from_enzyme` method also accepts a `nocut_n`
  parameter that determines an aminoacid that stops cleavage when at the
  N-terminal bond of the cut aminoacid.

### Fixed

- Fix the `nocut` parameter of `EnzymeSpecificGenerator` so that the regex will
  correctly identify any of the aminoacids given (before it treated the
  parameter as a literal string).
- Fix the type-hinting of the factory functions from `pydecoys.strategies`.

### Removed

- The `cut`, `nocut` and `sense` attributes of `EnzymeSpecificGenerator`.

## [0.2.1] - 2026-06-02

### Fixed

- Fix the version links in this CHANGELOG.
- Padronize file references in this CHANGELOG.
- Remove the redundant `[tool.poetry]` table from the pyproject file.
- Fix the installation guide.

## [0.2.0] - 2026-06-01

### Added

- Explicit error messages when trying to use functionality dependant on
  Biopython without Biopython installed.
- A `reset` method and `is_set` attribute to `strategies.ContextfulGenerator`
  protocol.
- A `get_contextualized_strategy` function that returns a contextualized
  strategy from a target dataset and strategy key.
- This CHANGELOG file.

### Changed

- Expand test coverage.
- Update how the main API handles ontext-dependant strategies (the API now
  sets, uses and then resets unset context-dependant strategies, but never
  modifies set instances; using a str-key guarantees an unset strategy).
- Improve documentation.

### Fixed

- Make `register` raise a `TypeError` when receiving a non callable.
- Make `seq_as_decoy` raise a `ValueError` upon an empty `str`.
- Padronize exception messages.

## [0.1.0] - 2026-05-27

### Added

- This project to GitHub to serve as a simple decoy generation tool for both
  CLI and Python workflows.
- A CLI application to generate decoy fasta files from target fasta files.
- An API to generate decoy protein sequences from target protein sequences in
  Python applications.
- The following decoy generation strategies:
  - reverse: Reverse the aminoacid order of the target protein
  - shuffle: Shuffle the target protein aminoacids
  - randomize: Randomize the target protein aminoacids using the proportions
    of the target database
  - reversepep: Reverse the aminoacid order of the enzymatic fragments, except
    the cleavage sites
  - shufflepep: Shuffle the aminoacids of the enzymatic fragments, except the
    cleavage sites
  - randomize: Randomize the aminoacids of the enzymatic fragments, except the
    cleavage sites, using the proportions of the target database (cleavage
    sites aren't counted)
- The following proteases:
  - Trypsin
  - Strict trypsin
  - Arg-C
  - Asp-N
  - Chymotrypsin
  - Glu-C
  - Lys-C
  - Lys-N
  - Pepsin-A
  - CNBr
- The following terminal-aminoacid options:
  - keepn: Preserve the N-terminal aminoacid
  - keepc: Preserve the C-terminal aminoacid
  - keepterm: Preserve the N- and C-terminal aminoacids
- The [documentation](https://cpbmf.github.io/pydecoys/).
- The [README](https://github.com/CPBMF/pydecoys/blob/main/README.md).

[unreleased]: https://github.com/CPBMF/pydecoys/compare/v0.3.0...HEAD
[0.3.0]: https:/github.com/CPBMF/pydecoys/compare/v0.2.1...v0.3.0
[0.2.1]: https:/github.com/CPBMF/pydecoys/compare/v0.2.0...v0.2.1
[0.2.0]: https:/github.com/CPBMF/pydecoys/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/CPBMF/pydecoys/releases/tag/v0.1.0

<!-- generated by git-cliff -->
