# Separate .pyi file to keep logic code clean while providing full type hinting

from __future__ import annotations

from typing import (
    Iterable,
    Literal,
    Protocol,
    TYPE_CHECKING,
    overload,
    runtime_checkable
)


if TYPE_CHECKING:
    from Bio.Seq import Seq, MutableSeq


type SeqLike = 'str | Seq | MutableSeq'


class DecoyGenerator(Protocol):
    @overload
    def __call__(self, sequence: Seq | MutableSeq) -> Seq:
        ...

    @overload
    def __call__(self, sequence: str) -> str:
        ...

    def __call__(self, sequence: SeqLike) -> SeqLike:
        ...


@runtime_checkable
class ContextfulGenerator(DecoyGenerator, Protocol):
    def learn_context(
        self,
        sequences: Iterable[SeqLike]
    ) -> None:
        ...


@overload
def reverse(sequence: Seq | MutableSeq) -> Seq: ...


@overload
def reverse(sequence: str) -> str: ...


@overload
def reverse_keep_n(sequence: Seq | MutableSeq) -> Seq: ...


@overload
def reverse_keep_n(sequence: str) -> str: ...


@overload
def reverse_keep_c(sequence: Seq | MutableSeq) -> Seq: ...


@overload
def reverse_keep_c(sequence: str) -> str: ...


@overload
def reverse_keep_term(sequence: Seq | MutableSeq) -> Seq: ...


@overload
def reverse_keep_term(sequence: str) -> str: ...


@overload
def shuffle(sequence: Seq | MutableSeq) -> Seq: ...


@overload
def shuffle(sequence: str) -> str: ...


@overload
def shuffle_keep_n(sequence: Seq | MutableSeq) -> Seq: ...


@overload
def shuffle_keep_n(sequence: str) -> str: ...


@overload
def shuffle_keep_c(sequence: Seq | MutableSeq) -> Seq: ...


@overload
def shuffle_keep_c(sequence: str) -> str: ...


@overload
def shuffle_keep_term(sequence: Seq | MutableSeq) -> Seq: ...


@overload
def shuffle_keep_term(sequence: str) -> str: ...


class PseudoReverseRule:
    cut: str
    nocut: str
    sense: Literal['C', 'N']

    def __init__(
        self,
        cut: str,
        sense: Literal['C', 'N'] = 'C',
        nocut: str | None = None
    ) -> None:
        ...

    @overload
    def __call__(self, sequence: Seq | MutableSeq) -> Seq: ...

    @overload
    def __call__(self, sequence: str) -> str: ...

    def decoy_from_str(self, sequence: str) -> str: ...

    def decoy_from_Seq(self, sequence: Seq | MutableSeq) -> Seq: ...


class PseudoShuffleRule:
    cut: str
    nocut: str
    sense: Literal['C', 'N']

    def __init__(
        self,
        cut: str,
        sense: Literal['C', 'N'] = 'C',
        nocut: str | None = None
    ) -> None:
        ...

    @overload
    def __call__(self, sequence: Seq | MutableSeq) -> Seq: ...

    @overload
    def __call__(self, sequence: str) -> str: ...

    def decoy_from_str(self, sequence: str) -> str: ...

    def decoy_from_Seq(self, sequence: Seq | MutableSeq) -> Seq: ...

    def _shuffle(self, frag: str) -> str: ...


pseudoreverse_trypsin: DecoyGenerator
pseudoreverse_stricttrypsin: DecoyGenerator
pseudoreverse_argc: DecoyGenerator
pseudoreverse_aspn: DecoyGenerator
pseudoreverse_chymo: DecoyGenerator
pseudoreverse_gluc: DecoyGenerator
pseudoreverse_lysc: DecoyGenerator
pseudoreverse_lysn: DecoyGenerator
pseudoreverse_stricttrypsin_keepn: DecoyGenerator

pseudoshuffle_trypsin: DecoyGenerator
pseudoshuffle_stricttrypsin: DecoyGenerator
pseudoshuffle_argc: DecoyGenerator
pseudoshuffle_aspn: DecoyGenerator
pseudoshuffle_chymo: DecoyGenerator
pseudoshuffle_gluc: DecoyGenerator
pseudoshuffle_lysc: DecoyGenerator
pseudoshuffle_lysn: DecoyGenerator
pseudoshuffle_stricttrypsin_keepn: DecoyGenerator
