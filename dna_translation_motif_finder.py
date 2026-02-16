#!/usr/bin/env python3
"""DNA translation and motif finder tool."""

from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass
from typing import List, Sequence

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

IUPAC_TO_REGEX = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "R": "[AG]",
    "Y": "[CT]",
    "S": "[GC]",
    "W": "[AT]",
    "K": "[GT]",
    "M": "[AC]",
    "B": "[CGT]",
    "D": "[AGT]",
    "H": "[ACT]",
    "V": "[ACG]",
    "N": "[ACGT]",
}


@dataclass
class InteractiveOptions:
    sequence: str
    frame: int = 1
    reverse_complement: bool = False
    stop_at_stop: bool = False
    motif: str | None = None


def normalize_dna(sequence: str) -> str:
    sequence = "".join(sequence.split()).upper()
    invalid = sorted({base for base in sequence if base not in {"A", "C", "G", "T", "N"}})
    if invalid:
        raise ValueError(f"DNA contains invalid bases: {', '.join(invalid)}")
    return sequence


def reverse_complement(sequence: str) -> str:
    table = str.maketrans("ACGTN", "TGCAN")
    return sequence.translate(table)[::-1]


def translate_dna(sequence: str, frame: int = 1, stop_at_stop: bool = False) -> str:
    if frame not in {1, 2, 3}:
        raise ValueError("Frame must be 1, 2, or 3.")

    sequence = normalize_dna(sequence)
    start = frame - 1
    amino_acids: List[str] = []

    for i in range(start, len(sequence) - 2, 3):
        codon = sequence[i : i + 3]
        aa = CODON_TABLE.get(codon, "X")
        if stop_at_stop and aa == "*":
            break
        amino_acids.append(aa)

    return "".join(amino_acids)


def motif_to_regex(motif: str) -> str:
    motif = motif.upper().strip()
    if not motif:
        raise ValueError("Motif cannot be empty.")

    try:
        return "".join(IUPAC_TO_REGEX[char] for char in motif)
    except KeyError as exc:
        raise ValueError(f"Unsupported motif symbol: {exc.args[0]}") from exc


def find_motifs(sequence: str, motif: str) -> Sequence[int]:
    sequence = normalize_dna(sequence)
    regex = motif_to_regex(motif)

    pattern = re.compile(f"(?=({regex}))")
    return [match.start() for match in pattern.finditer(sequence)]


def prompt_for_sequence() -> str | None:
    """Prompt for a DNA sequence."""
    try:
        value = input("Enter DNA sequence: ").strip()
    except EOFError:
        return None
    return value or None


def prompt_yes_no(question: str, default: bool = False) -> bool:
    """Prompt for yes/no and return bool."""
    suffix = " [Y/n]: " if default else " [y/N]: "

    while True:
        try:
            raw = input(question + suffix).strip().lower()
        except EOFError:
            return default

        if not raw:
            return default
        if raw in {"y", "yes"}:
            return True
        if raw in {"n", "no"}:
            return False

        print("Please answer y or n.")


def prompt_frame(default: int = 1) -> int:
    """Prompt for translation frame 1-3."""
    while True:
        try:
            raw = input(f"Reading frame [1/2/3] (default {default}): ").strip()
        except EOFError:
            return default

        if not raw:
            return default
        if raw in {"1", "2", "3"}:
            return int(raw)

        print("Please enter 1, 2, or 3.")


def interactive_menu() -> InteractiveOptions | None:
    """Collect options interactively (ideal for IDLE users)."""
    print("Interactive mode: enter options (press Enter to accept defaults).")
    sequence = prompt_for_sequence()
    if not sequence:
        return None

    frame = prompt_frame(default=1)
    reverse = prompt_yes_no("Use reverse complement?", default=False)
    stop_at_stop = prompt_yes_no("Stop translation at first stop codon?", default=False)

    try:
        motif_raw = input("Motif (optional, supports IUPAC codes): ").strip()
    except EOFError:
        motif_raw = ""

    return InteractiveOptions(
        sequence=sequence,
        frame=frame,
        reverse_complement=reverse,
        stop_at_stop=stop_at_stop,
        motif=motif_raw or None,
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Translate DNA and find motifs.")
    parser.add_argument("sequence", nargs="?", help="DNA sequence (A,C,G,T,N)")
    parser.add_argument("--frame", type=int, choices=[1, 2, 3], default=1, help="Reading frame")
    parser.add_argument("--reverse-complement", action="store_true", help="Use reverse complement before analysis")
    parser.add_argument("--stop-at-stop", action="store_true", help="Stop translation at the first stop codon")
    parser.add_argument("--motif", help="Motif to find (supports IUPAC ambiguity codes)")
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    if not args.sequence:
        if len(sys.argv) == 1:
            interactive = interactive_menu()
            if interactive:
                args.sequence = interactive.sequence
                args.frame = interactive.frame
                args.reverse_complement = interactive.reverse_complement
                args.stop_at_stop = interactive.stop_at_stop
                args.motif = interactive.motif

        if not args.sequence:
            parser.print_usage()
            print("error: sequence is required (example: ATGAAATGA)", file=sys.stderr)
            return 2

    try:
        sequence = normalize_dna(args.sequence)
        if args.reverse_complement:
            sequence = reverse_complement(sequence)

        protein = translate_dna(sequence, frame=args.frame, stop_at_stop=args.stop_at_stop)
        print(f"DNA:      {sequence}")
        print(f"Protein:  {protein}")

        if args.motif:
            positions = find_motifs(sequence, args.motif)
            if positions:
                joined = ", ".join(str(p) for p in positions)
                print(f"Motif '{args.motif.upper()}' found at positions (0-based): {joined}")
            else:
                print(f"Motif '{args.motif.upper()}' not found.")

    except ValueError as error:
        parser.error(str(error))

    return 0


if __name__ == "__main__":
    if "idlelib" in sys.modules:
        main()
    else:
        raise SystemExit(main())
