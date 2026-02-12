# DNA Translation and Motif Finder

A lightweight Python command-line tool that:

1. Translates DNA sequences into protein sequences.
2. Finds DNA motifs (including IUPAC ambiguity codes) in a sequence.

## Usage

```bash
python3 dna_translation_motif_finder.py "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG" --motif "GGN"
```

### Options

- `--frame {1,2,3}`: Select translation reading frame (default `1`)
- `--reverse-complement`: Analyze reverse complement instead of input sequence
- `--stop-at-stop`: Stop translation at first stop codon
- `--motif MOTIF`: Find motif and print all match positions (0-based)


If you run the script from **Python IDLE** (Run Module) without command-line arguments, it will prompt:

```text
Enter DNA sequence:
```

## Example

```bash
python3 dna_translation_motif_finder.py "ATGAAATGA" --stop-at-stop --motif "AAR"
```

Expected output shape:

- DNA sequence used for analysis
- Translated protein sequence
- Motif hit positions (if motif supplied)
