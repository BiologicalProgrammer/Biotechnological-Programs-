# DNA Translation and Motif Finder

This project provides two interfaces:

1. A **CLI tool** for translating DNA and finding motifs.
2. A **web app** (FastAPI + HTML/JS) for browser-based analysis.

## Features

- DNA normalization/validation (`A`, `C`, `G`, `T`, `N`)
- Translation in frame `1`, `2`, or `3`
- Reverse-complement option
- Optional stop-at-first-stop behavior
- Motif search with overlapping matches
- IUPAC motif ambiguity support (e.g., `AAR`, `GGN`)

## CLI usage

```bash
python3 dna_translation_motif_finder.py "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG" --motif "GGN"
```

### CLI options

- `--frame {1,2,3}`: Select translation reading frame (default `1`)
- `--reverse-complement`: Analyze reverse complement instead of input sequence
- `--stop-at-stop`: Stop translation at first stop codon
- `--motif MOTIF`: Find motif and print all match positions (0-based)

If you run from **Python IDLE** with no command-line arguments, an interactive menu collects:

- DNA sequence
- reading frame
- reverse-complement choice
- stop-at-stop choice
- optional motif

## Website usage

### 1) Install dependencies

```bash
python3 -m pip install -r requirements.txt
```

### 2) Run the web app

```bash
uvicorn app:app --reload
```

### 3) Open in browser

Visit: `http://127.0.0.1:8000`

Use the form to submit sequence/options and view translation + motif results.

## API endpoint

`POST /api/analyze`

Example JSON payload:

```json
{
  "sequence": "ATGAAATGA",
  "frame": 1,
  "reverse_complement": false,
  "stop_at_stop": true,
  "motif": "AAR"
}
```

Example response fields:

- `dna`
- `protein`
- `motif_positions`
- plus selected option flags
