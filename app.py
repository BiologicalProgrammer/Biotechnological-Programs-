from __future__ import annotations

from fastapi import FastAPI, HTTPException
from fastapi.responses import HTMLResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
from fastapi import Request
from pydantic import BaseModel

from dna_translation_motif_finder import (
    find_motifs,
    normalize_dna,
    reverse_complement,
    translate_dna,
)

app = FastAPI(title="DNA Translation and Motif Finder")
app.mount("/static", StaticFiles(directory="static"), name="static")
templates = Jinja2Templates(directory="templates")


class AnalyzeRequest(BaseModel):
    sequence: str
    frame: int = 1
    reverse_complement: bool = False
    stop_at_stop: bool = False
    motif: str | None = None


@app.get("/", response_class=HTMLResponse)
def index(request: Request) -> HTMLResponse:
    return templates.TemplateResponse("index.html", {"request": request})


@app.post("/api/analyze")
def analyze(payload: AnalyzeRequest) -> dict:
    try:
        sequence = normalize_dna(payload.sequence)
        if payload.reverse_complement:
            sequence = reverse_complement(sequence)

        protein = translate_dna(
            sequence=sequence,
            frame=payload.frame,
            stop_at_stop=payload.stop_at_stop,
        )

        positions: list[int] = []
        if payload.motif:
            positions = list(find_motifs(sequence, payload.motif))

        return {
            "dna": sequence,
            "protein": protein,
            "frame": payload.frame,
            "reverse_complement": payload.reverse_complement,
            "stop_at_stop": payload.stop_at_stop,
            "motif": payload.motif.upper() if payload.motif else None,
            "motif_positions": positions,
        }
    except ValueError as error:
        raise HTTPException(status_code=400, detail=str(error)) from error
