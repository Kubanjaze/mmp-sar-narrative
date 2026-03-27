# Phase 97 — MMP Analysis + Claude SAR Narrative Generator

## Version: 1.1 (Final as-built)

## Goal
Combine matched molecular pair (MMP) analysis (Phase 08) with Claude structured generation (Phase 63). Identify MMPs from the CETP inhibitor library using Tanimoto similarity-based pairing, compute activity cliffs (delta pIC50), and ask Claude to generate SAR narratives for the top-5 pairs.

## CLI
```bash
PYTHONUTF8=1 python main.py --top 5 --model claude-haiku-4-5-20251001
```

## Outputs
- `output/mmps.csv` — all matched molecular pairs with Tanimoto similarity and delta pIC50
- `output/sar_narratives.json` — Claude-generated SAR narratives for top-5 pairs
- Console summary

## Logic
1. Load `compounds.csv`
2. Compute Morgan fingerprints (radius=2, nBits=2048) for all valid compounds
3. Compute pairwise Tanimoto similarities
4. Filter pairs with Tanimoto >= 0.5 (structurally related) and |delta pIC50| > 0
5. Rank by |delta pIC50| descending — these are the most informative activity cliffs
6. Select top-N pairs
7. Send to Claude in single prompt: for each pair, generate a SAR narrative explaining what structural change caused the activity difference
8. Save outputs

## Key Concepts
- Matched molecular pairs as structure-activity relationship tool
- Tanimoto similarity for structural relatedness
- Claude generating domain-specific scientific narratives
- Activity cliffs: high structural similarity + large activity difference

## Verification Checklist
- [x] `--help` works
- [x] MMP computation runs without error
- [x] Top-5 pairs are correctly ranked by |delta pIC50|
- [x] Claude API call succeeds
- [x] Output JSON and CSV are well-formed

## Results
- 116 matched molecular pairs found (Tanimoto >= 0.5) from 45 compounds
- Top activity cliff: benz_004_CF3 vs benz_010_NMe2 (delta pIC50 = 2.050, Tanimoto = 0.543)
- Claude generated 5 SAR narratives via single API call (claude-haiku-4-5-20251001)
- All outputs well-formed (mmps.csv, sar_narratives.json)

## Key Findings
- Benzamide scaffold dominates the top activity cliffs (all top-5 are benz_ pairs)
- CF3 and CN substituents appear in the most potent compounds, NMe2 and OH in the least
- Lean API: one Claude call for all 5 narrative blocks

## Deviations
- Fixed `RDLogger` import for rdkit 2025.09.6 (moved from `rdkit.Chem` to `rdkit`)

## Risks
- Tanimoto-based pairing is approximate (not true MMP via SMILES fragmentation) — acceptable for demonstration
- Small dataset limits number of meaningful pairs
