"""Phase 97 — MMP Analysis + Claude SAR Narrative Generator.

Identifies matched molecular pairs from the CETP inhibitor library,
ranks by activity cliff magnitude, and generates Claude SAR narratives.
"""
import sys
import os

if sys.platform == "win32":
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8", errors="replace")

import argparse
import json
from pathlib import Path
from itertools import combinations

import pandas as pd
import numpy as np
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem
from anthropic import Anthropic
from dotenv import load_dotenv

RDLogger.DisableLog("rdApp.*")

env_path = Path(__file__).resolve().parent / ".env"
if env_path.exists():
    load_dotenv(env_path)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="MMP analysis + Claude SAR narrative generator",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--data", default="compounds.csv", help="Input CSV")
    p.add_argument("--top", type=int, default=5, help="Number of top MMPs for SAR narrative")
    p.add_argument("--sim-threshold", type=float, default=0.5, help="Min Tanimoto for MMP")
    p.add_argument("--model", default="claude-haiku-4-5-20251001", help="Claude model")
    p.add_argument("--output", default="output", help="Output directory")
    return p.parse_args()


def compute_mmps(df: pd.DataFrame, sim_threshold: float) -> pd.DataFrame:
    """Compute Tanimoto-based matched molecular pairs."""
    # Generate fingerprints
    fps = {}
    for _, row in df.iterrows():
        mol = Chem.MolFromSmiles(row["smiles"])
        if mol is not None:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
            fps[row["compound_name"]] = (fp, row["smiles"], row["pic50"])

    names = list(fps.keys())
    pairs = []
    for i, j in combinations(range(len(names)), 2):
        n1, n2 = names[i], names[j]
        fp1, smi1, pic1 = fps[n1]
        fp2, smi2, pic2 = fps[n2]
        sim = DataStructs.TanimotoSimilarity(fp1, fp2)
        if sim >= sim_threshold:
            delta = abs(pic1 - pic2)
            pairs.append({
                "compound_a": n1, "smiles_a": smi1, "pic50_a": pic1,
                "compound_b": n2, "smiles_b": smi2, "pic50_b": pic2,
                "tanimoto": round(sim, 3),
                "delta_pic50": round(delta, 3),
            })

    mmp_df = pd.DataFrame(pairs)
    if len(mmp_df) > 0:
        mmp_df = mmp_df.sort_values("delta_pic50", ascending=False).reset_index(drop=True)
    return mmp_df


def generate_sar_narratives(top_pairs: pd.DataFrame, model: str) -> list[dict]:
    """Ask Claude to generate SAR narratives for top MMP pairs."""
    client = Anthropic()

    pair_block = "\n".join(
        f"Pair {i+1}: {row['compound_a']} ({row['smiles_a']}, pIC50={row['pic50_a']:.2f}) vs "
        f"{row['compound_b']} ({row['smiles_b']}, pIC50={row['pic50_b']:.2f}) | "
        f"Tanimoto={row['tanimoto']:.3f}, delta_pIC50={row['delta_pic50']:.3f}"
        for i, row in top_pairs.iterrows()
    )

    prompt = f"""You are a medicinal chemist analyzing structure-activity relationships (SAR).

Below are matched molecular pairs from a CETP inhibitor library. Each pair consists of
structurally similar compounds (high Tanimoto similarity) with different activities.

{pair_block}

For each pair, generate a SAR narrative as a JSON array. Each element must have:
- "pair_index": int (1-based)
- "compound_a": string
- "compound_b": string
- "structural_change": string (describe the key structural difference)
- "activity_impact": string ("increase" or "decrease" from A to B)
- "narrative": string (2-3 sentences explaining why this structural change likely affects activity)
- "confidence": "high" | "medium" | "low"

Return ONLY the JSON array, no other text."""

    print(f"\n[2] Calling Claude ({model}) for SAR narratives...")
    response = client.messages.create(
        model=model,
        max_tokens=4096,
        messages=[{"role": "user", "content": prompt}],
    )

    text = response.content[0].text.strip()
    if text.startswith("```"):
        text = text.split("\n", 1)[1]
        if text.endswith("```"):
            text = text[: text.rfind("```")]
        text = text.strip()

    narratives = json.loads(text)
    print(f"  Generated {len(narratives)} SAR narratives")
    return narratives


def main() -> None:
    args = parse_args()
    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)

    # 1. Compute MMPs
    print(f"[1] Loading {args.data} and computing matched molecular pairs...")
    df = pd.read_csv(args.data)
    mmp_df = compute_mmps(df, args.sim_threshold)
    print(f"  Found {len(mmp_df)} pairs with Tanimoto >= {args.sim_threshold}")

    mmp_csv = out_dir / "mmps.csv"
    mmp_df.to_csv(mmp_csv, index=False)
    print(f"  Saved all MMPs to {mmp_csv}")

    if len(mmp_df) == 0:
        print("  No pairs found. Try lowering --sim-threshold.")
        return

    top_pairs = mmp_df.head(args.top).copy()
    print(f"\n  Top-{args.top} activity cliffs:")
    for _, row in top_pairs.iterrows():
        print(f"    {row['compound_a']} vs {row['compound_b']}: "
              f"delta={row['delta_pic50']:.3f}, sim={row['tanimoto']:.3f}")

    # 2. Claude SAR narratives
    narratives = generate_sar_narratives(top_pairs, args.model)

    # 3. Save
    report = {
        "phase": 97,
        "description": "MMP analysis + Claude SAR narrative generator",
        "model": args.model,
        "n_pairs": len(narratives),
        "narratives": narratives,
    }
    report_path = out_dir / "sar_narratives.json"
    with open(report_path, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2, ensure_ascii=False)
    print(f"\n[3] SAR narratives saved to {report_path}")
    print("Done.")


if __name__ == "__main__":
    main()
