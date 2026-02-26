# 02_choose_batch_column.py
# Chooses best batch column for ComBat from meta_final.csv

import pandas as pd
from pathlib import Path

IN_META = Path("meta_final.csv")
OUT_TABLE = Path("batch_diagnostics.csv")
OUT_PICK  = Path("chosen_batch_column.txt")

meta = pd.read_csv(IN_META)

candidates = [c for c in ["SMGEBTCH", "SMNABTCH"] if c in meta.columns]
if not candidates:
    raise ValueError("Expected SMGEBTCH and/or SMNABTCH in meta_final.csv, but none found.")

def batch_stats(col: str):
    s = meta[col].astype("string")
    valid = s.dropna()
    counts = valid.value_counts()

    return {
        "batch_col": col,
        "n_samples_total": int(meta.shape[0]),
        "n_samples_nonmissing": int(valid.shape[0]),
        "n_batches": int(counts.shape[0]),
        "min_batch_size": int(counts.min()) if len(counts) else 0,
        "median_batch_size": float(counts.quantile(0.50)) if len(counts) else 0.0,
        "max_batch_size": int(counts.max()) if len(counts) else 0,
        "n_batches_lt5": int((counts < 5).sum()),
        "n_batches_lt10": int((counts < 10).sum()),
    }

diag = pd.DataFrame([batch_stats(c) for c in candidates])

# Prefer fewer tiny batches; if tie, prefer larger min/median batch size
diag = diag.sort_values(
    by=["n_batches_lt5", "n_batches_lt10", "min_batch_size", "median_batch_size"],
    ascending=[True, True, False, False],
)

chosen = diag.iloc[0]["batch_col"]

diag.to_csv(OUT_TABLE, index=False)
OUT_PICK.write_text(str(chosen), encoding="utf-8")

print("Saved:", OUT_TABLE.resolve())
print("Chosen batch column:", chosen)
print("Saved:", OUT_PICK.resolve())
print("\nDiagnostics:")
print(diag.to_string(index=False))
