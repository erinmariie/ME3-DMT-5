import pandas as pd
import ast
from pathlib import Path


def import_xl():
    # === CONFIGURATION ===
    OUTPUT_DIR = Path("y2025_2026") # Existing output folder
    EXCEL_FILE = OUTPUT_DIR / "variables.xlsx"       # Input Excel file
    OUTPUT_FILE = OUTPUT_DIR / "variables.py"  # Output file path

    # === STEP 1: Load Excel ===
    # Expecting first column = variable name, second column = value
    df = pd.read_excel(EXCEL_FILE, header=None, names=["variable", "value"])

    # === STEP 2: Try to infer types ===
    def infer_type(value):
        """
        Try to convert Excel string values into Python literals (int, float, bool, etc.)
        Falls back to string if literal_eval fails.
        """
        if pd.isna(value):
            return None
        if isinstance(value, (int, float, bool)):
            return value
        if isinstance(value, str):
            value = value.strip()
            try:
                return ast.literal_eval(value)
            except (ValueError, SyntaxError):
                return value  # Leave as string
        return value

    df["value"] = df["value"].apply(infer_type)

    # === STEP 3: Generate 2025-26/variables.py ===
    # Ensure output folder exists
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    with open(OUTPUT_FILE, "w", encoding="utf-8") as f:
        f.write("# This file was auto-generated from variables.xlsx\n")
        f.write("# Do not edit manually — edit the Excel instead.\n\n")

        for _, row in df.iterrows():
            name = str(row["variable"]).strip()
            value = row["value"]
            if name and not pd.isna(name):
                f.write(f"{name} = {repr(value)}\n")

    print(f"Generated '{OUTPUT_FILE}' with {len(df)} variables.")
