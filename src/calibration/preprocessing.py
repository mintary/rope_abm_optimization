import csv
import pandas as pd
from pathlib import Path

RANKING_METHODS = ["morris", "random_forest"]

def process_parameters_from_csv(file_path: Path) -> pd.DataFrame:
    """
    Reads a CSV file containing parameter names and values, and returns a DataFrame.
    """
    df = pd.read_csv(file_path)
    
    df.columns = [c.strip() for c in df.columns]
    
    cols = [
        "Parameter Number",
        "Parameter Name",
        "Morris",
        "Random Forest",
        "Lower bound",
        "Upper bound",
        "Default Value"
    ]

    # Some columns may have extra spaces or slightly different names, so use fuzzy matching if needed
    df = df[[col for col in df.columns if any(key in col for key in cols)]]
    df = df.rename(columns=lambda x: x.strip())
    
    # Rename columns for consistency
    df = df.rename(columns={
        "Parameter Number": "parameter_number",
        "Parameter Name": "parameter_name",
        "Morris": "morris_ranking",
        "Random Forest": "random_forest_ranking",
        "Lower bound": "lower_bound",
        "Upper bound": "upper_bound",
        "Default Value": "default_value"
    })
    
    # Fill missing rankings with 9999 and convert to numeric
    for col in ["morris_ranking", "random_forest_ranking"]:
        df[col] = pd.to_numeric(df[col], errors="coerce").fillna(9999).astype(int)
    
    # Return only the requested columns
    return df[["parameter_number", "parameter_name", "morris_ranking", "random_forest_ranking", "lower_bound", "upper_bound", "default_value"]]

def extract_n_parameters(df: pd.DataFrame, ranking_method: str, n: int) -> pd.DataFrame:
    """
    Extracts the top n parameters from the processed DataFrame based on the ranking method.
    """
    if ranking_method not in RANKING_METHODS:
        raise ValueError(f"Invalid ranking method: {ranking_method}. Choose from {RANKING_METHODS}.")
    ranking_col = f"{ranking_method}_ranking"
    df = df.sort_values(by=ranking_col).reset_index(drop=True)
    # Select only the top n parameters
    df = df.head(n)
    return df[["parameter_number", "parameter_name", ranking_col, "lower_bound", "upper_bound", "default_value"]]