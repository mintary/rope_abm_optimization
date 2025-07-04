import platform
from pathlib import Path
from typing import Optional
import csv
import numpy as np
import subprocess
import logging
import sys
import os

def prepare_sample(
    values: list[float],
    sample_output_path: Path,
    logger: logging.Logger = logging.getLogger(__name__),
) -> None:
    """
    Prepare a sample input file for the ABM simulation.
    This contains the sampled parameters in a single row.
    Ensures the parent directory exists before writing.
    """
    logger.info(f"Preparing sample with values: {values}")
    sample_output_path = sample_output_path.resolve()
    sample_output_path.parent.mkdir(parents=True, exist_ok=True)
    arr = np.array(values, dtype=np.float32)
    arr = arr.reshape(1, -1) # Written as a single row
    np.savetxt(sample_output_path, arr, delimiter=",", fmt="%.6f")

# Mapping from CSV column names to field names
CSV_TO_FIELD = {
    "clock": "tick",
    "TNF": "TNF",
    "TGF": "TGF",
    "FGF": "FGF",
    "IL6": "IL6",
    "IL8": "IL8",
    "IL10": "IL10",
    "Tropocollagen": "Tropocollagen",
    "Collagen": "Collagen",
    "FragentedCollagen": "FragentedCollagen",
    "Tropoelastin": "Tropoelastin",
    "Elastin": "Elastin",
    "FragmentedElastin": "FragmentedElastin",
    "HA": "HA",
    "FragmentedHA": "FragmentedHA",
    "Damage": "Damage",
    "ActivatedFibroblast": "ActivatedFibroblast",
    "Fibroblast": "Fibroblast",
    "Elastic Mod (Pa)": "ElasticMod",
    "Swelling Ratio": "SwellingRatio",
    "Mass Loss (%)": "MassLoss",
}

def run_abm_simulation(
    bin_path: Path,
    bin_stderr_path: Path,
    bin_stdout_path: Path,
    biomarker_output_dir: Path,
    sample_path: Path,
    config_path: Path,
    cwd: Optional[Path] = None,
    num_ticks: int = 289,
    tracked_ticks: list[int] = [144, 288],
    output_biomarkers: Optional[list[str]] = None,
) -> list[dict]:
    """
    Run the ABM simulation (which is provided as a binary executable).
    Returns a list of dicts with keys: 'tick', 'Fibroblast', 'ActivatedFibroblast', 'Collagen'.
    """
    # Initialize logger
    logging.basicConfig(level=logging.INFO)
    logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s')
    logger = logging.getLogger(__name__)
    
    # Ensure the sample was created and exists
    if not sample_path.exists():
        logger.error(f"Sample file not found: {sample_path}")
        raise FileNotFoundError(f"Sample file not found: {sample_path}")
    
    # Determine the executable name based on the platform
    if platform.system() == "Windows":
        exe_name = "testRun.exe"
    else:
        exe_name = "testRun"
    test_run_path = bin_path / exe_name

    logger.info(f"Full binary path: {test_run_path}")
    
    if not test_run_path.exists():
        logger.error(f"ABM binary not found: {test_run_path}")
        raise FileNotFoundError(f"ABM binary not found: {test_run_path}")
    
    logger.info(f"Running ABM simulation with binary: {test_run_path}")
    results: list[dict]= []
    
    # Ensure output directories exist
    biomarker_output_dir.mkdir(parents=True, exist_ok=True)
    bin_stderr_path.parent.mkdir(parents=True, exist_ok=True)
    bin_stdout_path.parent.mkdir(parents=True, exist_ok=True)

    # Build the executable path for the subprocess command
    if cwd:
        exe_path = str((Path(cwd).resolve() / "bin" / exe_name))
    else:
        exe_path = os.fspath(test_run_path)
    command = [
        exe_path,
        "--numticks", str(num_ticks),
        "--inputfile", str(config_path),
        "--wxw", "0.6",
        "--wyw", "0.6",
        "--wzw", "0.6"
    ]
    subprocess.run(
        command,
        stdout=open(bin_stdout_path, "w"),
        stderr=open(bin_stderr_path, "w"),
        cwd=str(Path(cwd).resolve()) if cwd else None
    )
    # Check if the output biomarkers file exists
    output_biomarkers_path = biomarker_output_dir / "Output_Biomarkers.csv"
    logger.debug(f"Looking for output biomarkers file at: {output_biomarkers_path}")
    
    if not output_biomarkers_path.exists():
        logger.error(f"Output biomarkers file not found: {output_biomarkers_path}")
        logger.error(f"Directory contents: {list(biomarker_output_dir.glob('*'))}")
        raise FileNotFoundError(f"Output biomarkers file not found: {output_biomarkers_path}")
        
    with output_biomarkers_path.open("r") as f:
        reader = csv.reader(f)
        header = next(reader)
        logger.debug(f"Header from output biomarkers: {header}")
        fieldnames = [CSV_TO_FIELD.get(h.strip(), h.strip()) for h in header]
        results = []
        for row in reader:
            row_dict = {field: value for field, value in zip(fieldnames, row)}
            tick = int(row_dict["tick"])
            if tick in tracked_ticks:
                try:
                    fib = float(row_dict.get("Fibroblast", 0.0))
                    act_fib = float(row_dict.get("ActivatedFibroblast", 0.0))
                    collagen = float(row_dict.get("Collagen", 0.0))
                    result = {"tick": tick, "Fibroblast": fib, "ActivatedFibroblast": act_fib, "Collagen": collagen}
                    results.append(result)
                except (ValueError, TypeError):
                    logger.warning(f"Could not convert values for tick {tick} to float.")
        for result in results:
            logger.info(f"Tick {result['tick']}: Fibroblast: {result['Fibroblast']:.4f}, ActivatedFibroblast: {result['ActivatedFibroblast']:.4f}, Collagen: {result['Collagen']:.4f}")
    return results
