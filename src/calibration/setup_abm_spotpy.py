from spotpy.parameter import Uniform
from spotpy.objectivefunctions import rmse
import os
import spotpy
import uuid
import numpy as np
import shutil
import logging
from pathlib import Path
from typing import Optional, Callable
from src.calibration.simulation import prepare_sample, run_abm_simulation
import csv

CONFIG_FILE_NAMES = [
    "config_Scaffold_GH2.txt",
    "config_Scaffold_GH5.txt",
    "config_Scaffold_GH10.txt"
]

TICKS_PER_DAY = 48

class setup_abm_spotpy(object):
    def __init__(self, 
                 run_dir_parent_path: Path,
                 bin_dir: Path,
                 config_file_dir: Path,
                 params: list[spotpy.parameter.Uniform],
                 num_ticks,
                 tracked_ticks: list[int],
                 experimental_data_csv: Path,
                 error_function: Callable[[list[float], list[float]], float] = lambda sim, obs: float(np.mean((np.array(sim) - np.array(obs)) ** 2)) ** 0.5,
                 logger: logging.Logger = logging.getLogger(__name__)
            ):
        """
        Initialize the RopeSetup.
        We use the run_dir_parent_path to store subdirectories for each run.        
        """
        self.params = params
        self.run_dir_parent_path = run_dir_parent_path
        self.bin_dir = bin_dir
        self.config_file_dir = config_file_dir
        self.num_ticks = num_ticks
        self.tracked_ticks = tracked_ticks
        self.error_function = error_function
        self.experimental_data_csv = experimental_data_csv
        self.logger = logger

    def parameters(self):
        return spotpy.parameter.generate(self.params)

    def _extract_total_fibroblast_and_collagen(self, csv_path: Path, tracked_ticks: list[int]) -> list[float]:
        """
        Reads the CSV and returns a flat list: [total_fib_day3, total_collagen_day3, total_fib_day6, total_collagen_day6, ...]
        where total_fib = Fibroblast + ActivatedFibroblast.
        Only includes rows for tracked_ticks.
        """
        results = []
        with csv_path.open(newline='') as csvfile:
            reader = csv.DictReader(csvfile)
            reader.fieldnames = [h.strip() for h in reader.fieldnames] if reader.fieldnames else None
            for row in reader:
                clean_row = {k.strip(): v for k, v in row.items()}
                tick = int(clean_row['Day']) if 'Day' in clean_row else int(clean_row['Tick'])
                if 'Day' in clean_row:
                    tick = int(float(clean_row['Day']) * TICKS_PER_DAY)
                if tick in tracked_ticks:
                    try:
                        fib = float(clean_row.get('Fibroblast', 0.0))
                        collagen = float(clean_row.get('Collagen', 0.0))
                        total_fib = fib
                        results.extend([total_fib, collagen])
                    except (ValueError, TypeError):
                        self.logger.warning(f"Could not convert values for tick {tick} to float.")
        return results

    def simulation(self, 
                   vector, 
                   run_id: Optional[str] = None, 
                   config_path_names: list[str] = CONFIG_FILE_NAMES
                   ) -> list[float]:
        """
        Run the ABM simulation with the provided parameter vector.
        This prepares the sample input file, runs the simulation, and collects the results.
        We run the ABM simulation with different config files to ensure stability and average the results.

        Because the simulation can be run in parallel, we can use run_id to identify the specific run.
        This creates a new directory for each run_id, where the results will be stored.
        If run_id is None, a new UUID will be generated.
        """
        # Ensure logger is configured in each process (for multiprocessing)
        if not logging.getLogger().hasHandlers():
            logging.basicConfig(
                level=self.logger.level,
                format='%(asctime)s %(levelname)s %(name)s: %(message)s',
                handlers=[logging.StreamHandler()]
            )

        # Update the logger to include the run_id in the log messages, preserving original settings
        for handler in self.logger.handlers:
            fmt = handler.formatter._fmt if handler.formatter else "%(levelname)s:%(name)s:%(message)s"
            handler.setFormatter(logging.Formatter(f"[{run_id}] {fmt}"))

        self.logger.debug(f"Running simulation with vector: {vector} and run_id: {run_id}")
        if run_id is None:
            run_id = uuid.uuid4().hex # Or use index from the sampler
        
        # We create a directory for the run_id under the parent path
        # This allows us to keep the results organized and separate for each run
        run_id_dir = self.run_dir_parent_path / f"run_{run_id}"
        run_id_dir.mkdir(parents=True, exist_ok=True)

        self.logger.debug(f"Created directory {run_id_dir}")

        # Copy the binary executable to the run directory
        run_id_dir_bin = run_id_dir / "bin"
        run_id_dir_bin.mkdir(parents=True, exist_ok=True)
        
        # Log the bin directory contents before copying
        self.logger.debug(f"Source bin directory {self.bin_dir} contains: {list(self.bin_dir.glob('*'))}")
        
        # Copy the binary files to the run directory
        for item in self.bin_dir.glob('*'):
            if item.is_file():
                dest = run_id_dir_bin / item.name
                shutil.copy2(item, dest)
                self.logger.debug(f"Copied {item} to {dest}")
        
        # Log the bin directory contents after copying
        self.logger.debug(f"Run bin directory {run_id_dir_bin} now contains: {list(run_id_dir_bin.glob('*'))}")

        # Create configFiles directory in the run directory
        config_dir = run_id_dir / "configFiles"
        config_dir.mkdir(parents=True, exist_ok=True)
        
        # Copy the config files to the run directory's configFiles subdirectory
        for config_file in config_path_names:
            source = self.config_file_dir / config_file
            dest = config_dir / config_file
            shutil.copy2(source, dest)
            self.logger.debug(f"Copied config file from {source} to {dest}")

        sample_output_path = run_id_dir / "Sample.txt"
        prepare_sample(vector, sample_output_path, self.logger)

        # Prepare aggregate results
        aggregate_results: list[list[dict]] = []

        runs = len(config_path_names)

        for i in range(runs):
            self.logger.debug(f"Running ABM simulation for iteration: {i+1}")
            biomarker_output_dir = run_id_dir / "output" # The bin executable will write to this directory
            biomarker_output_dir.mkdir(parents=True, exist_ok=True)  # Ensure output dir exists
            bin_stderr_path = biomarker_output_dir / "ABM_simulation_stderr.txt"
            bin_stdout_path = biomarker_output_dir / "ABM_simulation_stdout.txt"

            self.logger.debug(f"Current working directory set to {run_id_dir}")
            run_results = run_abm_simulation(
                bin_path=run_id_dir_bin,
                bin_stderr_path=bin_stderr_path,
                bin_stdout_path=bin_stdout_path,
                biomarker_output_dir=biomarker_output_dir,
                sample_path=sample_output_path,
                config_path=config_dir / config_path_names[i],
                num_ticks=self.num_ticks,
                tracked_ticks=self.tracked_ticks,
                cwd=run_id_dir
            )
            self.logger.debug(f"ABM simulation completed for run_id: {run_id}, iteration: {i+1}")
            aggregate_results.append(run_results)
        self.logger.debug(f"All ABM simulations completed for run_id: {run_id}")

        print(aggregate_results)
        # Average the results across the runs for Fibroblast, ActivatedFibroblast, Collagen
        averaged_rows = []
        for i in range(len(aggregate_results[0])):
            tick = aggregate_results[0][i]["tick"]
            fibs = [run[i]["Fibroblast"] for run in aggregate_results]
            act_fibs = [run[i]["ActivatedFibroblast"] for run in aggregate_results]
            collagens = [run[i]["Collagen"] for run in aggregate_results]
            avg_fib = float(np.mean(fibs))
            avg_act_fib = float(np.mean(act_fibs))
            avg_collagen = float(np.mean(collagens))
            averaged_rows.append({'Tick': tick, 'Fibroblast': avg_fib, 'ActivatedFibroblast': avg_act_fib, 'Collagen': avg_collagen})
        # Write averaged results to a temp CSV and extract as floats
        import tempfile, csv as pycsv
        with tempfile.NamedTemporaryFile(mode='w+', newline='', delete=False) as tmpfile:
            writer = pycsv.DictWriter(tmpfile, fieldnames=['Tick', 'Fibroblast', 'ActivatedFibroblast', 'Collagen'])
            writer.writeheader()
            for row in averaged_rows:
                writer.writerow(row)
            tmpfile.flush()
            return self._extract_total_fibroblast_and_collagen(Path(tmpfile.name), self.tracked_ticks)

    def evaluation(self) -> list[float]:
        """
        Read the observed data from the provided path and return a flat list:
        [total_fib_day3, total_collagen_day3, total_fib_day6, total_collagen_day6, ...]
        """
        experimental_data_csv = self.experimental_data_csv
        if not experimental_data_csv.exists():
            self.logger.error(f"Observed data file not found: {experimental_data_csv}")
            raise FileNotFoundError(f"Observed data file not found: {experimental_data_csv}")
        return self._extract_total_fibroblast_and_collagen(experimental_data_csv, self.tracked_ticks)

    def objectivefunction(self, simulation: list[float], evaluation: list[float]) -> float:
        if len(simulation) != len(evaluation):
            raise ValueError(f"Simulation and evaluation lists must have the same length. Got {len(simulation)} and {len(evaluation)}.")
        result = self.error_function(simulation, evaluation)
        # Defensive: reduce to scalar if needed
        if isinstance(result, (np.ndarray, list)):
            result = np.asarray(result).flatten()
            if result.size == 1:
                return float(result.item())
            else:
                return float(result.mean())
        return float(result)