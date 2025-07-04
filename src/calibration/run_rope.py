import click
import spotpy
import logging
from pathlib import Path
from src.calibration.setup_abm_spotpy import setup_abm_spotpy
from src.calibration.preprocessing import process_parameters_from_csv, extract_n_parameters
from typing import Any, Dict
from dataclasses import dataclass

@dataclass
class LoggerContext:
    logger: logging.Logger

def create_context(log_level: str = "INFO") -> LoggerContext:
    # Setup logger
    logging.basicConfig(
        level=getattr(logging, log_level.upper(), logging.INFO),
        format='%(asctime)s %(levelname)s %(name)s: %(message)s',
        handlers=[logging.StreamHandler()]
    )
    logger: logging.Logger = logging.getLogger(__name__)
    return LoggerContext(logger=logger)

@click.command()
@click.option('--log-level', '-l', default='INFO', type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'], case_sensitive=False), 
              help="Set the logging level.")
@click.option('--param-ranking', '-pr', type=click.Choice(choices=['random_forest', 'morris'], case_sensitive=False), default='random_forest',
              help="Method to rank parameters for sensitivity analysis.")
@click.option('--param-num', '-pn', default=5, type=int, 
              help="Number of parameters to rank.")
@click.option('--num-iterations', '-i', default=800, type=int, 
              help="Number of iterations for SPOTPY optimization.")
@click.option('--run-dir-parent', '-r', type=click.Path(exists=True, file_okay=False, dir_okay=True),
              help="Parent directory for simulation run subdirectories.")
@click.option('--sensitivity-analysis-csv', '-s', type=click.Path(exists=True, dir_okay=False),
              help="CSV file containing sensitivity analysis results.")
@click.option('--experimental-data-csv', '-e', required=True, type=click.Path(exists=True, dir_okay=False),
             help="CSV file containing experimental/observed data for evaluation.")
@click.option('--config-file-dir', '-c', default=Path("configFiles"), type=click.Path(exists=True, file_okay=False, dir_okay=True),
              help="Directory containing configuration files.")
@click.option('--bin-dir', '-b', default=Path("bin"), type=click.Path(exists=True, file_okay=False, dir_okay=True),
              help="Directory containing the ABM simulation binary files.")
@click.option('--parallel', '-p', type=click.Choice(['mpc', 'mpi', 'seq'], case_sensitive=False), default='mpc',
              help="Parallelization method to use.")
@click.pass_context
def run(ctx, 
        log_level: str,
        param_ranking: str, 
        param_num: int, 
        num_iterations: int, 
        run_dir_parent: Path,
        sensitivity_analysis_csv: Path,
        config_file_dir: Path,
        bin_dir: Path,
        experimental_data_csv: Path,
        parallel: str
        ):
    """
    Command line interface for running the ABM simulation with spotpy.
    """
    # Ensure all paths are Path objects (Click may pass them as str)
    run_dir_parent = Path(run_dir_parent)
    sensitivity_analysis_csv = Path(sensitivity_analysis_csv)
    config_file_dir = Path(config_file_dir)
    bin_dir = Path(bin_dir)
    experimental_data_csv = Path(experimental_data_csv)

    logger_ctx = create_context(log_level)
    ctx.obj = logger_ctx
    logger = ctx.obj.logger

    # Prepare the parameter vector based on the ranking method
    df = process_parameters_from_csv(sensitivity_analysis_csv)
    params = extract_n_parameters(df, ranking_method=param_ranking, n=param_num)

    # Convert params into spotpy parameter objects
    spotpy_params = [
        spotpy.parameter.Uniform(
            name=row['parameter_number'],
            low=row['lower_bound'],
            high=row['upper_bound'],
            optguess=row['default_value']
        ) for _, row in params.iterrows()
    ]

    spotpy_setup = setup_abm_spotpy(
        run_dir_parent_path=run_dir_parent,
        bin_dir=bin_dir,
        config_file_dir=config_file_dir,
        params=spotpy_params,
        experimental_data_csv=experimental_data_csv,
        num_ticks=289,
        tracked_ticks=[144, 288],
        logger=logger
    )

    sampler = spotpy.algorithms.rope(
        spotpy_setup,
        dbname='rope_abm_optimization',
        dbformat='csv',
        parallel=parallel,
        save_sim=True,
    )
    sampler.sample(
        repetitions= num_iterations,
    )

if __name__ == "__main__":
    run()