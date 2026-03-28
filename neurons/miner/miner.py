import os
os.environ['CUBLAS_WORKSPACE_CONFIG'] = ':4096:8'

import sys
import json
import traceback
import time
import gc
import psutil

import bittensor as bt
import pandas as pd
from rdkit import Chem
from pathlib import Path

BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__)))
PARENT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(PARENT_DIR)

from .random_sampler import run_sampler
from combinatorial_db.reactions import get_smiles_from_reaction

# Lazy import scoring to avoid memory spike during startup
# Will be imported only if needed (scoring enabled in config)
_scoring_loaded = False
_psichic_wrapper = None

def _load_scoring_module():
    """Lazy load scoring module only when needed"""
    global _scoring_loaded, _psichic_wrapper
    if _scoring_loaded:
        return
    
    try:
        from validator.scoring import score_molecules_json
        from PSICHIC.wrapper import PsichicWrapper
        _psichic_wrapper = PsichicWrapper()
        _scoring_loaded = True
        bt.logging.info("Scoring module loaded successfully")
    except Exception as e:
        bt.logging.error(f"Failed to load scoring module: {e}")
        _scoring_loaded = True  # Don't retry

#DB_PATH = str(Path(nova_ph2.__file__).resolve().parent / "combinatorial_db" / "molecules.sqlite")
DB_PATH = str(Path(PARENT_DIR).resolve().parent / "combinatorial_db" / "molecules.sqlite")

def get_config(input_file: os.path = os.path.join(PARENT_DIR, "..", "input.json")):
    """
    Get config from input file
    """
    with open(input_file, "r") as f:
        d = json.load(f)
    config = {**d.get("config", {}), **d.get("challenge", {})}
    return config

def iterative_sampling_loop(
    db_path: str,
    sampler_file_path: str,
    output_path: str,
    config: dict,
    save_all_scores: bool = False
) -> None:
    """
    Infinite loop, runs until orchestrator kills it:
      1) Sample n molecules
      2) Score them
      3) Merge with previous top x, deduplicate, sort, select top x
      4) Write top x to file (overwrite) each iteration
    """
    n_samples = config["num_molecules"] * 5
    MAX_MEMORY_MB = 2048
    GC_INTERVAL = 50  # Run GC every 50 iterations (optimized)

    # Lazy-load PSICHIC only if scoring is enabled in config
    if config.get("enable_scoring", False):
        bt.logging.info("[Miner] Initializing PSICHIC model (once for full cycle)...")
        _load_scoring_module()
        if _scoring_loaded:
            import validator.scoring as scoring_module
            scoring_module.psichic = _psichic_wrapper
            bt.logging.info("[Miner] PSICHIC model ready.")
    else:
        bt.logging.info("[Miner] Scoring disabled (testnet mode)")

    top_pool = pd.DataFrame(columns=["name", "smiles", "InChIKey", "score"])

    iteration = 0
    while True:
        iteration += 1
        
        # Memory check & GC every N iterations
        if iteration % GC_INTERVAL == 0:
            gc.collect()
            try:
                mem_mb = psutil.Process(os.getpid()).memory_info().rss / (1024 * 1024)
                bt.logging.info(f"[Miner] Memory: {mem_mb:.1f}MB")
                if mem_mb > MAX_MEMORY_MB:
                    bt.logging.warning(f"[Miner] Memory exceeded {MAX_MEMORY_MB}MB — forcing GC")
                    gc.collect()
            except:
                pass
        
        bt.logging.info(f"[Miner] Iteration {iteration}: sampling {n_samples} molecules")

        sampler_data = run_sampler(n_samples=n_samples, 
                        subnet_config=config, 
                        output_path=sampler_file_path,
                        save_to_file=True,
                        db_path=db_path,
                        )
        
        if not sampler_data:
            bt.logging.warning("[Miner] No valid molecules produced; continuing")
            continue

        # Skip scoring if module not loaded (testnet mode)
        if not _scoring_loaded:
            bt.logging.warning("[Miner] Scoring module not available; using default scores")
            score_dict = {name: 0.5 for name in sampler_data.get("molecules", [])}
        else:
            from validator.scoring import score_molecules_json
            score_dict = score_molecules_json(sampler_file_path, 
                                             list(config["target_sequences"].keys()), 
                                             list(config["antitarget_sequences"].keys()), 
                                             config)
        
        if not score_dict:
            bt.logging.warning("[Miner] Scoring failed or mismatched; continuing")
            continue

        # Calculate final scores per molecule
        batch_scores = calculate_final_scores(score_dict, sampler_data, config, save_all_scores, iteration)
        
        # Explicit cleanup: clear stale references after scoring
        del score_dict
        del sampler_data
        gc.collect()

        # Merge, deduplicate, sort and take top x
        top_pool = pd.concat([top_pool, batch_scores])
        top_pool = top_pool.drop_duplicates(subset=["InChIKey"], keep="first")
        top_pool = top_pool.sort_values(by="score", ascending=False)
        top_pool = top_pool.head(config["num_molecules"])

        # format to accepted format
        top_entries = {"molecules": top_pool["name"].tolist()}

        # write to file (atomic write to prevent corruption on timeout)
        try:
            output_dir = os.path.dirname(output_path)
            if output_dir:
                os.makedirs(output_dir, exist_ok=True)
            
            # Write to temporary file first, then rename (atomic on POSIX)
            import tempfile
            temp_fd, temp_path = tempfile.mkstemp(dir=output_dir or None)
            try:
                with os.fdopen(temp_fd, 'w') as f:
                    json.dump(top_entries, f, ensure_ascii=False, indent=2)
                # Atomic rename
                import shutil
                shutil.move(temp_path, output_path)
            except:
                os.close(temp_fd)
                os.unlink(temp_path)
                raise
            
            bt.logging.info(f"[Miner] Wrote {config['num_molecules']} top molecules to {output_path}")
        except Exception as e:
            bt.logging.error(f"[Miner] Error writing output to {output_path}: {e}")
            raise
        bt.logging.info(f"[Miner] Average score: {top_pool['score'].mean()}")
        
        # Explicit cleanup: clear batch data after writing
        del batch_scores
        gc.collect()

def calculate_final_scores(score_dict: dict, 
        sampler_data: dict, 
        config: dict, 
        save_all_scores: bool = True,
        current_epoch: int = 0) -> pd.DataFrame:
    """
    Calculate final scores per molecule
    """

    names = sampler_data["molecules"]
    smiles = [get_smiles_from_reaction(name) for name in names]

    # Calculate InChIKey for each molecule to deduplicate molecules after merging
    inchikey_list = []
    
    for s in smiles:
        try:
            inchikey_list.append(Chem.MolToInchiKey(Chem.MolFromSmiles(s)))
        except Exception as e:
            bt.logging.error(f"Error calculating InChIKey for {s}: {e}")
            inchikey_list.append(None)

    # Calculate final scores for each molecule
    targets = score_dict[0]['ps_target_scores']
    antitargets = score_dict[0]['ps_antitarget_scores']
    final_scores = []
    for mol_idx in range(len(names)):
        # target average
        target_scores_for_mol = [target_list[mol_idx] for target_list in targets]
        avg_target = sum(target_scores_for_mol) / len(target_scores_for_mol)

        # antitarget average
        antitarget_scores_for_mol = [antitarget_list[mol_idx] for antitarget_list in antitargets]
        avg_antitarget = sum(antitarget_scores_for_mol) / len(antitarget_scores_for_mol)

        # final score
        score = avg_target - (config["antitarget_weight"] * avg_antitarget)
        final_scores.append(score)

    # Store final scores in dataframe
    batch_scores = pd.DataFrame({
        "name": names,
        "smiles": smiles,
        "InChIKey": inchikey_list,
        "score": final_scores
    })

    # NOTE: save_all_scores disabled — /workspace is read-only in Docker sandbox
    # Writing debug files there crashes the miner. Scores are tracked in top_pool instead.

    return batch_scores

def main(config: dict):
    # Use /tmp for intermediate files — /workspace is read-only in Docker sandbox
    iterative_sampling_loop(
        db_path=DB_PATH,
        sampler_file_path="/tmp/sampler_file.json",
        output_path=os.environ.get('NOVA_OUTPUT_PATH', '/output/result.json'),
        config=config,
        save_all_scores=False,
    )
 

if __name__ == "__main__":
    config = get_config()
    main(config)
