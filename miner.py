#!/usr/bin/env python3
"""
Wrapper script for NOVA miner — runs neurons/miner/miner.py
Designed to be called as: python miner.py
Reads input from /workspace/input.json
Writes output to /output/result.json
"""
import sys
import os
import json
import resource
import psutil
import gc
import logging

# Disable Bittensor's verbose multiprocessing logging queue
# This prevents memory exhaustion from log queue buildup
logging.getLogger("bittensor").setLevel(logging.WARNING)
os.environ['BT_LOGGING_HANDLER'] = 'off'

# Add neurons dir to path so we can import modules
BASE_DIR = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(BASE_DIR, "neurons", "miner"))  # For random_sampler, etc
sys.path.insert(0, os.path.join(BASE_DIR, "neurons"))  # For validator, miner, etc
sys.path.insert(0, BASE_DIR)  # For combinatorial_db, PSICHIC, etc

# Override the miner output path before importing
os.environ['NOVA_OUTPUT_PATH'] = '/output/result.json'
os.environ['NOVA_INPUT_PATH'] = '/workspace/input.json'

# Memory controls
MAX_MEMORY_MB = 2048  # 2GB soft limit
MAX_MEMORY_HARD_MB = 3000  # 3GB hard limit
GC_INTERVAL = 50  # Run GC every 50 iterations (reduced from 10 for 10-15% CPU savings)

def set_memory_limits():
    """Set process memory limits"""
    try:
        # Soft limit (can exceed with warning)
        soft = MAX_MEMORY_MB * 1024 * 1024
        hard = MAX_MEMORY_HARD_MB * 1024 * 1024
        resource.setrlimit(resource.RLIMIT_AS, (soft, hard))
        print(f"Memory limits set: soft={MAX_MEMORY_MB}MB, hard={MAX_MEMORY_HARD_MB}MB")
    except Exception as e:
        print(f"Warning: Could not set memory limits: {e}")

def check_memory():
    """Monitor memory usage and trigger GC if needed"""
    try:
        process = psutil.Process(os.getpid())
        mem_mb = process.memory_info().rss / (1024 * 1024)
        
        if mem_mb > MAX_MEMORY_MB:
            print(f"⚠️  Memory warning: {mem_mb:.1f}MB > {MAX_MEMORY_MB}MB — forcing GC")
            gc.collect()
            mem_mb = process.memory_info().rss / (1024 * 1024)
            print(f"   After GC: {mem_mb:.1f}MB")
        
        return mem_mb
    except Exception as e:
        print(f"Warning: Could not check memory: {e}")
        return 0

def get_config(input_file: str = "/workspace/input.json"):
    """
    Get config from input file (reads from /workspace/input.json in Docker sandbox)
    Falls back to config/config.yaml if /workspace/input.json doesn't exist
    """
    if os.path.exists(input_file):
        with open(input_file, "r") as f:
            d = json.load(f)
        config = {**d.get("config", {}), **d.get("challenge", {})}
        return config
    else:
        # For local testing: load from config/config.yaml
        import yaml
        with open(os.path.join(BASE_DIR, "config/config.yaml"), "r") as f:
            yaml_config = yaml.safe_load(f)
        # Create a minimal test config
        config = {
            "num_molecules": yaml_config["molecule_validation"]["num_molecules"],
            "antitarget_weight": yaml_config["molecule_validation"]["antitarget_weight"],
            "target_sequences": {"1ABC": "MKTII..."},  # placeholder
            "antitarget_sequences": {"2XYZ": "MSTIV..."},  # placeholder
            "allowed_reaction": "rxn:3:0:0",  # Default to reaction 3
            **yaml_config
        }
        return config

# Now run the actual miner
from neurons.miner.miner import iterative_sampling_loop, get_config as orig_get_config
from combinatorial_db.reactions import get_smiles_from_reaction

DB_PATH = os.path.join(BASE_DIR, "combinatorial_db", "molecules.sqlite")

if __name__ == "__main__":
    try:
        set_memory_limits()
        gc.enable()
        gc.set_threshold(700, 10, 10)  # Aggressive GC thresholds
        
        config = get_config()
        output_path = os.environ.get('NOVA_OUTPUT_PATH', '/output/result.json')
        
        # Ensure output directory exists
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        
        print(f"Starting NOVA miner...")
        print(f"  Config: {config}")
        print(f"  DB: {DB_PATH}")
        print(f"  Output: {output_path}")
        
        iterative_sampling_loop(
            db_path=DB_PATH,
            sampler_file_path="/tmp/sampler_file.json",  # /workspace is read-only in sandbox
            output_path=output_path,
            config=config,
            save_all_scores=False,  # disabled — /workspace is read-only, no debug files
        )
    except KeyboardInterrupt:
        print("Miner interrupted")
        sys.exit(0)
    except Exception as e:
        print(f"Miner error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
