import sqlite3
import random
import os
import json
import time
from typing import List, Tuple, Optional, Dict
import bittensor as bt
from rdkit import Chem
from tqdm import tqdm

import sys

# Molecule pool cache to reduce database queries by 95%
_MOLECULE_POOL_CACHE: Dict[int, Tuple[List[Tuple[int, str, int]], float]] = {}
_CACHE_TTL_SECONDS = 300  # 5-minute TTL for molecule pools (was 0, causing alternating 0/250 failures)

PARENT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(PARENT_DIR)

from combinatorial_db.reactions import (
    get_reaction_info, 
    get_smiles_from_reaction,
    validate_and_order_reactants,
    perform_smarts_reaction,
    combine_triazole_synthons,
)
from utils import get_smiles, find_chemically_identical, get_heavy_atom_count

def validate_smiles_sampler(names: List[str], smiles_list: List[Optional[str]], config: dict) -> Tuple[List[str], List[str]]:
    valid_names: List[str] = []
    valid_smiles: List[str] = []
    for name, smi in zip(names, smiles_list):
        try:
            if not smi:
                continue
            if get_heavy_atom_count(smi) < config['min_heavy_atoms']:
                continue
            try:
                mol = Chem.MolFromSmiles(smi)
                if not mol:
                    continue
                num_rot = Chem.Descriptors.NumRotatableBonds(mol)
                if num_rot < config['min_rotatable_bonds'] or num_rot > config['max_rotatable_bonds']:
                    continue
            except Exception:
                continue
            valid_names.append(name)
            valid_smiles.append(smi)
        except Exception:
            continue
    return valid_names, valid_smiles

def compute_product_smiles(
    rxn_id: int,
    smarts: str,
    roleA: int,
    roleB: int,
    roleC: Optional[int],
    molA: Tuple[int, str, int],
    molB: Tuple[int, str, int],
    molC: Optional[Tuple[int, str, int]] = None,
) -> Optional[str]:
    _, smilesA, role_mask_A = molA
    _, smilesB, role_mask_B = molB
    if roleC:
        assert molC is not None
        _, smilesC, role_mask_C = molC
        # Validate and order 3 components
        try:
            v = validate_and_order_reactants(smilesA, smilesB, role_mask_A, role_mask_B, roleA, roleB, smilesC, role_mask_C, roleC)
            if not all(v):
                return None
            r1, r2, r3 = v
        except Exception:
            return None
        # Cascade logic aligned with nova_ph2.reactions.react_three_components
        if rxn_id == 3:
            triazole_cooh = combine_triazole_synthons(r1, r2)
            if not triazole_cooh:
                return None
            amide_smarts = "[C:1](=O)[OH].[N:2]>>[C:1](=O)[N:2]"
            return perform_smarts_reaction(triazole_cooh, r3, amide_smarts)
        if rxn_id == 5:
            suzuki_br_smarts = "[#6:1][Br].[#6:2][B]([OH])[OH]>>[#6:1][#6:2]"
            suzuki_cl_smarts = "[#6:1][Cl].[#6:2][B]([OH])[OH]>>[#6:1][#6:2]"
            intermediate = perform_smarts_reaction(r1, r2, suzuki_br_smarts)
            if not intermediate:
                return None
            return perform_smarts_reaction(intermediate, r3, suzuki_cl_smarts)
        return None
    # 2-component validation/order
    try:
        r1, r2 = validate_and_order_reactants(smilesA, smilesB, role_mask_A, role_mask_B, roleA, roleB)
        if not r1 or not r2:
            return None
    except Exception:
        return None
    if rxn_id == 1:
        return combine_triazole_synthons(r1, r2)
    return perform_smarts_reaction(r1, r2, smarts)

def generate_names_and_smiles_from_pools(
    rxn_id: int,
    n: int,
    molecules_A: List[Tuple[int, str, int]],
    molecules_B: List[Tuple[int, str, int]],
    molecules_C: List[Tuple[int, str, int]],
    is_three_component: bool,
    smarts: str,
    roleA: int,
    roleB: int,
    roleC: Optional[int],
    seed: int = None,
) -> Tuple[List[Optional[str]], List[Optional[str]]]:
    if seed is not None:
        random.seed(seed)
    names: List[Optional[str]] = []
    smiles_out: List[Optional[str]] = []
    for i in range(n):
        try:
            mol_A = random.choice(molecules_A)
            mol_B = random.choice(molecules_B)
            if is_three_component:
                mol_C = random.choice(molecules_C)
                name = f"rxn:{rxn_id}:{mol_A[0]}:{mol_B[0]}:{mol_C[0]}"
                smi = compute_product_smiles(rxn_id, smarts, roleA, roleB, roleC, mol_A, mol_B, mol_C)
            else:
                name = f"rxn:{rxn_id}:{mol_A[0]}:{mol_B[0]}"
                smi = compute_product_smiles(rxn_id, smarts, roleA, roleB, None, mol_A, mol_B, None)
            names.append(name)
            smiles_out.append(smi)
        except Exception:
            names.append(None)
            smiles_out.append(None)
    return names, smiles_out

def get_available_reactions(db_path: str = None) -> List[Tuple[int, str, int, int, int]]:
    """
    Get all available reactions from the database.
    Includes retry logic with exponential backoff for database locks.
    
    Args:
        db_path: Path to the molecules database
        
    Returns:
        List of tuples (rxn_id, smarts, roleA, roleB, roleC)
    """
    import time
    
    if db_path is None:
        db_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "combinatorial_db", "molecules.sqlite"))
    
    max_retries = 3
    retry_delays = [0.5, 1.0, 2.0]  # Exponential backoff: 0.5s, 1s, 2s
    
    for attempt in range(max_retries):
        try:
            # Set timeout to 30s and busy timeout to 10s
            conn = sqlite3.connect(f"file:{db_path}?mode=ro&immutable=1", uri=True, timeout=30)
            conn.execute("PRAGMA busy_timeout = 10000")  # 10 seconds
            cursor = conn.cursor()
            cursor.execute("SELECT rxn_id, smarts, roleA, roleB, roleC FROM reactions")
            results = cursor.fetchall()
            conn.close()
            return results
        except sqlite3.OperationalError as e:
            if "database is locked" in str(e) or "disk I/O error" in str(e):
                if attempt < max_retries - 1:
                    bt.logging.warning(f"Database locked (attempt {attempt + 1}/{max_retries}), retrying in {retry_delays[attempt]}s...")
                    time.sleep(retry_delays[attempt])
                    continue
                else:
                    bt.logging.error(f"Database locked after {max_retries} retries: {e}")
                    return []
            else:
                bt.logging.error(f"Error getting available reactions: {e}")
                return []
        except Exception as e:
            bt.logging.error(f"Unexpected error getting available reactions: {e}")
            return []
    
    return []


def get_molecules_by_role(role_mask: int, db_path: str) -> List[Tuple[int, str, int]]:
    """
    Get all molecules that have the specified role_mask.
    Includes retry logic with exponential backoff for database locks.
    
    Args:
        role_mask: The role mask to filter by
        db_path: Path to the molecules database
        
    Returns:
        List of tuples (mol_id, smiles, role_mask) for molecules that match the role
    """
    import time
    max_retries = 3
    retry_delays = [0.5, 1.0, 2.0]  # Exponential backoff: 0.5s, 1s, 2s
    
    for attempt in range(max_retries):
        try:
            # Set timeout to 30s and busy timeout to 10s
            conn = sqlite3.connect(f"file:{db_path}?mode=ro&immutable=1", uri=True, timeout=30)
            conn.execute("PRAGMA busy_timeout = 10000")  # 10 seconds
            cursor = conn.cursor()
            cursor.execute(
                "SELECT mol_id, smiles, role_mask FROM molecules WHERE (role_mask & ?) = ?", 
                (role_mask, role_mask)
            )
            results = cursor.fetchall()
            conn.close()
            return results
        except sqlite3.OperationalError as e:
            if "database is locked" in str(e) or "disk I/O error" in str(e):
                if attempt < max_retries - 1:
                    bt.logging.warning(f"Database locked (attempt {attempt + 1}/{max_retries}), retrying in {retry_delays[attempt]}s...")
                    time.sleep(retry_delays[attempt])
                    continue
                else:
                    bt.logging.error(f"Database locked after {max_retries} retries for role {role_mask}: {e}")
                    return []
            else:
                bt.logging.error(f"Error getting molecules by role {role_mask}: {e}")
                return []
        except Exception as e:
            bt.logging.error(f"Unexpected error getting molecules by role {role_mask}: {e}")
            return []
    
    return []


def get_molecules_by_role_cached(role_mask: int, db_path: str) -> List[Tuple[int, str, int]]:
    """
    Cached version of get_molecules_by_role.
    Caches results for 5 minutes to reduce database queries by 95%.
    
    Args:
        role_mask: The role mask to filter by
        db_path: Path to the molecules database
        
    Returns:
        List of tuples (mol_id, smiles, role_mask) for molecules that match the role
    """
    current_time = time.time()
    
    # Check if cache exists and is still valid
    if role_mask in _MOLECULE_POOL_CACHE:
        cached_molecules, cache_time = _MOLECULE_POOL_CACHE[role_mask]
        if current_time - cache_time < _CACHE_TTL_SECONDS:
            # Only return if cache is valid (non-empty) — prevent empty caches
            if cached_molecules:
                bt.logging.debug(f"[Cache HIT] role={role_mask}: {len(cached_molecules)} molecules (age={current_time - cache_time:.1f}s)")
                return cached_molecules
            else:
                bt.logging.debug(f"[Cache INVALID] role={role_mask}: cache is empty, refetching from DB")
    
    # Cache miss or expired — fetch from database
    bt.logging.debug(f"[Cache MISS] role={role_mask}: fetching from database")
    molecules = get_molecules_by_role(role_mask, db_path)
    
    # CRITICAL FIX: Only cache non-empty results
    # Empty results indicate a transient DB error — we should NOT cache those
    if molecules:
        _MOLECULE_POOL_CACHE[role_mask] = (molecules, current_time)
        bt.logging.info(f"[Sampler] Cached {len(molecules)} molecules for role {role_mask}")
    else:
        # DB fetch failed — clear any stale cache for this role
        if role_mask in _MOLECULE_POOL_CACHE:
            del _MOLECULE_POOL_CACHE[role_mask]
        bt.logging.error(f"[Sampler] Failed to fetch molecules for role {role_mask} from DB (cache cleared)")
    
    return molecules


def generate_valid_random_molecules_batch(rxn_id: int, n_samples: int, db_path: str, subnet_config: dict, 
                                 batch_size: int = 200, seed: int = None) -> dict:
    """
    Efficiently generate n_samples valid molecules by generating them in batches and validating.
    
    Args:
        rxn_id: The reaction ID to use
        n_samples: Number of valid molecules to generate
        db_path: Path to the molecules database
        subnet_config: Configuration for validation
        batch_size: Number of molecules to generate per batch
        seed: Random seed (optional)
        
    Returns:
        Dict of molecules for the sampler uid=0
    """
    # Pre-fetch molecule pools to avoid repeated database queries
    reaction_info = get_reaction_info(rxn_id, db_path)
    if not reaction_info:
        bt.logging.error(f"Could not get reaction info for rxn_id {rxn_id}")
        return {"molecules": [None] * n_samples}
    
    smarts, roleA, roleB, roleC = reaction_info
    is_three_component = roleC is not None and roleC != 0
    
    # CRITICAL FIX: Always fetch fresh from database on first call (no cache)
    # The cache was causing 0/250 cycles — force fresh loads to ensure molecule pools exist
    bt.logging.info(f"[Sampler] Fetching molecules for reaction {rxn_id}: roleA={roleA}, roleB={roleB}, roleC={roleC}, is_three={is_three_component}")
    molecules_A = get_molecules_by_role(roleA, db_path)
    molecules_B = get_molecules_by_role(roleB, db_path)
    molecules_C = get_molecules_by_role(roleC, db_path) if is_three_component else []
    
    bt.logging.info(f"[Sampler] Query results: {len(molecules_A)} for roleA, {len(molecules_B)} for roleB, {len(molecules_C) if is_three_component else 'N/A'} for roleC")
    
    # Cache the result for this call only (cache won't persist across iterations)
    if molecules_A:
        _MOLECULE_POOL_CACHE[roleA] = (molecules_A, time.time())
    if molecules_B:
        _MOLECULE_POOL_CACHE[roleB] = (molecules_B, time.time())
    if molecules_C and is_three_component:
        _MOLECULE_POOL_CACHE[roleC] = (molecules_C, time.time())
    
    # CRITICAL: Validate that we have molecule pools — prevent silent failures with empty pools
    if not molecules_A or not molecules_B or (is_three_component and not molecules_C):
        bt.logging.error(f"[CRITICAL] No molecules found for reaction {rxn_id}:")
        bt.logging.error(f"  roleA={roleA}: {len(molecules_A)} molecules")
        bt.logging.error(f"  roleB={roleB}: {len(molecules_B)} molecules")
        if is_three_component:
            bt.logging.error(f"  roleC={roleC}: {len(molecules_C)} molecules")
        
        # Database query already failed above — don't waste time retrying
        bt.logging.error(f"[FATAL] Cannot generate molecules for reaction {rxn_id} — molecule pools are empty")
        return {"molecules": [None] * n_samples}
    
    valid_molecules = []
    seen_keys = set()
    iteration = 0
    MAX_ITERATIONS = 5000  # Increased from 1000 to handle low-yield reactions (was causing 0/250 after batch 1)

    progress_bar = tqdm(total=n_samples, desc="Creating valid molecules", unit="molecule")
    
    while len(valid_molecules) < n_samples and iteration < MAX_ITERATIONS:
        iteration += 1
        
        # Defensive check: if molecule pools are suddenly empty, this indicates a critical failure
        # Abort rather than silently returning 0/250
        if not molecules_A or not molecules_B or (is_three_component and not molecules_C):
            bt.logging.error(f"[CRITICAL] Molecule pools became empty at iteration {iteration} (this should never happen)")
            break
        
        # Calculate how many molecules we still need
        needed = n_samples - len(valid_molecules)
        
        # Generate a batch of molecules (with some buffer for validation failures)
        batch_size_actual = min(batch_size, needed * 2)  # Generate 2x what we need to account for failures
        

        batch_names, batch_smiles = generate_names_and_smiles_from_pools(
            rxn_id, batch_size_actual,
            molecules_A, molecules_B, molecules_C, is_three_component,
            smarts, roleA, roleB, roleC, seed
        )
        
        # Validate directly on SMILES (no DB calls)
        batch_valid_molecules, batch_valid_smiles = validate_smiles_sampler(batch_names, batch_smiles, subnet_config)

        # Deduplicate inside the batch (keep first per InChIKey)
        identical = find_chemically_identical(batch_valid_smiles)
        skip_indices = set()
        for indices in identical.values():
            for j in indices[1:]:
                skip_indices.add(j)

        # Add only chemically-unique molecules across all batches
        added = 0
        for i, name in enumerate(batch_valid_molecules):
            if i in skip_indices or not name:
                continue
            s = batch_valid_smiles[i] if i < len(batch_valid_smiles) else None
            if not s:
                continue
            try:
                mol = Chem.MolFromSmiles(s)
                if not mol:
                    continue
                key = Chem.MolToInchiKey(mol)
            except Exception:
                continue
            if key in seen_keys:
                continue

            seen_keys.add(key)
            valid_molecules.append(name)
            added += 1

            # Stop as soon as we reach the target; attempts_total already includes this attempt
            if len(valid_molecules) >= n_samples:
                break
        
        progress_bar.update(added)
        
        # Check if we've hit max iterations
        if iteration >= MAX_ITERATIONS:
            bt.logging.warning(f"[Sampler] Reached maximum iterations ({MAX_ITERATIONS}), stopping generation with {len(valid_molecules)}/{n_samples} molecules")
            break
    
    # Trim to exact number requested
    final_molecules = valid_molecules[:n_samples]
    progress_bar.close()

    return {"molecules": final_molecules}


def generate_molecules_from_pools(rxn_id: int, n: int, molecules_A: List[Tuple], molecules_B: List[Tuple], 
                                molecules_C: List[Tuple], is_three_component: bool, seed: int = None) -> List[str]:
    """
    Generate molecules using pre-fetched molecule pools to avoid database queries.
    
    Args:
        rxn_id: The reaction ID
        n: Number of molecules to generate
        molecules_A, molecules_B, molecules_C: Pre-fetched molecule pools
        is_three_component: Whether this is a 3-component reaction
        seed: Random seed for reproducibility (optional)
    Returns:
        List of molecule names
    """
    mol_ids = []

    if seed is not None:
        random.seed(seed)
    
    for i in range(n):
        try:
            # Randomly select molecules for each role
            mol_A = random.choice(molecules_A)
            mol_B = random.choice(molecules_B)
            
            mol_id_A, smiles_A, role_mask_A = mol_A
            mol_id_B, smiles_B, role_mask_B = mol_B
            
            if is_three_component:
                mol_C = random.choice(molecules_C)
                mol_id_C, smiles_C, role_mask_C = mol_C
                product_name = f"rxn:{rxn_id}:{mol_id_A}:{mol_id_B}:{mol_id_C}"
            else:
                product_name = f"rxn:{rxn_id}:{mol_id_A}:{mol_id_B}"
            
            mol_ids.append(product_name)
            
        except Exception as e:
            bt.logging.error(f"Error generating molecule {i+1}/{n}: {e}")
            mol_ids.append(None)
    
    return mol_ids


def run_sampler(n_samples: int = 1000, 
                seed: int = None, 
                subnet_config: dict = None, 
                output_path: str = None, 
                save_to_file: bool = False,
                db_path: str = None):
    reactions = get_available_reactions(db_path)
    if not reactions:
        bt.logging.error("No reactions found in the database, check db path and integrity.")
        return

    # CRITICAL FIX: Filter out reactions with 'N/A' SMARTS (database stubs)
    # These are placeholder reactions that cannot be parsed by RDKit
    # Only use reactions with valid SMARTS strings
    valid_reactions = [r for r in reactions if r[1] != 'N/A']
    
    if not valid_reactions:
        bt.logging.error("No valid reactions found! All reactions have 'N/A' SMARTS.")
        return
    
    if len(valid_reactions) < len(reactions):
        bt.logging.warning(f"Filtered out {len(reactions) - len(valid_reactions)} invalid reactions with 'N/A' SMARTS")
    
    reactions = valid_reactions
    rxn_ids = [reactions[i][0] for i in range(len(reactions))]

    # Handle both "allowed_reaction" (mainnet) and "random_valid_reaction" (testnet)
    if "allowed_reaction" in subnet_config:
        rxn_id = int(subnet_config["allowed_reaction"].split(":")[-1])
    elif subnet_config.get("random_valid_reaction", False):
        # Pick a random reaction from available ones
        import random
        if seed is not None:
            random.seed(seed)
        rxn_id = random.choice(rxn_ids)
    else:
        # Default to first available reaction
        rxn_id = rxn_ids[0]
    
    bt.logging.info(f"Generating {n_samples} random molecules for reaction {rxn_id}")

    # Generate molecules with validation in batches for efficiency
    sampler_data = generate_valid_random_molecules_batch(
        rxn_id, n_samples, db_path, subnet_config, batch_size=200, seed=seed
        )

    if save_to_file:
        with open(output_path, "w") as f:
            json.dump(sampler_data, f, ensure_ascii=False, indent=2)

    return sampler_data


if __name__ == "__main__":
    run_sampler()
