#!/usr/bin/env python3
"""
Diagnostic test for NOVA miner 0/250 cycle bug
Tests for: database issues, SMARTS validation failures, reaction info corruption
"""
import sys
import os
import json
import time
import logging

# Minimal setup
logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] %(message)s')
logger = logging.getLogger("DIAGNOSTIC")

# Add paths
BASE_DIR = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(BASE_DIR, "neurons", "miner"))
sys.path.insert(0, os.path.join(BASE_DIR, "neurons"))
sys.path.insert(0, BASE_DIR)

DB_PATH = os.path.join(BASE_DIR, "combinatorial_db", "molecules.sqlite")

def test_reaction_info_consistency():
    """Test if get_reaction_info returns consistent results across multiple calls"""
    from combinatorial_db.reactions import get_reaction_info
    
    logger.info("\n=== TEST 1: Reaction Info Consistency ===")
    
    for rxn_id in [1, 2, 3, 4, 5]:
        logger.info(f"\nTesting rxn_id={rxn_id}...")
        results = []
        
        # Call get_reaction_info 20 times
        for i in range(20):
            result = get_reaction_info(rxn_id, DB_PATH)
            results.append(result)
            if i == 0 or result != results[0]:
                logger.info(f"  Call {i+1}: {result}")
        
        # Check consistency
        all_same = all(r == results[0] for r in results)
        if all_same:
            logger.info(f"✅ rxn_id={rxn_id} CONSISTENT across 20 calls")
        else:
            logger.error(f"❌ rxn_id={rxn_id} INCONSISTENT!")
            for i, r in enumerate(results):
                if r != results[0]:
                    logger.error(f"  Call {i+1} DIFFERS: {r} vs {results[0]}")

def test_molecule_pool_freshness():
    """Test if molecule pools are fresh and valid each time"""
    from neurons.miner.random_sampler import get_molecules_by_role, get_available_reactions
    from rdkit import Chem
    
    logger.info("\n=== TEST 2: Molecule Pool Freshness ===")
    
    reactions = get_available_reactions(DB_PATH)
    logger.info(f"Available reactions: {len(reactions)}")
    
    for i, (rxn_id, smarts, roleA, roleB, roleC) in enumerate(reactions):
        logger.info(f"\nReaction {rxn_id}: roles=({roleA},{roleB},{roleC}), smarts={smarts[:30]}...")
        
        # Test roleA molecules
        mols_a = get_molecules_by_role(roleA, DB_PATH)
        logger.info(f"  roleA={roleA}: {len(mols_a)} molecules")
        
        if mols_a:
            # Check first 3 are valid SMILES
            valid_count = 0
            for mol_id, smiles, role_mask in mols_a[:3]:
                try:
                    mol = Chem.MolFromSmiles(smiles)
                    if mol:
                        valid_count += 1
                except:
                    pass
            logger.info(f"    Sampled 3 molecules: {valid_count}/3 valid SMILES")

def test_smarts_reaction_validity():
    """Test if SMARTS reactions are valid and consistent"""
    from combinatorial_db.reactions import perform_smarts_reaction, get_reaction_info
    from neurons.miner.random_sampler import get_molecules_by_role
    from rdkit import Chem
    
    logger.info("\n=== TEST 3: SMARTS Reaction Validity ===")
    
    # Test each reaction with real molecules
    for rxn_id in [1, 2, 4]:  # Skip 3,5 as they require special logic
        reaction_info = get_reaction_info(rxn_id, DB_PATH)
        if not reaction_info:
            logger.error(f"❌ Could not get reaction info for rxn_id={rxn_id}")
            continue
        
        smarts, roleA, roleB, roleC = reaction_info
        logger.info(f"\nReaction {rxn_id}: smarts={smarts}")
        
        # Get molecules
        mols_a = get_molecules_by_role(roleA, DB_PATH)
        mols_b = get_molecules_by_role(roleB, DB_PATH)
        
        if not mols_a or not mols_b:
            logger.error(f"  No molecules found for roles")
            continue
        
        # Try 5 random reactions
        import random
        success_count = 0
        fail_count = 0
        
        for attempt in range(5):
            mol_a = random.choice(mols_a)
            mol_b = random.choice(mols_b)
            
            _, smiles_a, _ = mol_a
            _, smiles_b, _ = mol_b
            
            result = perform_smarts_reaction(smiles_a, smiles_b, smarts)
            if result:
                success_count += 1
                if attempt == 0:
                    logger.info(f"  Attempt 1 SMARTS result: {result[:50]}...")
            else:
                fail_count += 1
        
        logger.info(f"  5 attempts: {success_count} succeeded, {fail_count} failed")

def test_sampler_multiple_iterations():
    """Test if sampler works consistently across multiple iterations"""
    from neurons.miner.random_sampler import run_sampler
    
    logger.info("\n=== TEST 4: Sampler Multiple Iterations ===")
    
    config = {
        "num_molecules": 50,
        "antitarget_weight": 1.0,
        "min_heavy_atoms": 20,
        "min_rotatable_bonds": 2,
        "max_rotatable_bonds": 10,
        "target_sequences": {"1ABC": "MKTII"},
        "antitarget_sequences": {"2XYZ": "MSTIV"},
        "allowed_reaction": "rxn:2:0:0",
    }
    
    for iteration in range(10):
        logger.info(f"\nIteration {iteration + 1}/10...")
        
        start = time.time()
        try:
            result = run_sampler(
                n_samples=250,
                subnet_config=config,
                output_path="/tmp/sampler_test.json",
                save_to_file=False,
                db_path=DB_PATH,
            )
            elapsed = time.time() - start
            
            if result and "molecules" in result:
                count = len([m for m in result["molecules"] if m])
                logger.info(f"  ✅ Generated {count}/250 molecules in {elapsed:.2f}s")
                
                if count == 0:
                    logger.error(f"  ❌ BUG REPRODUCED: 0/250 molecules!")
                    # Print first few sampler results for debugging
                    for i in range(min(5, len(result["molecules"]))):
                        logger.error(f"    result[{i}] = {result['molecules'][i]}")
            else:
                logger.error(f"  ❌ No result or missing 'molecules' key")
        except Exception as e:
            logger.error(f"  ❌ Error: {e}")
            import traceback
            traceback.print_exc()

def test_reaction_info_state_isolation():
    """Test if calling get_reaction_info modifies any global state"""
    from combinatorial_db.reactions import get_reaction_info
    import gc
    
    logger.info("\n=== TEST 5: Reaction Info State Isolation ===")
    
    logger.info("Calling get_reaction_info repeatedly with forced GC...")
    
    for i in range(50):
        rxn = get_reaction_info(1, DB_PATH)
        if i % 10 == 0:
            gc.collect()
            logger.info(f"  Call {i+1}/50: {rxn}")
    
    logger.info("✅ State isolation test complete")

if __name__ == "__main__":
    try:
        test_reaction_info_consistency()
        test_molecule_pool_freshness()
        test_smarts_reaction_validity()
        test_sampler_multiple_iterations()
        test_reaction_info_state_isolation()
        
        logger.info("\n=== DIAGNOSTIC COMPLETE ===")
    except KeyboardInterrupt:
        logger.info("Interrupted")
        sys.exit(0)
    except Exception as e:
        logger.error(f"Fatal error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
