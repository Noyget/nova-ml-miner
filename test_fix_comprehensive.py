#!/usr/bin/env python3
"""
Comprehensive test to verify the 0/250 bug is fixed
Tests multiple iterations to ensure no alternating failures
"""
import sys
import os
import json
import logging
import time

logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(name)s] %(message)s')
logger = logging.getLogger("TEST")

# Setup
BASE_DIR = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(BASE_DIR, "neurons", "miner"))
sys.path.insert(0, os.path.join(BASE_DIR, "neurons"))
sys.path.insert(0, BASE_DIR)

DB_PATH = os.path.join(BASE_DIR, "combinatorial_db", "molecules.sqlite")

def test_fix():
    """Run sampler for 10+ iterations checking for 0/250 patterns"""
    from neurons.miner.random_sampler import run_sampler
    
    logger.info("=" * 80)
    logger.info("TEST: Running sampler for 15 iterations to check for 0/250 cycles")
    logger.info("=" * 80)
    
    config = {
        "num_molecules": 100,
        "antitarget_weight": 1.0,
        "min_heavy_atoms": 20,
        "min_rotatable_bonds": 2,
        "max_rotatable_bonds": 10,
        "target_sequences": {"1ABC": "MKTII"},
        "antitarget_sequences": {"2XYZ": "MSTIV"},
        "allowed_reaction": "rxn:2",  # Format: "rxn:N" where N is the reaction ID
        "enable_scoring": False,
    }
    
    results = []
    for iteration in range(15):
        logger.info(f"\n[Iteration {iteration + 1}/15]")
        
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
                num_molecules = len(result["molecules"])
                num_valid = len([m for m in result["molecules"] if m])
                
                logger.info(f"  Generated: {num_valid}/{num_molecules} molecules in {elapsed:.2f}s")
                results.append((iteration + 1, num_valid, elapsed))
                
                if num_valid == 0:
                    logger.error(f"  ❌ BUG REPRODUCED: 0/250 molecules!")
                    # Print first few for debugging
                    for i in range(min(5, num_molecules)):
                        logger.error(f"    molecules[{i}] = {result['molecules'][i]}")
                elif num_valid < 100:
                    logger.warning(f"  ⚠️  Low yield: only {num_valid}/250 molecules")
                else:
                    logger.info(f"  ✅ Good yield: {num_valid}/250 molecules")
            else:
                logger.error(f"  ❌ No result or missing 'molecules' key")
                results.append((iteration + 1, 0, elapsed))
        except Exception as e:
            logger.error(f"  ❌ Error: {e}")
            import traceback
            traceback.print_exc()
            results.append((iteration + 1, -1, 0))
    
    # Summary
    logger.info("\n" + "=" * 80)
    logger.info("SUMMARY")
    logger.info("=" * 80)
    logger.info(f"{'Iter':<6} {'Molecules':<12} {'Time (s)':<10}")
    logger.info("-" * 80)
    
    zero_count = 0
    low_count = 0
    good_count = 0
    
    for iter_num, count, elapsed in results:
        status = "✅ GOOD" if count > 100 else "⚠️  LOW" if count > 0 else "❌ ZERO"
        logger.info(f"{iter_num:<6} {count:<12} {elapsed:<10.2f}  {status}")
        
        if count == 0:
            zero_count += 1
        elif count < 100:
            low_count += 1
        else:
            good_count += 1
    
    logger.info("-" * 80)
    logger.info(f"RESULTS: {good_count} good, {low_count} low, {zero_count} zero")
    logger.info("=" * 80)
    
    if zero_count > 0:
        logger.error(f"\n❌ BUG STILL EXISTS: {zero_count} iterations had 0/250 molecules")
        return False
    elif low_count > 2:
        logger.warning(f"\n⚠️  Some instability detected: {low_count} iterations had low yields")
        return False
    else:
        logger.info(f"\n✅ SUCCESS: Bug appears to be fixed!")
        logger.info(f"   All {len(results)} iterations produced >100 molecules")
        return True

if __name__ == "__main__":
    try:
        success = test_fix()
        sys.exit(0 if success else 1)
    except KeyboardInterrupt:
        logger.info("Interrupted")
        sys.exit(0)
    except Exception as e:
        logger.error(f"Fatal error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
