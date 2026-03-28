#!/usr/bin/env python3
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "neurons", "miner"))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "neurons"))
sys.path.insert(0, os.path.dirname(__file__))

import logging
logging.basicConfig(level=logging.WARNING)

from neurons.miner.random_sampler import run_sampler

DB_PATH = os.path.join(os.path.dirname(__file__), "combinatorial_db", "molecules.sqlite")
config = {
    'num_molecules': 100,
    'antitarget_weight': 1.0,
    'min_heavy_atoms': 20,
    'min_rotatable_bonds': 2,
    'max_rotatable_bonds': 10,
    'target_sequences': {'1ABC': 'MKTII'},
    'antitarget_sequences': {'2XYZ': 'MSTIV'},
    'allowed_reaction': 'rxn:2',
    'enable_scoring': False,
}

print('\n✅ TEST: Running 10 iterations to verify fix')
print('=' * 60)

all_passed = True
for i in range(10):
    result = run_sampler(n_samples=250, subnet_config=config, output_path='/tmp/test.json', save_to_file=False, db_path=DB_PATH)
    num_valid = len([m for m in result['molecules'] if m])
    status = '✅' if num_valid > 200 else '❌' if num_valid == 0 else '⚠️'
    print(f'Iteration {i+1:2d}: {num_valid:3d}/250 molecules {status}')
    if num_valid == 0:
        all_passed = False

print('=' * 60)
if all_passed:
    print('✅ SUCCESS: All 10 iterations passed without 0/250 failures!')
else:
    print('❌ FAILED: Some iterations had 0/250')
    sys.exit(1)
