#!/usr/bin/env python3
"""
Fix database corruption — replace 'N/A' SMARTS with valid reaction templates or remove them.

PROBLEM:
- Reactions 1, 3, 5 have 'N/A' as SMARTS values (placeholders)
- RDKit cannot parse 'N/A' as SMARTS
- Miner crashes in infinite error loop trying to parse 'N/A'

SOLUTION:
- Keep the 'N/A' reactions in the database
- Mark them as "special" reactions that only use hardcoded logic
- Modify random_sampler.py to SKIP 'N/A' reactions in the random sampling
- This way only reactions 2 and 4 (valid SMARTS) are used for normal sampling
- Reactions 1, 3, 5 can still work if validators explicitly request them

Actually, wait — let me check if there are any molecules in the molecules table first.
"""

import sqlite3
import os

db_path = '/home/openclaw/.openclaw/workspace/nova_ml_build/combinatorial_db/molecules.sqlite'
conn = sqlite3.connect(db_path)
cursor = conn.cursor()

# Check molecules table
cursor.execute("SELECT COUNT(*) FROM molecules;")
mol_count = cursor.fetchone()[0]
print(f"Total molecules in database: {mol_count}")

if mol_count == 0:
    print("\n⚠️  CRITICAL: Molecules table is EMPTY!")
    print("   This explains why the miner fails — no molecules to sample from.")
    print("\n   Options:")
    print("   1. Seed the molecules table with valid molecules")
    print("   2. Use a different database with pre-populated molecules")
    print("   3. Load molecules from an external source")
else:
    print(f"✅ Molecules table has {mol_count} entries")
    cursor.execute("SELECT mol_id, smiles, role_mask FROM molecules LIMIT 5;")
    for row in cursor.fetchall():
        print(f"   {row}")

conn.close()

# Now let's fix the sampler to skip 'N/A' reactions
print("\n" + "="*60)
print("APPLYING FIX: Skip 'N/A' reactions in sampler")
print("="*60)

sampler_path = '/home/openclaw/.openclaw/workspace/nova_ml_build/neurons/miner/random_sampler.py'

with open(sampler_path, 'r') as f:
    code = f.read()

# Find the section in iterative_sampling_loop where reactions are sampled
# We need to filter out 'N/A' reactions

if "smarts = reaction_info[0]" in code or "smarts, roleA, roleB" in code:
    print("✅ Found reaction sampling code")
    
    # Add a filter after get_available_reactions()
    if "def iterative_sampling_loop" in code:
        # Find and modify get_available_reactions to filter
        print("✅ Located iterative_sampling_loop")
        
        # We need to add: available_reactions = [r for r in available_reactions if r[1] != 'N/A']
        
        if "available_reactions = get_available_reactions" in code:
            old_line = "available_reactions = get_available_reactions"
            new_code = code.replace(
                "available_reactions = get_available_reactions",
                """# Filter out 'N/A' placeholder reactions (not valid SMARTS)
available_reactions = get_available_reactions"""
            )
            
            # Add the filter right after the get_available_reactions call
            if "available_reactions = get_available_reactions" in new_code:
                # Find the next line and insert our filter
                lines = new_code.split('\n')
                new_lines = []
                for i, line in enumerate(lines):
                    new_lines.append(line)
                    if 'available_reactions = get_available_reactions' in line and i+1 < len(lines):
                        indent = len(line) - len(line.lstrip())
                        # Add filter on next line
                        new_lines.append(' ' * indent + '# Skip reactions with invalid \'N/A\' SMARTS')
                        new_lines.append(' ' * indent + 'available_reactions = [r for r in available_reactions if r[1] != \'N/A\']')
                        new_lines.append(' ' * indent + 'if not available_reactions:')
                        new_lines.append(' ' * (indent + 4) + 'bt.logging.error("No valid reactions available!")')
                        new_lines.append(' ' * (indent + 4) + 'raise RuntimeError("Database has no valid reactions")')
                
                new_code = '\n'.join(new_lines)
                
                with open(sampler_path, 'w') as f:
                    f.write(new_code)
                
                print("✅ Applied filter to skip 'N/A' reactions")
else:
    print("❌ Could not find reaction sampling code")
