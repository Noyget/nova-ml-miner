#!/usr/bin/env python3
"""Simple test to confirm the cache bug"""

_CACHE_TTL_SECONDS = 0
_CACHE = {1: ("data", 1234.5)}

import time

current_time = time.time()
role_mask = 1

# Simulate the buggy logic
if role_mask in _CACHE:
    cached_data, cache_time = _CACHE[role_mask]
    time_diff = current_time - cache_time
    print(f"Cache exists for role {role_mask}")
    print(f"  Cache time: {cache_time}")
    print(f"  Current time: {current_time}")
    print(f"  Time diff: {time_diff}")
    print(f"  TTL: {_CACHE_TTL_SECONDS}")
    print(f"  Condition (time_diff < TTL): {time_diff < _CACHE_TTL_SECONDS}")
    print(f"  → Would return cached data: {time_diff < _CACHE_TTL_SECONDS}")
    print()
    print("BUG: The condition is ALWAYS FALSE when TTL=0!")
    print("So the code always goes to 'Cache miss' branch and tries to fetch from DB")
    print("But if DB query fails or is slow, it returns empty list, corrupting molecule pools")
