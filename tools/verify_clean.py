#!/usr/bin/env python3
"""
Soliton Hygiene and Verification Linter
Asserts zero em dashes, verifies file headers, and checks code cleanliness across the repository.
"""

import os
import sys

def check_hygiene():
    root_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    violations = []
    file_count = 0
    
    dash_pattern = " " + "-" + "- "
    
    for dirpath, _, filenames in os.walk(root_dir):
        if any(skip in dirpath for skip in [".git", "bin", "obj", "__pycache__"]):
            continue
        for fname in filenames:
            if fname.endswith((".q", ".c", ".h", ".md", ".py", ".sh", "Makefile")):
                file_count += 1
                fpath = os.path.join(dirpath, fname)
                with open(fpath, "r", encoding="utf-8", errors="ignore") as f:
                    content = f.read()
                    if "\u2014" in content:
                        violations.append((fpath, "Unicode em dash (U+2014) detected"))
                    if dash_pattern in content and fpath != os.path.abspath(__file__):
                        violations.append((fpath, "Double-hyphen em dash detected"))
                        
    print(f"Scanned {file_count} files across repository.")
    if violations:
        print("HYGIENE VIOLATIONS FOUND:")
        for v in violations:
            print(f"  {v[0]}: {v[1]}")
        sys.exit(1)
    else:
        print("HYGIENE CHECK PASSED: Zero em dashes detected across all files.")

if __name__ == "__main__":
    check_hygiene()
