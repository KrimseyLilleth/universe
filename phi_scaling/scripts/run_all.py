"""Run every fetch script in order. Each is independent; a failure in one
does not abort the others."""
from __future__ import annotations

import importlib
import sys
import traceback

MODULES = [
    "01_nuclear_binding",
    "02_ionization_energies",
    "03_shell_closures",
    "04_bond_dissociation",
    "05_orbital_resonances",
    "06_stellar_imf",
    "07_biological_scales",
]


def main() -> int:
    sys.path.insert(0, str(__import__("pathlib").Path(__file__).parent))
    failures = 0
    for name in MODULES:
        try:
            mod = importlib.import_module(name)
            mod.main()
        except Exception:  # noqa: BLE001
            failures += 1
            traceback.print_exc()
    print(f"\nDone. {len(MODULES) - failures}/{len(MODULES)} scripts succeeded.")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
