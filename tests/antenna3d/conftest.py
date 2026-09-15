import sys
from pathlib import Path

_HERE = Path(__file__).resolve().parent          # tests/antenna3d
_ROOT = _HERE.parent.parent                      # the repository root

# So `phantom` imports as a module, and so pytest works from a checkout that
# has not been installed.
sys.path.insert(0, str(_ROOT))
sys.path.insert(0, str(_HERE))
