import sys
from pathlib import Path

# The repository root, so `combine` imports. It is a single toolbox-level
# module rather than part of either package, so the editable install does not
# put it on the path.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
