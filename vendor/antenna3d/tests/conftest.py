import sys
from pathlib import Path

# So `pytest` works from the folder without installing, and `phantom` imports as a module.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))
