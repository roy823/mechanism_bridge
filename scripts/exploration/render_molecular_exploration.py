"""Draw molecular formulas, 2D structures and actual 3D reaction trajectories."""
import argparse
import json
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from mechbridge.molecular_visuals import render_visuals

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('run',type=Path)
    a=p.parse_args()
    print(json.dumps(render_visuals(a.run),indent=2))
