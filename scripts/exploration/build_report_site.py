"""Build the unified offline research portal and all actual molecular views."""
import json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from mechbridge.report_site import build_site

if __name__=='__main__':print(json.dumps(build_site(),indent=2))
