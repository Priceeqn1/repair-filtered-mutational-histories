"""Run archived Paper II research checks, each in a clean Python subprocess."""
from pathlib import Path
import argparse, subprocess, sys
ROOT=Path(__file__).resolve().parent
QUICK=[
    'reproduce/verify_generator.py',
    'reproduce/verify_population_switch.py',
    'reproduce/verify_balance.py',
    'reproduce/verify_common_green.py',
    'reproduce/verify_large_M_independent.py',
    'reproduce/verify_sign_classification.py',
    'reproduce/verify_heterogeneous.py',
]
EXTRA=[
    'reproduce/certify_M3_global.py',
    'checks/certify_M3_with_bracket.py',
    'reproduce/certify_M10_parity.py',
    'checks/verify_all_eight_published.py',
    'reproduce/sensitivity/certify_moderate.py',
    'reproduce/sensitivity/certify_lower.py',
]
parser=argparse.ArgumentParser()
parser.add_argument('--full',action='store_true',help='Also run exact symbolic/rational certificate scripts.')
args=parser.parse_args()
for file in QUICK + (EXTRA if args.full else []):
    print('\n===',file,'===',flush=True)
    subprocess.run([sys.executable,str(ROOT/file)],cwd=ROOT,check=True)
print('\nAll requested checks completed without assertion failures.')
