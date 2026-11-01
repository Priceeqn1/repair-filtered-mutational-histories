"""Recreate manuscript figures with their documented Python implementations."""
from pathlib import Path
import subprocess,sys
ROOT=Path(__file__).resolve().parent
SCRIPTS=[
  'reproduce/plot_switches.py',
  'reproduce/plot_opportunity.py',
  'reproduce/plot_reentrant.py',
  'reproduce/plot_threshold_regimes.py',
  'reproduce/plot_biological_sensitivity.py',
]
for file in SCRIPTS:
    print('\n===',file,'===',flush=True)
    subprocess.run([sys.executable,str(ROOT/file)],cwd=ROOT,check=True)
print('\nAll Python figures generated in figures/.')
print('The model schematic is provided separately as editable TikZ: figures/model_schematic.tex')
