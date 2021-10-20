import os
import glob

files = glob.glob('*.csv')
for f in files:
    os.remove(f)