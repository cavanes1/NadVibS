'''
Extract the line spectrum from output.dat

File format: 2 columns to be (x, y)
'''

import argparse
from pathlib import Path
import numpy

def parse_args() -> argparse.Namespace: # Command line input
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('input', type=Path, help='output.dat to extract')
    parser.add_argument('-o','--output', type=Path, default=Path('spectrum.txt'), help='output file (default = spectrum.txt)')
    parser.add_argument('-t','--threshold', type=float, default=10**-15, help='convergence threshold')
    parser.add_argument('-c','--check_convergence', action='store_true', help='print convergence information')
    parser.add_argument('-i','--invert', action='store_true', help='invert position axis')
    parser.add_argument('-sm','--shift_maximum', type=float, help='shift the maximum to')
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = parse_args()
    # Read
    with open(args.input,'r') as f: lines=f.readlines()
    for i in range(len(lines)):
        if "FIRST eigenvalue:" in lines[i]:
            temp = lines[i].split()
        if "Root  Index  Convergence     E(cm-1)        E(eV)     Intensity" in lines[i]:
            start = i + 1
        if "All requested eigenvalues computed." in lines[i]:
            stop  = i - 1
            break
    n = stop - start
    c = numpy.empty(n)
    x = numpy.empty(n)
    y = numpy.empty(n)
    for i in range(n):
        temp = lines[start + i].split()
        c[i] = float(temp[2].strip())
        x[i] = float(temp[4].strip())
        y[i] = float(temp[5].strip())
    # Modify
    if args.invert:                x  = -x
    if args.shift_maximum != None: x += args.shift_maximum - x[y.argmax()]
    # Output
    with open(args.output,'w') as f:
        for i in range(n):
            if c[i] < args.threshold and y[i] != 0:
                print("%12.5f%12.7f"%(x[i], y[i]), file=f)
                if args.check_convergence: print("%12.5f%12.7f"%(x[i], c[i]))