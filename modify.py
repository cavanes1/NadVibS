'''
Renormalize, rescale, invert, shift the spectrum
Throw out the lines with 0 intensity

File format: 2 columns to be (x, y)
'''

import argparse
from pathlib import Path
import numpy

def parse_args() -> argparse.Namespace: # Command line input
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('input', type=Path, help='spectrum to modify')
    parser.add_argument('-o','--output', type=Path, default=Path('output.txt'), help='output file (default = output.txt)')
    parser.add_argument('-n','--normalize', action='store_true', help='normalize intensity to let max = 1')
    parser.add_argument('-r','--rescale', type=float, default=1.0, help='rescale intensity (default = 1)')
    parser.add_argument('-i','--invert', action='store_true', help='invert position axis')
    parser.add_argument('-s','--shift', type=float, default=0.0, help='shift position, if(invert): shift after inversion (default = 0)')
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = parse_args()
    # Read
    with open(args.input,'r') as f: lines=f.readlines()
    n=len(lines); x=numpy.empty(n); y=numpy.empty(n)
    for i in range(n):
        temp=lines[i].split()
        x[i]=float(temp[0].strip()); y[i]=float(temp[1].strip())
    # Modify
    if args.normalize:      y /= numpy.max(y)
    if args.rescale != 1.0: y *= args.rescale
    if args.invert:         x  = -x
    if args.shift != 0.0:   x += args.shift
    # Output
    with open(args.output,'w') as f:
        for i in range(n):
            if(y[i]!=0): print(x[i],y[i],sep='    ',file=f)