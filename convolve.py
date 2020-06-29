'''
Convolve the line spectrum

File format: 2 columns to be (x, y)
'''

import argparse
from pathlib import Path
import math
import numpy
import matplotlib.pyplot as plt

def parse_args() -> argparse.Namespace: # Command line input
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('start', type=float, help='start of the convolutional spectrum')
    parser.add_argument('stop' , type=float, help=' stop of the convolutional spectrum')
    parser.add_argument('input', type=Path, help='spectrum to convolve')
    parser.add_argument('aperture', type=str, help='Gaussian or Lorentzian')
    parser.add_argument('sigma', type=float, help='width of aperture')
    parser.add_argument('-sm','--shift_maximum', type=float, help='shift the convolution maximum to')
    parser.add_argument('-o','--output', type=Path, default=Path('convolution.txt'), help='output file (default = convolution.txt)')
    parser.add_argument('-r','--reference', type=Path, help='reference spectrum file to compare with')
    parser.add_argument('-p','--plot', action='store_true', help='plot the convolutional spectrum')
    parser.add_argument('-pr','--plot_range', type=float, nargs=2, help='x-axis plot range')
    args = parser.parse_args()
    return args

def Gaussian(x:float, miu:float, sigma:float) -> float:
    temp = (x-miu) / sigma
    y = numpy.exp(-0.5 * temp*temp)
    return y

def Lorentzian(x:float, miu:float, sigma:float) -> float:
    temp = (x-miu) / sigma
    y = 1.0/sigma / (1.0 + temp*temp)
    return y

aperture_dict = {'Gaussian':Gaussian, 'Lorentzian':Lorentzian}

def convolve(aperture:'function(x,miu,sigma)',
x:numpy.ndarray, y:numpy.ndarray, sigma:float,
grid:numpy.ndarray) -> numpy.ndarray:
    envelope = numpy.zeros(grid.shape)
    for i in range(grid.shape[0]):
        for j in range(x.shape[0]):
            envelope[i] += y[j] * aperture(grid[i], x[j], sigma)
    envelope /= numpy.max(envelope)
    return envelope

if __name__ == "__main__":
    # Initialize
    args = parse_args()
    with open(args.input,'r') as f: lines = f.readlines()
    x = numpy.empty(len(lines)); y = numpy.empty(len(lines))
    for i in range(len(lines)):
        temp = lines[i].split()
        x[i] = float(temp[0].strip()); y[i] = float(temp[1].strip())
    grid = numpy.linspace(args.start, args.stop, num=1000)
    # Do the job
    envelope = convolve(aperture_dict[args.aperture], x, y, args.sigma, grid)
    if args.shift_maximum != None:
        shift = args.shift_maximum - grid[envelope.argmax()]
        grid += shift
        x    += shift
        print('The peaks have been shifted by ', shift)
    # Output
    with open(args.output,'w') as f:
        for i in range(grid.shape[0]):
            print(grid[i], envelope[i], sep='    ', file=f)
    # Optionally plot
    if args.reference != None or args.plot:
        if args.reference != None:
            with open(args.reference,'r') as f: lines = f.readlines()
            xc = numpy.empty(len(lines)); yc = numpy.empty(len(lines))
            for i in range(len(lines)):
                temp = lines[i].split()
                xc[i] = float(temp[0].strip()); yc[i] = float(temp[1].strip())
            plt.plot(xc, yc, color='black', label='reference')
        if args.plot_range != None:
            left  = args.plot_range[0]
            right = args.plot_range[1]
        else:
            if args.reference == None:
                left  = args.start
                right = args.stop
            else:
                left  = xc[0]
                right = xc[xc.shape[0]-1]
        plt.plot(grid, envelope, color='blue', label='convolution')
        for i in range(x.shape[0]):
            if x[i] >= left and x[i] <= right:
                plt.vlines(x[i], 0.0, y[i], color='red')
        plt.xlim(left, right)
        plt.legend(frameon=False)
        plt.show()
    