'''
There usually are 2 kinds of spectrums to plot:
1. Line spectrum, obtained from time independent Schrodinger equation
2. Continuous spectrum, obtained from experiment or broadening line spectrum or dynamics

File format: 2 columns to be (x, y)
'''

''' advanced input '''
# Useful external links:
#    https://matplotlib.org/api/text_api.html, for font
#    https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.pyplot.legend.html, especially for loc

# Plot style
InvertLineOrigin=False # Default origin is bottom, set to True for spectrum with lines originating from top (e.g. IR)
LineWidth=3.5
LineColor=['red']

ContinuousLineWidth=3.5
ContinuousLineStyle=['solid']
ContinuousColor=['black']
ContinuousLegend=['experiment']

TitleFontSize=32; TitleFontWeight='regular'
LabelFontSize=36; LabelFontWeight='bold'
AxisThickness=3
TickFontSize=18; TickFontWeight='bold'; TickDirection='in'; TickLength=10; TickWidth=2
LegendLocation=0; LegendFrame=False; LegendFontSize=27

xleft=0; xright=1.2; ylow=0; yup=1 # Plot range

import argparse
from pathlib import Path
import numpy
import matplotlib.pyplot as plt

def parse_args() -> argparse.Namespace: # Command line input
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('-line','--line_spectrum', nargs='*', help='line spectrum file')
    parser.add_argument('-lc','--line_color', nargs='*', help='color of each line spectrum')
    parser.add_argument('-cont','--continuous_spectrum', nargs='*', help='continuous spectrum file')
    parser.add_argument('-cc','--continuous_color', nargs='*', help='color of each continuous spectrum')
    parser.add_argument('-cl','--continuous_legend', nargs='*', help='legend of each continuous spectrum')
    parser.add_argument('-t','--title', type=str, default='Photoelectron spectrum', help='default = Photoelectron spectrum')
    parser.add_argument('-lx','--label_x', type=str, default='Electron kinetic energy / eV', help='default = x')
    parser.add_argument('-ly','--label_y', type=str, default='Relative intensity', help='default = y')
    args = parser.parse_args()
    # Sanity check and some preprocess
    if args.line_spectrum is not None:
        for i in range(len(args.line_spectrum)):
            args.line_spectrum[i] = Path(args.line_spectrum[i])
        if args.line_color is not None:
            assert len(args.line_color)==len(args.line_spectrum), "one color per line spectrum"
    if args.continuous_spectrum is not None:
        for i in range(len(args.continuous_spectrum)):
            args.continuous_spectrum[i] = Path(args.continuous_spectrum[i])
        if args.continuous_color is not None:
            assert len(args.continuous_color)==len(args.continuous_spectrum), "one color per continuous spectrum"
    return args

if __name__ == "__main__":
    args = parse_args()
    # Read data, create artists
    if args.line_spectrum is not None:
        # Adjust plot style control
        if args.line_color is not None:
            LineColor = args.line_color
        else:
            for j in range(1, len(args.line_spectrum)): LineColor.append(LineColor[0])
        # Plot each spectrum
        for j in range(len(args.line_spectrum)):
            with open(args.line_spectrum[j],'r') as f: lines=f.readlines()
            nline=len(lines); xline=numpy.empty(nline); yline=numpy.empty(nline)
            for i in range(nline):
                temp=lines[i].split()
                xline[i]=float(temp[0]); yline[i]=float(temp[1])
            if(InvertLineOrigin):
                for i in range(nline):
                    if xline[i]>=xleft and xline[i]<=xright:
                        plt.vlines(xline[i],yline[i],yup,lw=LineWidth,color=LineColor[j])
            else:
                for i in range(nline):
                    if xline[i]>=xleft and xline[i]<=xright:
                        plt.vlines(xline[i],ylow,yline[i],lw=LineWidth,color=LineColor[j])
    if args.continuous_spectrum is not None:
        # Adjust plot style control
        if args.continuous_color is not None:
            ContinuousColor = args.continuous_color
        else:
            for j in range(1, len(args.continuous_spectrum)): ContinuousColor.append(ContinuousColor[0])
        if args.continuous_legend is not None:
            ContinuousLegend = args.continuous_legend
        else:
            for j in range(1, len(args.continuous_spectrum)): ContinuousLegend.append(ContinuousLegend[0])
        # Plot each spectrum
        for j in range(len(args.continuous_spectrum)):
            with open(args.continuous_spectrum[j],'r') as f: lines=f.readlines()
            ncontinuous=len(lines); xcontinuous=numpy.empty(ncontinuous); ycontinuous=numpy.empty(ncontinuous)
            searchleft=True; indexleft=0; searchright=True; indexright=ncontinuous
            for i in range(ncontinuous):
                temp=lines[i].split()
                xcontinuous[i]=float(temp[0].strip()); ycontinuous[i]=float(temp[1].strip())
                if(searchleft and xcontinuous[i]>=xleft):
                    searchleft=False; indexleft=i
                if(searchright and xcontinuous[i]>xright):
                    searchright=False; indexright=i-1
            plt.plot(xcontinuous[indexleft:indexright],ycontinuous[indexleft:indexright],ls=ContinuousLineStyle[j],lw=ContinuousLineWidth,color=ContinuousColor[j],label=ContinuousLegend[j])
    
    ax=plt.gca() # Adjust plot style
    
    plt.title(args.title,fontsize=TitleFontSize,fontweight=TitleFontWeight)
    
    plt.xlabel(args.label_x,fontsize=LabelFontSize,fontweight=LabelFontWeight)
    plt.ylabel(args.label_y,fontsize=LabelFontSize,fontweight=LabelFontWeight)
    
    plt.xlim(xleft,xright); plt.ylim(ylow,yup)
    
    ax.spines['left'].set_linewidth(AxisThickness); ax.spines['right'].set_linewidth(AxisThickness)
    ax.spines['bottom'].set_linewidth(AxisThickness); ax.spines['top'].set_linewidth(AxisThickness)
    
    plt.xticks(fontsize=TickFontSize,fontweight=TickFontWeight)
    plt.yticks(fontsize=TickFontSize,fontweight=TickFontWeight)
    plt.minorticks_on()
    ax.tick_params(which='major',direction=TickDirection,length=  TickLength  ,width=TickWidth,top=True,right=True)
    ax.tick_params(which='minor',direction=TickDirection,length=0.5*TickLength,width=TickWidth,top=True,right=True)
    xloc=ax.get_xticks(minor=True); temp=xloc[1]-xloc[0]; plt.xlim(xleft-temp,xright+temp) # Add a minor tick
    yloc=ax.get_yticks(minor=True); temp=yloc[1]-yloc[0]; plt.ylim(ylow -temp,yup   +temp) # to look prettier
    
    plt.legend(loc=LegendLocation,frameon=LegendFrame,fontsize=LegendFontSize)
    
    plt.show()
