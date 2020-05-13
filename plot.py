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

# Plot object control
LineWidth=[3.5]
LineColor=['red']
InvertLineOrigin=False # Default origin is bottom, set to True for spectrum with lines originating from top (e.g. IR)

ContinuousLineStyle=['solid']
ContinuousLineWidth=[3.5]
ContinuousColor=['blue']
ContinuousLegend=['experiment']

# Plot style
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
    parser.add_argument('-l','--LineSpectrum', nargs='*', help='line spectrum file')
    parser.add_argument('-c','--ContinuousSpectrum', nargs='*', help='line spectrum file')
    parser.add_argument('-t','--title', type=str, default='Photoelectron spectrum', help='default = Photoelectron spectrum')
    parser.add_argument('-lx','--label_x', type=str, default='Electron kinetic energy / eV', help='default = x')
    parser.add_argument('-ly','--label_y', type=str, default='Relative intensity', help='default = y')
    args = parser.parse_args()
    if args.LineSpectrum is not None:
        for i in range(len(args.LineSpectrum)):
            args.LineSpectrum[i] = Path(args.LineSpectrum[i])
    if args.ContinuousSpectrum is not None:
        for i in range(len(args.ContinuousSpectrum)):
            args.ContinuousSpectrum[i] = Path(args.ContinuousSpectrum[i])
    return args

if __name__ == "__main__":
    args = parse_args()
    # Read data, create artists
    if args.LineSpectrum is not None:
        for j in range(len(args.LineSpectrum)):
            with open(args.LineSpectrum[j],'r') as f:
                data=f.readlines()
                nline=len(data); xline=numpy.empty(nline); yline=numpy.empty(nline)
                for i in range(nline):
                    temp=data[i].split()
                    xline[i]=float(temp[0].strip()); yline[i]=float(temp[1].strip())
            if(InvertLineOrigin):
                for i in range(1,nline):
                    if xline[i]>=xleft and xline[i]<=xright:
                        plt.vlines(xline[i],yline[i],yup,lw=LineWidth[j],color=LineColor[j])
            else:
                for i in range(1,nline):
                    if xline[i]>=xleft and xline[i]<=xright:
                        plt.vlines(xline[i],ylow,yline[i],lw=LineWidth[j],color=LineColor[j])
    if args.ContinuousSpectrum is not None:
        for j in range(len(args.ContinuousSpectrum)):
            with open(args.ContinuousSpectrum[j],'r') as f:
                data=f.readlines()
                ncontinuous=len(data); xcontinuous=numpy.empty(ncontinuous); ycontinuous=numpy.empty(ncontinuous)
                searchleft=True; indexleft=0; searchright=True; indexright=ncontinuous
                for i in range(ncontinuous):
                    temp=data[i].split()
                    xcontinuous[i]=float(temp[0].strip()); ycontinuous[i]=float(temp[1].strip())
                    if(searchleft and xcontinuous[i]>=xleft):
                        searchleft=False; indexleft=i
                    if(searchright and xcontinuous[i]>xright):
                        searchright=False; indexright=i
            plt.plot(xcontinuous[indexleft:indexright],ycontinuous[indexleft:indexright],ls=ContinuousLineStyle[j],lw=ContinuousLineWidth[j],color=ContinuousColor[j],label=ContinuousLegend[j])
    
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
    yloc=ax.get_xticks(minor=True); temp=yloc[1]-yloc[0]; plt.ylim(ylow -temp,yup   +temp) # to look prettier
    
    plt.legend(loc=LegendLocation,frameon=LegendFrame,fontsize=LegendFontSize)
    
    plt.show()