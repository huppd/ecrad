#!/usr/bin/env python3

import click
from numpy.lib.npyio import mafromtxt
import pandas as pd
from numpy import genfromtxt
from statistics import mean, stdev

@click.command()
@click.argument('outputfile')
def cli(outputfile):
    ref = pd.read_csv(
        'omptiming.txt',
        delim_whitespace=True,
        names=["routine", "time", "var", "ncalls"],
        # index_col=0,
    )
    
    print()
    mean_time = ref['time'].median()
    stde_time = ref['time'].mad()
    runs = ref['time'].size

    print('median[s]:                  ', mean_time)
    print('mean absolute deviation[s]: ', stde_time)
    print('mean absolute deviation[%]: ', stde_time/mean_time*100)
    print('runs:                       ', runs)
    f = open(outputfile, "w")
    f.write('Solver MCICA part                                                               '+str(mean_time)+ ' '+str(stde_time)+' '+str(runs)+'\n')

if __name__ == '__main__':
    cli()