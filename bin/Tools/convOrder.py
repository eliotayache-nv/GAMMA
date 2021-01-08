# -*- coding: utf-8 -*-
# @Author: Eliot Ayache
# @Date:   2020-01-06 14:05:36
# @Last Modified by:   Eliot Ayache
# @Last Modified time: 2020-01-06 17:37:30


"""
This program fits the convergence order of a given test.

Returns:
- Prints value of order of convergence
"""


# ----------------------------------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
import argparse
import astropy.io.ascii as ascii
from scipy import optimize

# ----------------------------------------------------------------------------------------
def parseArguments():
    """Parses the optional arguments to the program"""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-f","--filename", 
        help="data file (error vs resolution). This file must be a 2 column text file \
        with no header. Left column contains values of resolution and right column \
        contains values of L1 error.", type=str)

    parser.add_argument("-p","--plot", 
        help="plot option. Plots the data and the fit.", 
        action='store_true')

    parser.add_argument("-r","--rates", 
        help="individual rates instead of global fit.", 
        action='store_true')

    # Parse arguments
    args = parser.parse_args()
    return args

# ----------------------------------------------------------------------------------------
def fitFunc(x, *params):
    """
        Linear fit function

        Arguments:
        - x : position to evaluate
        - a : coefficient
        - b : offset

        Returns:
        - ax + b
    """
    # print params
    a, b = params
    return a*x+b

# ----------------------------------------------------------------------------------------
def plotResults(res, err, rates=[], *params):
    """
    Plots the data and the fit.

    Arguments:
    - res : resolution data array (size N)
    - err : error data array (size N)
    - rates=[] ; list of convergence rates.
    - *params : fit function parameters.

    Returns:
    - N/A
    """

    fig, ax1 = plt.subplots()

    color = 'tab:blue'

    ax1.set_xlabel('Resolution (number of zones)')
    ax1.set_ylabel('L1 error', color=color)
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.scatter(res, err, color=color, zorder=1)
    ax1.tick_params(axis='y', labelcolor=color)


    if rates == []:
        x = np.logspace(np.min(res), np.max(res), 100)
        y = np.exp(fitFunc(np.log(x), *params))
        plt.title("Convergence rate = %f" %np.abs(params[0])) 
        ax1.plot(x,y, zorder=1)
    else:
        ax2 = ax1.twinx()
        color = 'tab:red'
        ax2.set_ylabel('Convergence rate', color=color)
        ax2.bar(res[1:], -rates, width=res[1:]/2., color=color, zorder=3, alpha=0.3)
        ax2.tick_params(axis='y', labelcolor=color)

    plt.tight_layout()
    plt.show()

# ----------------------------------------------------------------------------------------
def measureRates(res, err):
    """
    Returns the convergence rates for every gap in the data.

    Arguments:
    - res : resolution data array (size N)
    - err : error data array (size N)
    
    Returns:
    - rates : Array of rates corresponding to res[1:]
    """

    res = np.log(res)
    err = np.log(err)

    rates = (err[1:] - err[:-1]) / (res[1:] - res[:-1]) 
    return rates


# ----------------------------------------------------------------------------------------
# Script starts here
# ----------------------------------------------------------------------------------------
args = parseArguments()

# Check for file:
if args.filename == None:
    print "ERROR! No data file specified!"
    print "Please specify using -f, --filename."
    print "Exiting with error code 1."
    exit(1)

f = ascii.read(args.filename)
print("Loaded file successfully!")
res_data = np.array(f['col1'])
err_data = np.array(f['col2'])

if args.rates:
    rates = measureRates(res_data, err_data)
    print(rates)
else:
    params, params_cov = optimize.curve_fit(fitFunc, np.log(res_data), np.log(err_data),
                                            p0=[-2, 0])

    print("Convergence rate = %f ± %f") %(np.abs(params[0]), params_cov[0,0])

if args.plot:
    if args.rates:
        plotResults(res_data, err_data, rates=rates)
    else:
        plotResults(res_data, err_data, *params)


print("Done!")










