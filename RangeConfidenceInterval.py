"""
 Standardized time-series confidence interval method
 based on range scaling function; see "Using excursions to
 analyze simulation output", J. Calvin, in Probability in the
 Engineering and Informational sciences, Vol. 24, 2010 pp.201-218.
 Assume: path passed to confInt is the output of a discrete-event
 simulation that satisfies a functional central limit theorem.
"""

from scipy.special import kn
from scipy.special import iv
import math
import numpy as np

#mylist = [1.0, 3.0, 4.0, 5.0]
#myarray = np.asarray(mylist)
#print ("The new created array is : ", end =" ")
# for i in range (0, 4):
#print (myarray[i], end =" ")


def getCritical(percent):
    lower = 0.01
    upper = 9.0
    while (upper - lower > 0.00000000001):
        x = 0.5 * (lower + upper)
        y = eval(x)
        if (y > percent):
            upper = x
        else:
            lower = x
    return 0.5 * (lower + upper)


def eval(x):
    res = 0.0
    for n in list(range(10000)):
        res = res + (n+1) * kn(1, math.pi * (n+1) * x)
    res = 1.0 - math.pi * x * x * res
    return res

# def static void confInt(double cl, double[] cr, double[] path) {


def confInt(cl, lpath):
    path = np.asarray(lpath)
    max = len(lpath)
    Y = np.asarray(lpath)

    # First normalize the path; cumulative sums tied down at 1.
    Y[0] = path[0]
    for i in range(1, max):
        Y[i] = path[i] + Y[i - 1]
    for i in range(0, max):
        Y[i] = Y[i] / (max - 1.0)
    mu = Y[max - 1]
    for i in range(1, max):
        Y[i] = Y[i] - (i / (max - 1.0)) * Y[max - 1]
    g = 0.0
    sup = 0.0
    inf = 0.0
    for i in range(1, max):
        if (Y[i] > sup):
            sup = Y[i]
        if (Y[i] < inf):
            inf = Y[i]
    g = sup - inf  # Range of normalized path.
    contCor = 2.0 * 0.5826 * math.sqrt(2.0 / math.pi) * g / math.sqrt(max)
    g = g + contCor  # Apply continuity correction.
    confRange = getCritical(cl)

    #System.out.println("95% confidence interval for sojourn time:");
    #cr[0] = mu;
    #cr[1] = g * confRange;
    return mu, g * confRange

    #System.out.println((mu - g * confRange) + ", " + (mu + g * confRange));
