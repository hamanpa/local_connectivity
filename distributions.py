import numpy as np


def norm_hyperbolic_pdf(x, a, t, b):
    return b*np.exp(-a*np.sqrt(t**2 + x**2))

def hyperbolic_pdf(x, a, t):
    return np.exp(-a*np.sqrt(t**2 + x**2))

def gaussian_pdf(x, a, sigma):
    return a*np.exp(-x**2/(2*sigma**2))

def exponential_pdf(x, a, b):
    return a*np.exp(-b*x)

def long_range_pdf(x, a,b):

    pass
