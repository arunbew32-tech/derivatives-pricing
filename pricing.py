import numpy as np
from scipy.stats import norm

def d1_d2(S_0, r,T,X,sigma,q):
    d_1 = ((np.log(S_0/X)) + (r - q + (sigma**(2))/2)*T)/(sigma*(np.sqrt(T)))
    d_2 = d_1 - sigma*np.sqrt(T)
    return d_1, d_2


def black_scholes(S_0, r,T,X,sigma,q, option_type):
    d_1, d_2 = d1_d2(S_0, r,T,X,sigma,q)

    if option_type == "call":
        C = S_0*np.exp(-q * T)*norm.cdf(d_1) - X * np.exp(-r * T) * norm.cdf(d_2)
        return C
    elif option_type == "put":
        P = X*np.exp(-r * T) * norm.cdf(-d_2) - S_0*np.exp(-q * T)*norm.cdf(-d_1)
        return P
    else:
        return("Error")

def put_call_parity(S_0, r,T,X,q):
    return S_0*np.exp(-q* T) - X*np.exp(-r*T)
