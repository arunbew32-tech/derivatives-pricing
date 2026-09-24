import numpy as np
from scipy.stats import norm
from pricing import d1_d2

# greeks introduction

# delta - change in option price per £1 underlying move
# gamma - change in delta per £1 underlying move
# theta - option price change per day
# vega - option price change per 1% volatiltiy change
# rho -  option price change per 1% rate change 

# changed here to the merton model for future reference to make things more general

def delta(S_0, r,T,X,sigma,q, option_type):
    d_1, d_2 = d1_d2(S_0, r,T,X,sigma,q)
    if option_type == "call":
        return np.exp(-q * T)*norm.cdf(d_1)
    elif option_type == "put":
        return -np.exp(-q * T)*norm.cdf(-d_1)
    else:
        return("Error")

def gamma(S_0, r,T,X,sigma,q):
    d_1, d_2 = d1_d2(S_0, r,T,X,sigma,q)
    return np.exp(-q* T)*((norm.pdf(d_1))/(S_0 * sigma * np.sqrt(T)))

def theta(S_0, r,T,X,sigma,q, option_type): 
    d_1, d_2 = d1_d2(S_0, r,T,X,sigma,q)

    if option_type == "call":
        return (-np.exp(-q *T)*((S_0 * norm.pdf(d_1)* sigma)/(2* np.sqrt(T))) - r*X*np.exp(-r *T)*norm.cdf(d_2) + q*S_0*np.exp(-q * T)*norm.cdf(d_1))*(1/365)
    elif option_type == "put":
        return (-np.exp(-q *T)*((S_0 * norm.pdf(d_1)* sigma)/(2* np.sqrt(T))) + r*X*np.exp(-r *T)*norm.cdf(-d_2) - q*S_0*np.exp(-q * T)*norm.cdf(-d_1))*(1/365)
    else:
        return("Error")

def vega(S_0, r,T,X,sigma,q): 
    d_1, d_2 = d1_d2(S_0, r,T,X,sigma,q)
    return (S_0*np.exp(-q*T)*norm.pdf(d_1)*np.sqrt(T))

def rho(S_0, r,T,X,sigma,q, option_type): 
    d_1, d_2 = d1_d2(S_0, r,T,X,sigma,q)
    if option_type == "call":
        return X*T*np.exp(-r * T)*norm.cdf(d_2)
    elif option_type == "put":
        return -X*T*np.exp(-r * T)*norm.cdf(-d_2)
    else:
        return("Error")
