import numpy as np
from scipy.stats import norm

# assumptions the model is making:
# volatility is constant 
# risk free rate is constant 
# no divideds during the options time frame
# all motions are relatively within normal processes and no large jumps
# No transaction costs or taxes

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

call_price = black_scholes(1.2, 0.02, 0.5, 1.25, 0.25, 0,"call")
put_price = black_scholes(1.2,0.02,0.5,1.25,0.25,0,"put")

discrepancy = (call_price - put_price) - put_call_parity(1.2, 0.02, 0.5, 1.25,0)

print(f"Discrepancy:  {discrepancy:.5f}  ")


# greeks introduction

# delta - how much the option will increase if the stock increased by one unit 
# (if delta = 0.2 then a £1 stock incerase see's a £0.2 option increase) 
# gamma - the rate of change of the delta 
# theta - the negative change in the value of the option over the time to expiry 
# vega - the change of the option price with respect to the volatiltiy 
# rho -  change of the options value respect to the risk free rate 

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
        return -np.exp(-q *T)*((S_0 * norm.pdf(d_1)* sigma)/(2* np.sqrt(T))) - r*X*np.exp(-r *T)*norm.cdf(d_2) + q*S_0*np.exp(-q * T)*norm.cdf(d_1)
    elif option_type == "put":
        return -np.exp(-q *T)*((S_0 * norm.pdf(d_1)* sigma)/(2* np.sqrt(T))) + r*X*np.exp(-r *T)*norm.cdf(-d_2) - q*S_0*np.exp(-q * T)*norm.cdf(-d_1)
    else:
        return("Error")

def vega(S_0, r,T,X,sigma,q):
    d_1, d_2 = d1_d2(S_0, r,T,X,sigma,q)
    return S_0*np.exp(-q*T)*norm.pdf(d_1)*np.sqrt(T)

def rho(S_0, r,T,X,sigma,q, option_type):
    d_1, d_2 = d1_d2(S_0, r,T,X,sigma,q)
    if option_type == "call":
        return X*T*np.exp(-r * T)*norm.cdf(d_2)
    elif option_type == "put":
        return -X*T*np.exp(-r * T)*norm.cdf(-d_2)
    else:
        return("Error")

print(f"call price : {black_scholes(1.2, 0.02, 0.5, 1.25, 0.25, 0,"call"):.4f}")
print(f"delta:  {delta(1.2, 0.02, 0.5, 1.25, 0.25, 0,"call"):.4f}")
print(f"gamma:  {gamma(1.2, 0.02, 0.5, 1.25, 0.25, 0):.4f}")
print(f"theta:  {theta(1.2, 0.02, 0.5, 1.25, 0.25, 0,"call"):.4f}")
print(f"vega:  {vega(1.2, 0.02, 0.5, 1.25, 0.25, 0):.4f}")
print(f"rho:    {rho(1.2, 0.02, 0.5, 1.25, 0.25, 0,"call"):.4f}")
