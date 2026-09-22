import numpy as np
from scipy.stats import norm

# assumptions the model is making:
# volatility is constant 
# risk free rate is constant 
# underlying follows GBM
# no jumps in underlying price
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
    return (S_0*np.exp(-q*T)*norm.pdf(d_1)*np.sqrt(T))*0.01

def rho(S_0, r,T,X,sigma,q, option_type): 
    d_1, d_2 = d1_d2(S_0, r,T,X,sigma,q)
    if option_type == "call":
        return X*T*np.exp(-r * T)*norm.cdf(d_2)*0.01
    elif option_type == "put":
        return -X*T*np.exp(-r * T)*norm.cdf(-d_2)*0.01
    else:
        return("Error")

# implementing finite difference approximation for the greeks

def finite_difference_delta_test(S_0, r,T,X,sigma,q, option_type):
    approx_delta = (black_scholes(S_0+ S_0*0.001, r,T,X,sigma,q, option_type) - black_scholes(S_0 - S_0*0.001, r,T,X,sigma,q, option_type))/(2*S_0*0.001)
    assert abs(delta(S_0, r,T,X,sigma,q, option_type)- approx_delta) < 1e-4 

print(finite_difference_delta_test(1.2, 0.02, 0.5, 1.25, 0.25, 0,"put"))

def finite_difference_gamma_test(S_0, r,T,X,sigma,q, option_type):
    approx_gamma = (black_scholes(S_0+ S_0*0.0001, r,T,X,sigma,q, option_type)- 2*black_scholes(S_0, r,T,X,sigma,q, option_type) + black_scholes(S_0 - S_0*0.0001, r,T,X,sigma,q, option_type))/((S_0*0.0001)**2)
    assert abs(gamma(S_0, r,T,X,sigma,q)- approx_gamma) < 1e-4 

def finite_difference_theta_test(S_0, r,T,X,sigma,q, option_type):
    approx_theta = -(black_scholes(S_0, r,T+1/365,X,sigma,q, option_type) - black_scholes(S_0, r,T-1/365,X,sigma,q, option_type))/2
    assert abs(theta(S_0, r,T,X,sigma,q, option_type)- approx_theta) < 1e-4 

def finite_difference_vega_test(S_0, r,T,X,sigma,q, option_type):
    approx_vega = (black_scholes(S_0, r,T,X,sigma+0.01,q, option_type) - black_scholes(S_0, r,T,X,sigma-0.01,q, option_type))/2
    assert abs(vega(S_0, r,T,X,sigma,q)- approx_vega) < 1e-4 

def finite_difference_rho_test(S_0, r,T,X,sigma,q, option_type):
    approx_rho = 0.01*(black_scholes(S_0, r+0.001,T,X,sigma,q, option_type) - black_scholes(S_0, r-0.001,T,X,sigma,q, option_type))/(2*0.001)
    assert abs(rho(S_0, r,T,X,sigma,q,option_type)- approx_rho) < 1e-4 


print(f"price : {black_scholes(1.2, 0.02, 0.5, 1.25, 0.25, 0,"put"):.4f}")
print(f"delta:  {delta(1.2, 0.02, 0.5, 1.25, 0.25, 0,"put"):.4f}    approx = {finite_difference_delta_test(1.2, 0.02, 0.5, 1.25, 0.25, 0,"put")}")
print(f"gamma:  {gamma(1.2, 0.02, 0.5, 1.25, 0.25, 0):.4f}          approx = {finite_difference_gamma_test(1.2, 0.02, 0.5, 1.25, 0.25, 0,"put")}")
print(f"theta:  {theta(1.2, 0.02, 0.5, 1.25, 0.25, 0,"put"):.4f}    approx = {finite_difference_theta_test(1.2, 0.02, 0.5, 1.25, 0.25, 0,"put")}")
print(f"vega:   {vega(1.2, 0.02, 0.5, 1.25, 0.25, 0):.4f}           approx = {finite_difference_vega_test(1.2, 0.02, 0.5, 1.25, 0.25, 0,"put")}")
print(f"rho:    {rho(1.2, 0.02, 0.5, 1.25, 0.25, 0,"put"):.4f}      approx = {finite_difference_rho_test(1.2, 0.02, 0.5, 1.25, 0.25, 0,"put")}")
