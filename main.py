import numpy as np
from scipy.stats import norm
# values


#defining a function to run the whole operation in one go
def Black_Scholes_call(S_0, r,T,X,sigma):
    d_1 = ((np.log(S_0/X)) + (r + (sigma**(2))/2)*T)/(sigma*(np.sqrt(T)))
    d_2 = d_1 - sigma*np.sqrt(T)
    
    C = S_0*norm.cdf(d_1) - X * np.exp(-r * T) * norm.cdf(d_2)

    return C

call_price = Black_Scholes_call(1.2, 0.02, 0.5, 1.25, 0.25)

print(f"The price of your call option is ${call_price:.3f}")