import numpy as np
from scipy.stats import norm
# assumptions the model is making:
# volatility is constant 
# risk free rate is constant 
# no divideds during the options time frame
# all motions are relatively within normal processes and no large jumps
# No transaction costs or taxes

def black_scholes(S_0, r,T,X,sigma, option_type):
    d_1 = ((np.log(S_0/X)) + (r + (sigma**(2))/2)*T)/(sigma*(np.sqrt(T)))
    d_2 = d_1 - sigma*np.sqrt(T)
    

    if option_type == "call":
        C = S_0*norm.cdf(d_1) - X * np.exp(-r * T) * norm.cdf(d_2)
        return C
    elif option_type == "put":
        P = X*np.exp(-r * T) * norm.cdf(-d_2) - S_0*norm.cdf(-d_1)
        return P
    else:
        return("Error")
    
def put_call_parity(S_0, r,T,X):
    return S_0 - X*np.exp(-r*T)

call_price = black_scholes(1.2, 0.02, 0.5, 1.25, 0.25, "call")
put_price = black_scholes(1.2,0.02,0.5,1.25,0.25,"put")

discrepancy = (call_price - put_price) - put_call_parity(1.2, 0.02, 0.5, 1.25)

print(f"Discrepancy:  {discrepancy:.5f}  ")
