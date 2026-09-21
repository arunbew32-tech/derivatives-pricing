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


# assumptions the model is making:
# volatility is constant 
# risk free rate is constant 
# no divideds during the options time frame
# all motions are relatively within stochatic processes and no large jumps
# No transaction costs or taxes


# making a put option 
def Black_Scholes_put(S_0, r,T,X,sigma):
    d_1 = ((np.log(S_0/X)) + (r + (sigma**(2))/2)*T)/(sigma*(np.sqrt(T)))
    d_2 = d_1 - sigma*np.sqrt(T)

    P =  X * np.exp(-r * T) * norm.cdf(-d_2)- S_0*norm.cdf(-d_1) 

    return P

put_price = Black_Scholes_put(1.2, 0.02, 0.5, 1.25, 0.25)

#choice between the 2
Option = input("call or put?   ")

if Option == "call":
    print(f"Your call price is{call_price :.3f}")
elif Option == "put":
    print(f"Your put price is {put_price :.3f}")
else:
    print("ERROR")