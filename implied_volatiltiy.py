from greeks import vega
from pricing import black_scholes

def implied_newton_raphson(S_0, r,T,X,sigma_n,q ,option_type ,C_price):
    i =0 

    tolerance = 1e-6
    error = float("inf")
    while error > tolerance and i < 1000:
        bs_price = black_scholes(S_0, r,T,X,sigma_n,q, option_type)
        sigma_n = sigma_n - (bs_price-C_price)/(vega(S_0, r,T,X,sigma_n,q)) 
        error = abs(bs_price - C_price)
        i += 1

        if error <= tolerance:
            return sigma_n
        
    raise ValueError("Failed to converge")


price = black_scholes(1.2, 0.02, 0.5, 1.25, 0.35, 0, "put")
print(implied_newton_raphson(S_0 = 1.2, r = 0.02,T =  0.5,X = 1.25,sigma_n= 0.05,q= 0,option_type="put",C_price = price ))
# is not directly observable but can be calculated using the black-scholes model
# using an iterative model to estimate the implied volatility, newton raphson
