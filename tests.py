from pricing import black_scholes
from greeks import delta, gamma, theta, vega, rho

# implementing finite difference approximation for the greeks

def finite_difference_delta_test(S_0, r,T,X,sigma,q, option_type):
    approx_delta = (black_scholes(S_0+ S_0*0.001, r,T,X,sigma,q, option_type) - black_scholes(S_0 - S_0*0.001, r,T,X,sigma,q, option_type))/(2*S_0*0.001)
    assert abs(delta(S_0, r,T,X,sigma,q, option_type)- approx_delta) < 1e-4 


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