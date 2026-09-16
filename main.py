import numpy as np
import math
from scipy.stats import norm
# values

S_0 = 1.2 # stock price
r = 0.02 # risk free rate
T = 0.5 # time to maturity (y)
X = 1.25 # strike price
v = 0.25 # volatility

d_1 = ((np.log(S_0/X)) + (r + (v**(2))/2)*T)/(v*(np.sqrt(T)))

d_2 = d_1 - v*np.sqrt(T)

C = S_0*norm.cdf(d_1) - X * np.exp(-r * T) * norm.cdf(d_2)

print(f"The price of your call option is ${C:.3f}")