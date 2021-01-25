import numpy as np

maxError = .001 #Sets the maximum error allowed in the root
maxN = 100  #Sets the maximum number of iterations allowed
initGuess = .5  #Sets the inital guess
epsilon = 10**(-14) #Prevents division by too small a number

def func(x): #Function to find root of
    return np.cos(x)-x**3

def funcPrime(x): #1st derivative of func(x)
    return -np.sin(x)-3*x**2

m = initGuess
for i in range(maxN):
    if (abs(func(m)) <= maxError):
        print("The root is at x = ", m)
        break
    else:
        if (epsilon <= abs(funcPrime(m))):
            m = m - func(m)/funcPrime(m)
        else:
            print("Denominator too small")
            break
