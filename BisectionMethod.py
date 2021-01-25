import numpy as np 

maxError = .001  #Sets the maximum error allowed in the root
maxN = 100  #Sets the maximum number of iterations
a = 0.01    #Sets the left endpoint of the interval
b = 0.99    #Sets the right endpoint of the interval

def func(x):    #Defines the function that we seek a root for
    return (x**2*(3-x))/(1-x)**3-62.4   

if (func(a)*func(b) > 0):
    print("Inconclusive, there may or may not be root in the interval")
elif (func(a) == 0):
    print("The root is at x = ", a)
elif (func(b) == 0):
    print("The root is at x = ", b)
else:
    for n in range(maxN):
        m = (a+b)/2
        if (abs(func(m)) <= abs(maxError)):
             print("The root is approximately x = ", m)
             break
        else:
            if (np.sign(func(a)) == np.sign(func(m))):
                a = m
            elif (np.sign(func(b)) == np.sign(func(m))):
                b = m
        
        
