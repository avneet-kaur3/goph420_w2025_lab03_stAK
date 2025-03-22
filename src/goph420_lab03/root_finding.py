# Rearranging for the root-finding equation and plotting it.
'''
The following section of code rearranges equation (1) provided in the lab manual to a root finding equation that equates to 0.
To initialize the method of finding the roots and to supplement it, the first part plots this function. The frequency can be varied.
'''
import numpy as np
import matplotlib.pyplot as plt

#defining the constants for g(zeta).
freq = 1.0 #hertz
rho_1 = 1800 #kg/m^3
rho_2 = 2500 #kg/m^3
beta_1 = 1900 #m/s
beta_2 = 3200 #m/s
H = 4000 #km

#defining the function g(ζ) that we will be finding the roots for.
def g(zeta):
    y = np.tan(2*np.pi*freq*zeta) - ((rho_2/rho_1)*(np.sqrt((H**2)*((beta_1**-2) - (beta_2**-2)) - zeta**2)/zeta))
    return y
#defining a maximum value for the ζ based on lab methods.
zeta_max = np.sqrt((H**2)*((beta_1**-2) - (beta_2**-2)))
#choosing a very small value of the minimum ζ to avoid division by zero error.
zeta_lower = 1e-3
k = 0

#defining the graph.
plt.figure()
plt.xlabel('ζ')
plt.ylabel('g(ζ)')
plt.grid(True)

#defining a loop to detect asymptotes and modify the graph accordingly.
while zeta_lower <= zeta_max:
    zeta_upper = ((2*k)+1)/(4*freq)
    interval = np.linspace(zeta_lower, zeta_upper, 500)
    y = g(interval)
    
    #setting a value for the asymptotes to take on NaN value after a certain y value has been reached to see the roots easier.
    y[np.abs(y) > 1e2] = np.nan
    #plotting g(ζ) as a function of ζ for the frequency chosen.
    plt.plot(interval, y)

    zeta_lower = zeta_upper
    k += 1

plt.show()

#Iterative Newton-Raphson method for finding the roots of g(ζ).
'''
The following section of code defines a function for the first derivative of g(ζ) and the iterative method for Newton-Raphson root finding. 
It returns the root, number of iterations and an array of relative error given an initial guess of the root.'
'''
#defining the first derivative of the function. the frequency has been kept variable so that it can be easily modified.
def dgdx(zeta):
    y_prime = ((2*np.pi*freq)/(np.cos(2*np.pi*freq*zeta)**2)) + (138125/(456*(zeta**2)*np.sqrt(16575-(5776*(zeta**2)))))
    return y_prime

#defining the newton-raphson method.
def root_newton_raphson(x0, f, dfdx):
    '''
    The function finds roots of a given function f with first derivative dfdx based on the Newton Raphson method.
    Input: initial guess x0 (float), function f, function dfdx
    Output: final estimate of the root (float), number of iterations (int), one-dimensional vector of the relative error at each
            iteration (numpy.ndarray)
    '''
    #keeping a count of the number of iterations.
    iterations = 0
    #defining a list to append the relative errors into.
    rel_error = []
    #basing the while loop off a large number of maximum iterations. most iterations for this problem fall well below half of this number.
    while iterations < 50:
        x1 = x0 - (f(x0)/dfdx(x0))
        error = abs(x1 -x0)
        rel_error.append(error) 
        iterations += 1
        #if the relative error is below a tolerance of 1e-4, where 4 digits have converged, the method stops and identifies a root.
        if error < 1e-4:
            return x1, iterations, np.array(rel_error)
        else:
            x0 = x1

# Identifying the roots based on four frequencies: 0.25, 0.5, 0.75 and 1.0 Hz.
'''
The first two sections of code can be individually run to find all roots based on the graph plotted and our iterative Newton-Raphson method. 
This section of the code identifies all roots (or modes) of ζ to further find the velocities and wavelengths. 
At the very end, the velocities, c_L, and wavelengths, λ_L, are plotted against frequency for each mode.'
'''
#for a frequency f = 0.25 Hz.
#based on the graph plotted when the variable freq is changed to a value of 0.25, there appears to be 1 root.
freq = 0.25
root_1_freq1, iter_1_freq1, error_1_freq1 = root_newton_raphson(0.2, g, dgdx)
print("For frequency of 0.25 Hz, the first root is:", root_1_freq1)

print("\n")

#for a frequency f = 0.5 Hz.
#based on the graph plotted when the variable freq is changed to a value of 0.5, there appear to be 2 roots.
freq = 0.5
root_1_freq2, iter_1_freq2, error_1_freq2 = root_newton_raphson(0.4, g, dgdx)
print("For frequency of 0.5 Hz, the first root is:", root_1_freq2)
root_2_freq2, iter_2_freq2, error_2_freq2 = root_newton_raphson(1.1, g, dgdx)
print("For frequency of 0.5 Hz, the second root is:", root_2_freq2)

print("\n")

#for a frequency f = 0.75 Hz.
#based on the graph plotted when the variable freq is changed to a value of 0.75, there appear to be 3 roots.
freq = 0.75
root_1_freq3, iter_1_freq3, error_1_freq3 = root_newton_raphson(0.3, g, dgdx)
print("For frequency of 0.75 Hz, the first root is:", root_1_freq3)
root_2_freq3, iter_2_freq3, error_2_freq3 = root_newton_raphson(0.6, g, dgdx)
print("For frequency of 0.75 Hz, the second root is:", root_2_freq3)
root_3_freq3, iter_3_freq3, error_3_freq3 = root_newton_raphson(1.2, g, dgdx)
print("For frequency of 0.75 Hz, the third root is:", root_3_freq3)

print("\n")

#for a frequency f = 1.0 Hz.
#based on the graph plotted when the variable freq is changed to a value of 1.0, there appear to be 4 roots.
freq = 1.0
root_1_freq4, iter_1_freq4, error_1_freq4 = root_newton_raphson(0.24, g, dgdx)
print("For frequency of 1.0 Hz, the first root is:", root_1_freq4)
root_2_freq4, iter_2_freq4, error_2_freq4 = root_newton_raphson(0.7, g, dgdx)
print("For frequency of 1.0 Hz, the second root is:", root_2_freq4)
root_3_freq4, iter_3_freq4, error_3_freq4 = root_newton_raphson(1.0, g, dgdx)
print("For frequency of 1.0 Hz, the third root is:", root_3_freq4)
root_4_freq4, iter_4_freq4, error_4_freq4 = root_newton_raphson(1.5, g, dgdx)
print("For frequency of 1.0 Hz, the fourth root is:", root_4_freq4)

#defining the modes and the frequencies they correspond to.
mode_1 = np.array([0.7743642822513618, 0.4392224915420927, 0.30548243171235145, 0.23406653286159174])
mode_1_freq = np.array([0.25, 0.50, 0.75, 1.00])
mode_2 = np.array([1.2795167697299594, 0.9088987590779034, 0.6997603830453311])
mode_2_freq = np.array([0.50, 0.75, 1.00])
mode_3 = np.array([1.4746031844344674, 1.1558320051102862])
mode_3_freq = np.array([0.75, 1.00])
mode_4 = 1.5788173798566245
mode_4_freq = 1.00

#calculating c_L and lambda_L from equations 2 and 3 for each mode.
c_L_1 = np.sqrt(1/((1/beta_1**2)-((mode_1/H)**2)))
lambda_L_1 = c_L_1/mode_1_freq

c_L_2 = np.sqrt(1/((1/beta_1**2)-((mode_2/H)**2)))
lambda_L_2 = c_L_2/mode_2_freq

c_L_3 = np.sqrt(1/((1/beta_1**2)-((mode_3/H)**2)))
lambda_L_3 = c_L_3/mode_3_freq

c_L_4 = np.sqrt(1/((1/beta_1**2)-((mode_4/H)**2)))
lambda_L_4 = c_L_4/mode_4_freq

#plotting c_L as a function of frequency for each mode.
plt.figure(figsize=(10, 5))
plt.plot(mode_1_freq, c_L_1, label="Mode 1")
plt.plot(mode_2_freq, c_L_2, label="Mode 2")
plt.plot(mode_3_freq, c_L_3, label="Mode 3")
plt.plot(mode_4_freq, c_L_4, label="Mode 4", marker='o')
plt.xlabel("Frequency (Hz)")
plt.ylabel("Velocity, c_L (m/s)")
plt.title("c_L vs Frequency")
plt.legend()
plt.grid(True)
plt.show()

#plotting lambda_L as a function of frequency for each mode.
plt.figure(figsize=(10, 5))
plt.plot(mode_1_freq, lambda_L_1, label="Mode 1")
plt.plot(mode_2_freq, lambda_L_2, label="Mode 2")
plt.plot(mode_3_freq, lambda_L_3, label="Mode 3")
plt.plot(mode_4_freq, lambda_L_4, label="Mode 4", marker='o')
plt.xlabel("Frequency (Hz)")
plt.ylabel("Wavelength, λ_L (m)")
plt.title("λ_L vs Frequency")
plt.legend()
plt.grid(True)
plt.show()