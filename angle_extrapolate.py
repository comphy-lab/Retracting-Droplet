from scipy.optimize import fsolve
import numpy as np

h = 6e-5
V = 5e-10
theta = 23.5*np.pi/180

# Define the equation
def equation(thetaS):
    return np.pi*h**3*(np.sin(thetaS))**4*(2+np.cos(thetaS))/(3*V*(1+np.cos(thetaS))**2) - (abs(np.cos(theta)-np.cos(thetaS)))**3

# Initial guess for the root
initial_guess = 35*np.pi/180

# Use fsolve to find the root
root = fsolve(equation, initial_guess)
root = np.degrees(root)

print(f"{root}")