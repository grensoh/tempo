import numpy as np
from scipy.optimize import minimize_scalar

m_aimant = 0.000636
m_charge = 0.01

print(f"Moment magnétique : {(1.2*(np.pi* 0.003**3)/(4*np.pi* 10**(-7)))} Am²")

Rb = 1.72 * 10**(-8) * (25/(np.pi * 0.00015**2))

print(f"Résistance de la bobine : {Rb}")

I = 0.4

print(f"Courant par défaut dans la bobine : {I} A")

def force(I, texte):
    F = (3 * 0.003 * np.pi * 0.003**2 * 1.2 * 3 *  850 * I * 0.008**2 * 0.01) / (2 * (0.008**2 + 0.01**2)**(5/2))
    if texte == "PNP":
        F += -9.81*(3*m_aimant + m_charge)
    elif texte == "NPN":
        F += 9.81*(3*m_aimant + m_charge)
    
    return f"Force totale {texte}: {F} N"

print(force(I, ""))

while True:
    
    #I = float(input("Courant cible (A) : "))


    R2, R4, R5, R6, R7 = input("R2 R4 R5 R6 R7 : ").split(" ")
    R2, R4, R5, R6, R7 = int(R2), int(R4), int(R5), int(R6), int(R7)
    Ic_PNP = 6/((R2)+R5+R7+Rb)
    Ic_NPN = 6/((R4/2) + R5 +R6 +Rb)
    print(f"PNP = {Ic_PNP} A >>>>> {force(Ic_PNP, "PNP")}")
    print(f"NPN = {Ic_NPN/2} A par NPN >>>>> {force(Ic_NPN, "NPN")}")
    print("--------------------------------------------")

    

        
