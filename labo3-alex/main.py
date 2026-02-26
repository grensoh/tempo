import numpy as np
import matplotlib . pyplot as plt
import pandas as pd


# Fichier
#Saturation: C:/Users/alexa/OneDrive/Documents/UNIF INGE BAC1/Q2/Projet 2/Labo/labo3_saturation_data.csv
#V_in- << V_out: C:/Users/alexa/OneDrive/Documents/UNIF INGE BAC1/Q2/Projet 2/Labo/labo3_gain_data.csv
#expérience1: C:\Users\alexa\OneDrive\Documents\UNIF INGE BAC1\Q2\Projet 2\Labo\labo3_exp1_clean_data.csv

file_location = "labo3_exp1_clean_data.csv"
#Attention à changer les \ en / (du path)

# Lecture du fichier
with open(file_location, 'r') as f:
    lines = f.readlines()

data = pd.read_csv(file_location, sep=";", decimal=",", skiprows=2)

#Prend les données dans le fichier [lignes, colonnes]
#pour saturation: 1227:2249
#pour exp1: 1176:2245
time = data.iloc[1249:2245, 0] 
signal1 = data.iloc[1249:2245, 1]
signal2 = data.iloc[1249:2245, 2]

# Affichage

plt.plot( time , signal1 , label = r"$V_{in-}$" , color = "#744FC6" , linewidth=2) #Courbe de V_in-

plt.plot( time , signal2 , label = r"$V_{out}$" , color = "#4F86C6" , linewidth=2) #Courbe de V_out


plt.xlabel( " Temps [s] " ) 
plt.ylabel( " Tension [V] " ) 
plt.title( " Expérience 2 " ) # Titre
plt.legend(loc="upper right")
plt.xlim(time.min() , time.max()) # Limite de l'axe x
plt.ylim(-1,7)# Limite de l'axe y

plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)

plt.tight_layout()
plt.show()
