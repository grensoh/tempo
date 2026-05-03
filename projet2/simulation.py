# -------------------------------------------------------------------------
#
# PYTHON for DUMMIES 21-22
# Problème 7
#
# Script de test
#  Vincent Legat
#
# Largement inspiré du programme de Nicolas Roisin :-)
# Ou les méthodes numériques pour obtenir la solution du projet P2 !
#
# -------------------------------------------------------------------------

from numpy import *
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# ------------------------------------------------------------------------------------
#
# Intégration de la force de lorentz pour un courant unitaire dans la bobine
# (Modifié pour prendre en compte la géométrie/hauteur de la bobine)
#
def lorentzComputeForce(Xmagnet,Ymagnet,Zmagnet,Rcoil,Hcoil,triangles,Zshift,mu0,mu,nSpires) :
  
  nElem = len(triangles)
  m     = len(Zshift)
  
  surface = zeros(nElem)
  for iElem in range(nElem) :
      x = Xmagnet[triangles[iElem,:]]
      y = Ymagnet[triangles[iElem,:]]    
      surface[iElem] = ((x[1]-x[0])*(y[2]-y[0]) - (y[1]-y[0])*(x[2]-x[0]))/2.0
  surfaceMagnet = sum(surface)
 
  F = zeros(m)
  # Répartition des spires sur la hauteur de la bobine
  Z_spires = linspace(-Hcoil/2.0, Hcoil/2.0, nSpires)

  for i in range(m):  
      for jElem in range(nElem) :
          Xp = Rcoil - mean(Xmagnet[triangles[jElem,:]]) 
          Yp =       - mean(Ymagnet[triangles[jElem,:]])
          
          # On somme la contribution de chaque spire de la bobine
          for z_spire in Z_spires:
              Zp = z_spire - Zmagnet - Zshift[i]
              r  = sqrt(Xp*Xp + Yp*Yp + Zp*Zp)
              coeff = -(mu0*mu) / (4*pi*r**5)
              F[i] += coeff * (3*Zp*Xp) * surface[jElem] / surfaceMagnet
              
      # On multiplie par la circonférence de la spire (nSpires est géré par la boucle ci-dessus)
      F[i] *= (2 * pi * Rcoil)
      
  return F
  
# ------------------------------------------------------------------------------------
#
# Intégration du mouvement de l'aimant par la méthode classique de Runge-Kutta d'ordre 4
#
def lorentzRungeKutta(Xstart,Xend,Ustart,n,f):
  X = linspace(Xstart,Xend,n+1)
  U = zeros((n+1,2)); U[0,:] = Ustart
  h = (Xend - Xstart)/n
  for i in range(n):
    t = X[i]

    K1 = f(t, U[i,:]       )
    K2 = f(t + h/2, U[i,:]+K1*h/2)
    K3 = f(t + h/2, U[i,:]+K2*h/2)
    K4 = f(t + h, U[i,:]+K3*h  )
    U[i+1,:] = U[i,:] + h*(K1+2*K2+2*K3+K4)/6

    if (U[i+1,0] < -4) : # Ajusté pour descendre plus bas si besoin
        U[i+1:,0] = -4;
        U[i+1:,1] =  0; break
        
  return X,U
 
# ------------------------------------------------------------------------------------
#
# Interpolation linéaire par morceaux 
#
def lorentzInterpolate(x, X, U):
    x = clip(x, X[0], X[-1])   # 🔥 correction clé
    i = max(insert(where(X[:-1] <= x)[0],0,0))     
    return (U[i]   * (x-X[i+1]) / (X[i]-X[i+1]) + 
            U[i+1] * (x - X[i]) / (X[i+1]-X[i]))
 

# ------------------------------------------------------------------------------------ 
#
# Script de test - Paramètres matériels
#
# ------------------------------------------------------------------------------------

def main() :
  mu0     = 4*pi*10**(-7)     # permeabilité du vide en [H/cm] 
  mMagnet = 0.000636 + 0.002            # masse de l'aimant [kg]
  
  # --- NOUVELLES DIMENSIONS ---
  Rmagnet = 0.3             # rayon de l'aimant [cm]
  Hmagnet = 0.3             # épaisseur de l'aimant [cm]
  Rcoil   = 0.8              # rayon de la bobine [cm]
  Hcoil   = 1.8              # hauteur de la bobine [cm]
  # ----------------------------

  Zmagnet = 0.0              # position verticale de l'aimant en [cm]
  Br      = 1.2              # magnetisation residuelle du NdFeB en [T]
  mu      = Rmagnet**2*Hmagnet*pi*Br / mu0    
  nSpires = 459              # nombre de spires

# ------------------------------------------------------------------------------------
#
# -1- Construction d'un maillage de triangles
#
# ------------------------------------------------------------------------------------
  from scipy.spatial import Delaunay
 
  nR      = 6
  nTheta  = 5
  nNode   = 1 + sum(arange(1,nR))*nTheta
  R       = zeros(nNode)
  Theta   = zeros(nNode)
 
  index = 1; dR = 1.0/(nR-1)
  for i in range(1,nR):
      dTheta = 2*pi/(i*nTheta)
      for j in range(0,i*nTheta):
          R[index]     = i*dR
          Theta[index] = j*dTheta; index += 1
 
  X       = R*cos(Theta)
  Y       = R*sin(Theta)
 
  triangles = Delaunay(stack((X,Y),1)).simplices
  nElem = len(triangles)
  
# ------------------------------------------------------------------------------------
#
# -2- Calcul de la force de lorentz pour diverses hauteurs
#
# ------------------------------------------------------------------------------------
  m       = 61
  Zstart  = -4                        # Elargi pour voir toute la chute [cm]
  Zstop   =  4                        # [cm]
  Zshift  = linspace(Zstart,Zstop,m)
  Tstart  = 0                         # [s]
  Tstop   = 0.5                       # [s]
  T,delta = linspace(Tstart,Tstop,m,retstep=True)
 
  Xmagnet = Rmagnet*R*cos(Theta)
  Ymagnet = Rmagnet*R*sin(Theta)     
  
  # Ajout de Hcoil dans les arguments
  Florentz = lorentzComputeForce(Xmagnet,Ymagnet,Zmagnet,Rcoil,Hcoil,triangles,Zshift,mu0,mu,nSpires) 
  
# ------------------------------------------------------------------------------------
#
# -3- Calcul du mouvement de l'aimant
#
# ------------------------------------------------------------------------------------

  f_Hz = 50
  omega = 2*pi*f_Hz
  amplitude = 0.4


  def f(t, u):
    z = u[0]
    v = u[1]

    I_t = amplitude * sin(omega * t)

    Fcoil = lorentzInterpolate(z, Zshift, Florentz)

    dzdt = v

    dvdt = (1/mMagnet) * (I_t * Fcoil - 9.81 * mMagnet)
    return array([dzdt, dvdt])
 
  Tstart = 0
  Tend   = 5
  n      = 1000   # 👈 IMPORTANT : augmenter la résolution !

  Ustart = [2.0, 0.0]
  T, Uoscillation = lorentzRungeKutta(Tstart, Tend, Ustart, n, f)

# ------------------------------------------------------------------------------------
#
# -4- Quelques jolis plots et animation
#
# ------------------------------------------------------------------------------------
  plt.rcParams['toolbar'] = 'None'
  
  # Calcul des positions des dipoles pour l'affichage
  Xdipole = mean(Xmagnet[triangles[:,:]],axis=1)  
  Ydipole = mean(Ymagnet[triangles[:,:]],axis=1)  

  def frame(i):
    plt.clf()

    n_grid = 80 
    X_grid, Z_grid = meshgrid(linspace(-5,5,n_grid),linspace(-5,5,n_grid))
    Y_grid = zeros_like(X_grid)
    Bx = zeros(shape(X_grid))
    Bz = zeros(shape(X_grid))

    for iElem in range(nElem):
      Xp = X_grid - Xdipole[iElem] 
      Yp = Y_grid - Ydipole[iElem]
      Zp = Z_grid - Zmagnet - Zshift[i]
      r     = sqrt(Xp*Xp + Yp*Yp + Zp*Zp)
      coeff = -(mu0*mu) / (4*pi*r**5)
      Bx   += coeff * (3*Zp*Xp)
      Bz   += coeff * (3*Zp*Zp - r*r)
      
    # Affichage du champ magnétique (sans masque, les lignes traversent tout)
    plt.streamplot(X_grid, Z_grid, Bx, Bz, density=1.4, linewidth=None, color='blue')

    # Dessin de l'aimant rectangulaire
    x = array([-Rmagnet,Rmagnet,Rmagnet,-Rmagnet,-Rmagnet]) 
    y = array([0,0,Hmagnet,Hmagnet,0]) + Zmagnet - Hmagnet/2.0 + Zshift[i]
    plt.fill(x,y,facecolor='blue',alpha=1)

    # Dessin des spires de la bobine (Style "Mutual Induction")
    spires_visuelles = 20 
    z_positions = linspace(-Hcoil/2, Hcoil/2, spires_visuelles)
    
    for z_spire in z_positions:
        plt.plot([-Rcoil, Rcoil], [z_spire, z_spire], "o-r", linewidth=2, markersize=3)
  
    plt.xlim((-5,5)); plt.ylim((-5,5))
    plt.title('Electromagnetic Field')
 
# ------------------------------------------------------------------------------------
  
  fig=plt.figure("Maillage de l'aimant")
  plt.plot(Xmagnet,Ymagnet,'or')
  plt.triplot(Xmagnet,Ymagnet,triangles,'-k')
  plt.plot(Xdipole,Ydipole,'ob')  
  plt.axis("equal")
  plt.axis("off")
 
  plt.figure("Force de lorentz")
  plt.plot(Zshift,Florentz,'-r')
  plt.plot(1.0,lorentzInterpolate(1.0,Zshift,Florentz),'or')
  plt.text(-4.0,100,"Force de lorentz [N]",color='red',fontsize=12)
  plt.text(1.5,-100,"Distance [cm]",color='black',fontsize=12)
 
  fig=plt.figure("Mouvement de l'aimant : z(t) [cm]")
  plt.plot(T,Uoscillation[:,0],'b', label="I = +-0.4 [A]")
  plt.legend()
  
  
  movie = animation.FuncAnimation(plt.figure("Claude's project",figsize=(10,10)),frame,m,interval=20,repeat=False)
  plt.show()

main()
