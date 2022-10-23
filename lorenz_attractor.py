import math
import matplotlib.pyplot as plt

SIGMA = 10
BETA = 8/3
RHO = 28
ITERATION = 40
PAS = 0.01
I_PAS = math.floor(ITERATION / PAS)
INIT_X = 1.0
INIT_Y = 1.0
INIT_Z = 1.0

# Lorenz

def lorenz(x: float, y: float, z: float):
    return  SIGMA * (y - x), x * (RHO - z) - y, x * y - BETA * z
    
# RK4

def derniere_etape_rk4(v: float, k1: float, k2: float, k3: float, k4: float):
    return v + (k1 + 2 * k2 + 2 * k3 + k4) * (PAS / 6)
    
def runge_kutta_4(x: float, y: float, z: float):
    # for _ in range(0, I_PAS):
    "Formule RK4 pour trouver le prochain point"
    k1x, k1y, k1z = lorenz(x, y, z)
    
    P_2 = PAS / 2
    
    k2x, k2y, k2z = lorenz(x + k1x * P_2, y + k1y * P_2, z + k1z * P_2)
    
    k3x, k3y, k3z = lorenz(x + k2x * P_2, y + k2y * P_2, z + k2z * P_2)
    
    k4x, k4y, k4z = lorenz(x + PAS * k3x, y + PAS * k3y, z + PAS * k3z)
    
    x1 = derniere_etape_rk4(x, k1x, k2x, k3x, k4x)
    y1 = derniere_etape_rk4(y, k1y, k2y, k3y, k4y)
    z1 = derniere_etape_rk4(z, k1z, k2z, k3z, k4z)
    
    return x1, y1, z1
        
def lorenz_rk4(x0: float, y0: float, z0: float):
    x, y, z = [x0], [y0], [z0]
    for i in range(0, I_PAS):
        xi_1, yi_1, zi_1 = runge_kutta_4(x[i], y[i], z[i])
        x.append(xi_1)
        y.append(yi_1)
        z.append(zi_1)
    return x, y, z

# Euler

def lorenz_euler(x0: float, y0: float, z0: float):
    x, y, z = [x0], [y0], [z0]
    for i in range(0, I_PAS):
        xi_1, yi_1, zi_1 = lorenz(x[i], y[i], z[i])
        x.append(x[i] + xi_1 * PAS)
        y.append(y[i] + yi_1 * PAS)
        z.append(z[i] + zi_1 * PAS)
    return x, y, z

# Graph

def affiche_graphs(rk4, euler):
    # Unpack
    x_rk4, y_rk4, z_rk4 = rk4
    x_euler, y_euler, z_euler = euler
    
    fig=plt.figure()
    
    ax=fig.gca(projection='3d')
    
    ax.plot(x_rk4, y_rk4, z_rk4, 'r', linewidth=0.2) # RK4
    ax.plot(x_euler, y_euler, z_euler, 'b', linewidth=0.2) # Euler
    
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    
    plt.title("Modélisation de l'équation de Lorenz")
    
    plt.draw()
    plt.show()
    
# Main

rk4 = lorenz_rk4(INIT_X, INIT_Y, INIT_Z)

euler = lorenz_euler(INIT_X, INIT_Y, INIT_Z)

affiche_graphs(rk4, euler)