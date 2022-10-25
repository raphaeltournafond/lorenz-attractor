# Imports

import math # floor

# plots
import matplotlib.pyplot as plt
import matplotlib.figure as Figure

# GUI
import tkinter as tk
from tkinter import *
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

# Misc
import os


# Constants

sigma       = 10
beta        = 8/3
rho         = 23
pas         = 0.01

iteration = 40
x_init = 0
y_init = 0.1
z_init = 0
i_pas = math.floor(iteration / pas)

# Lorenz

def lorenz(x: float, y: float, z: float):
    return  sigma * (y - x), x * (rho - z) - y, x * y - beta * z
    
# RK4

def derniere_etape_rk4(v: float, k1: float, k2: float, k3: float, k4: float):
    return v + (k1 + 2 * k2 + 2 * k3 + k4) * (pas / 6)
    
def runge_kutta_4(x: float, y: float, z: float):
    "Formule RK4 pour trouver le prochain point"
    k1x, k1y, k1z = lorenz(x, y, z)
    
    P_2 = pas / 2
    
    k2x, k2y, k2z = lorenz(x + k1x * P_2, y + k1y * P_2, z + k1z * P_2)
    
    k3x, k3y, k3z = lorenz(x + k2x * P_2, y + k2y * P_2, z + k2z * P_2)
    
    k4x, k4y, k4z = lorenz(x + pas * k3x, y + pas * k3y, z + pas * k3z)
    
    x1 = derniere_etape_rk4(x, k1x, k2x, k3x, k4x)
    y1 = derniere_etape_rk4(y, k1y, k2y, k3y, k4y)
    z1 = derniere_etape_rk4(z, k1z, k2z, k3z, k4z)
    
    return x1, y1, z1
        
def lorenz_rk4():
    x, y, z = [x_init], [y_init], [z_init]
    for i in range(0, i_pas):
        xi_1, yi_1, zi_1 = runge_kutta_4(x[i], y[i], z[i])
        x.append(xi_1)
        y.append(yi_1)
        z.append(zi_1)
    return x, y, z

# Euler

def lorenz_euler():
    x, y, z = [x_init], [y_init], [z_init]
    for i in range(0, i_pas):
        xi_1, yi_1, zi_1 = lorenz(x[i], y[i], z[i])
        x.append(x[i] + xi_1 * pas)
        y.append(y[i] + yi_1 * pas)
        z.append(z[i] + zi_1 * pas)
    return x, y, z

# Graphes

fig_lorenz: Figure
fig_axes: Figure

def genere_rk4_euler(rk4: tuple, euler: tuple):
    # Unpack
    x_rk4, y_rk4, z_rk4 = rk4
    x_euler, y_euler, z_euler = euler
    
    fig = plt.figure()
    
    ax = fig.gca(projection='3d')
    
    ax.plot(x_rk4, y_rk4, z_rk4, 'r', linewidth = 0.5) # RK4
    ax.plot(x_euler, y_euler, z_euler, 'b', linewidth = 0.5) # Euler
    
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    plt.legend(["RK4","Euler"])
    
    plt.title('Attracteur de Lorenz :x0 = '
              + str(x_init)
              + ', y0 = ' + str(y_init)
              + ', z0 = ' + str(z_init)
              + '\nsigma = ' + str(sigma)
              + ', beta = ' + str(beta)
              + ', rho = ' + str(rho))
    
    return fig
    
def genere_diff_axes(rk4: tuple, euler: tuple):
    # Unpack
    x_rk4, y_rk4, z_rk4 = rk4
    x_euler, y_euler, z_euler = euler
    
    fig, axs = plt.subplots(3)
    
    fig.suptitle('Ecarts par axe RK4 - Euler')
    
    range_pas =  range(i_pas + 1)
    
    axs[0].plot(range_pas, x_rk4, 'r', linewidth=0.5)
    axs[0].plot(range_pas, x_euler, 'b', linewidth=0.5)
    axs[1].plot(range_pas, y_rk4, 'r', linewidth=0.5)
    axs[1].plot(range_pas, y_euler, 'b', linewidth=0.5)
    axs[2].plot(range_pas, z_rk4, 'r', linewidth=0.5)
    axs[2].plot(range_pas, z_euler, 'b', linewidth=0.5)
    
    axs[0].legend(["x rk4", "x euler"])
    axs[1].legend(["y rk4", "y euler"])
    axs[2].legend(["z rk4", "z euler"])
    
    return fig

def genere_graphs(rk4: tuple, euler: tuple):
    global fig_lorenz, fig_axes
    
    fig_lorenz = genere_rk4_euler(rk4, euler)
    fig_axes = genere_diff_axes(rk4, euler)
    
# GUI - main

lorenz_canvas: FigureCanvasTkAgg
axes_canvas: FigureCanvasTkAgg

def rafraichirGraphes(root: Tk, s: Scale, b: Scale, r:Scale, x: Scale, y: Scale, z: Scale, i: Scale, p: Scale):
    global sigma, beta, rho, x_init, y_init, z_init, iteration, pas, i_pas, lorenz_canvas, axes_canvas
    sigma = s.get()
    beta = b.get()
    rho = r.get()
    x_init = x.get()
    y_init = y.get()
    z_init = z.get()
    iteration = i.get()
    pas = p.get()
    i_pas = math.floor(iteration / pas)
    
    rk4 = lorenz_rk4()
    euler = lorenz_euler()
    
    genere_graphs(rk4, euler)
    
    lorenz_canvas.get_tk_widget().destroy()
    lorenz_canvas = FigureCanvasTkAgg(fig_lorenz, root)
    lorenz_canvas.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH)
    lorenz_canvas.draw()
    
    axes_canvas.get_tk_widget().destroy()
    axes_canvas = FigureCanvasTkAgg(fig_axes, root)
    axes_canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand = 1)
    axes_canvas.draw()
    
# Sauvegarde

def creer_dossier_images():
    #python program to check if a directory exists
    path = "./images"
    # Check whether the specified path exists or not
    isExist = os.path.exists(path)
    if not isExist:
        # Create a new directory because it does not exist
        os.makedirs(path)
        
def sauvegarde():
    global fig_lorenz, fig_axes
    
    creer_dossier_images()
    
    fig_lorenz.savefig('./images/attracteur.png')
    fig_axes.savefig('./images/ecarts.png')

def main():
    global lorenz_canvas, axes_canvas, fig_lorenz, fig_axes
    
    root = tk.Tk()
    
    controles = Frame(root)
    
    # Scales
    
    scales = Frame(controles)
    
    x = Scale(scales, label = 'x initial', from_ = -200.0, to = 200.0, resolution = 0.1, length = 600, orient = HORIZONTAL)
    x.set(x_init)
    x.pack()
    
    y = Scale(scales, label = 'y initial', from_ = -200.0, to = 200.0, resolution = 0.1, length = 600, orient = HORIZONTAL)
    y.set(y_init)
    y.pack()
    
    z = Scale(scales, label = 'z initial', from_ = -200.0, to = 200.0, resolution = 0.1, length = 600, orient = HORIZONTAL)
    z.set(z_init)
    z.pack()
    
    i = Scale(scales, label = 'Iterations', from_ = 1, to = 1000, resolution = 1 , length = 600, orient = HORIZONTAL)
    i.set(iteration)
    i.pack()
    
    s = Scale(scales, label = 'Sigma', from_ = 0, to = 50, resolution = 1.0, length = 600, orient = HORIZONTAL)
    s.set(sigma)
    s.pack()
    
    b = Scale(scales, label = 'Beta', from_ = 0, to = 10, resolution = 0.0001, length = 600, orient = HORIZONTAL)
    b.set(beta)
    b.pack()
    
    r = Scale(scales, label = 'Rho', from_ = 0, to = 100, resolution = 1.0, length = 600, orient = HORIZONTAL)
    r.set(rho)
    r.pack()
    
    p = Scale(scales, label = 'Pas', from_ = 0.01, to = 0.0001, resolution =  0.0001, length = 600, orient = HORIZONTAL)
    p.set(pas)
    p.pack()
    
    scales.pack(side = tk.TOP)
    
    # Bouton
    
    refresh = Button(controles, text = 'Rafraichir', bg = 'cyan', command =lambda: rafraichirGraphes(root, s, b, r, x, y, z, i, p))
    refresh.pack()
    
    save = Button(controles, text = 'Sauvegarder', bg = 'magenta', command =lambda: sauvegarde())
    save.pack(side = tk.RIGHT)
    
    
    controles.pack(side = tk.BOTTOM, padx=20, pady=20)
    
    # Graphes
    
    rk4 = lorenz_rk4()
    euler = lorenz_euler()
    
    genere_graphs(rk4, euler)
    
    lorenz_canvas = FigureCanvasTkAgg(fig_lorenz, root)
    lorenz_canvas.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH)
    
    axes_canvas = FigureCanvasTkAgg(fig_axes, root)
    axes_canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand = 1)
    
    # Affichage fenetre
    
    root.state('zoomed') # Plein ecran
    root.mainloop()
    
main()
