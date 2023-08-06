import numpy as np
import matplotlib.pyplot as plt

n = 101 # Number of grid points in x-direction
m = 101 # Number of grid points in y-direction
T = np.zeros((n, m)) # 2D array to store temperature at each grid point and time step
T_last = np.zeros((n, m))
T_center = np.zeros((n, m))
T_right = np.zeros((n, m))

x = [0.0] * n
y = [0.0] * m
     
a = np.zeros(n)
b = np.zeros(n)
c = np.zeros(n)
d = np.zeros(n)

time_values = []
temperature_last = []
temperature_center = []
temperature_right = []


Temp_x = np.zeros((n, m)) # Temperatures at next time step in x direction
Temp_y = np.zeros((n, m)) # Temperatures at next time step in y direction

# meter
l = 0.04
dx = l/(n - 1.0) # Grid spacing in x-direction
dy = l/(m - 1.0) # Grid spacing in y-direction
k = 50.0 # Thermal conductivity of material ( W/(m*K))
h = 400.0 # Convective heat transfer coefficient
rho = 7850 # kg/m^3
Cp = 500 # Cp is specific heat capacity(in (j/kg-K)
alpha = k / (rho*Cp) # Thermal diffusivity (in W/m^2-K)
dt = 0.50 # Time step
r = dt * alpha/(2*(dx*dx))# Ratio used in tridiagonal system
t_inf = 25 # Ambient temperature

for i in range(n):
    x[i]=dx*i
    
for j in range(m):
    y[i]=dy*j
# set initial temperature
for j in range(m):
    for i in range(n):
        T[i][j] = 100
        
# initialize t0 as 0
t = 0
for p in range(120):   # repeat the following 600 times
       
    # x-Sweep
    for j in range(1,m-1):
        for i in range(1, n - 1):
            a[i] = -r
            b[i] = 1 + 2*r
            c[i] = -r
            d[i] = r * (T[i][j + 1] + T[i][j - 1]) + (1 - 2*r) * T[i][j]
            
        a[0] = 0.0
        b[0] = -1.0
        c[0] = 1.0
        d[0] = 0.0
        
        a[n - 1] = -k
        b[n - 1] = k + h*dx
        c[n - 1] = 0.0
        d[n - 1] = h*dx*t_inf
        
        for i in range(1, n):
            
            b[i] = b[i] - c[i-1] * a[i] / b[i-1]
            d[i] = d[i] - d[i-1] * a[i] / b[i-1]
       
        Temp_x[n - 1][j] = d[n - 1] / b[n - 1]
       
        for i in range(n - 2, -1, -1):
            Temp_x[i][j] = (d[i] - c[i]*Temp_x[i + 1][j])/b[i]
        

    for i in range(n):
        Temp_x[i][0] = Temp_x[i][1]  # Boundary condition at the bottom
        Temp_x[i][m - 1] = (k * Temp_x[i][m - 2] + h * dy * t_inf) / (k + h * dy)  # Boundary condition at the top
    
    for i in range(n):
        for j in range(m):
            T[i][j]=Temp_x[i][j]
    
    # y_sweep
    for i in range(1, n-1):
        for j in range(1, m - 1):
            a[j] = -r
            b[j] = 1 + 2*r
            c[j] = -r 
            d[j] = r * (T[i+1][j] + T[i-1][j]) + (1 - 2 * r) * T[i][j]
         
        # Boundary conditions at Base
        a[0] = 0.0
        b[0] = -1.0
        c[0] = 1.0
        d[0] = 0.0
                
        a[m - 1] = -k
        b[m - 1] = k + h * dx
        c[m - 1] = 0.0
        d[m - 1] = h * dx * t_inf
        

        for j in range(1, m):
            b[j] = b[j] - c[j - 1] * a[j] / b[j - 1]
            d[j] = d[j] - d[j - 1] * a[j] / b[j - 1]
    
        Temp_y[i][m - 1] = d[m - 1] / b[m - 1]
        
        
        for j in range(m - 2, -1, -1):
            
            Temp_y[i][j] = (d[j] - c[j] * Temp_y[i][j + 1]) / b[j]
            
    for j in range(m):
        Temp_y[0][j] = Temp_y[1][j] # Boundary condition at the bottom
        Temp_y[n - 1][j] = (k * Temp_y[n - 2][j] + h * dx * t_inf) / (k + h * dx)  # Boundary condition at the top                       
   # print(Temp_y)  
    for i in range(n):
        for j in range(m):
            T[i][j]=Temp_y[i][j]  
    T_last[i][j] = T[100][100]
    T_center[i][j]=T[0][0]
    T_right[i][j]=T[100][0]
    
    # Append time and temperature values to the lists
    time_values.append(t)
    temperature_last.append(T_last[i][j])
    temperature_center.append(T_center[i][j])
    temperature_right.append(T_right[i][j])
    
    #print(f"t[{t}]  = {T_last[i][j]:.2f}")
    t += dt

# Plot the graph

plt.plot(time_values, temperature_right, marker='o', label='Temperature at right(4,0)')
plt.plot(time_values, temperature_center, marker='o', label=' Temperature at center(0,0)')
plt.plot(time_values, temperature_last, marker='o', label='Temperature at corner(4,4)')


plt.xlabel('Time')
plt.ylabel('Temperature')
plt.title('Temperature Change')

plt.legend()  # Adding legend to differentiate between initial and last temperature

plt.grid(True)

plt.show()




