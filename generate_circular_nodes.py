import numpy as np
import matplotlib.pyplot as plt
   
X = np.array([5])
Y = np.array([5])
nodes = 13


x,y = 5,5
for i in range(nodes):
    alpha = 2*i * np.pi / (nodes-1)
 
    x_i = x + 4 - 4*np.cos(alpha)
    y_i = y + 4*np.sin(alpha)
    
    X = np.append(X, x_i)
    Y = np.append(Y, y_i)
    print(i, ":",alpha, x_i,y_i)
    
for i in range(nodes):
    l_i = np.sqrt((X[i]-X[i+1])**2+(Y[i]-Y[i+1])**2)
    z = (X[i]-X[i+1])+(Y[i]-Y[i+1])
    theta_i = np.angle(z, deg=True)
    
    print(i, ":", l_i, theta_i)

plt.figure()
plt.plot(X,Y, ".")
plt.show()

