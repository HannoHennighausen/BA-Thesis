import matplotlib.pyplot as plt
import numpy as np

x, y = np.loadtxt('coor.txt', unpack=True)
n=np.size(x)
# =============================================================================
# x=np.array([])
# y=np.array([])
# with open('coor.txt', 'r') as f:
#     for line in f:
#         first, second = line.split()
#         x.append(first)
#         y.append(second)
# 
# =============================================================================
plt.figure()
plt.scatter(x,y)
# =============================================================================
# 
# for i in enumerate(n):
#     for txt in i:
#         plt.annotate(txt, (x[i], y[i]))
#     
# =============================================================================
plt.show()