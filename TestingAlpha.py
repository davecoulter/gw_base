import matplotlib.pyplot as plt
import numpy as np

c = 3*10**8
x = np.linspace(0,c*0.9,10000)
y = np.sqrt((c+x)/(c-x)) - 1
cap = np.linspace(0,1, 10000)

plt.figure(1)
plt.scatter(x,y, c = cap, cmap = "viridis",vmin=0, vmax=1)
plt.xlabel("Recessional Velocity [m/s]")
plt.ylabel("Redshift")

cbar = plt.colorbar()
cbar.set_label("hello")
plt.savefig("Redshift vs Velocity.png")