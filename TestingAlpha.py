import matplotlib.pyplot as plt
import numpy as np

plt.figure(1)
c = 3*10**8
x = np.linspace(0,c*0.9,10000)
y = np.sqrt((c+x)/(c-x)) - 1
plt.xlabel("Recessional Velocity [m/s]")
plt.ylabel("Redshift")

plt.scatter(x,y)
plt.savefig("Redshift vs Velocity.png")