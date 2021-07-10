import matplotlib.pyplot as plt
import numpy as np
import csv
from datetime import datetime

if False:
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

# x = np.zeros(10)
# index = 0
# for i in range(6):
#     x[index] = 1
#     index = index + 1
# print(x)
# print(x[:index])
with open("local_data/DES_allBands.csv", mode='r') as csv_file:
    galaxy_code = csv.DictReader(csv_file)
    index = 0
    i = 0
    perc = 5
    perc_now = perc
    total_file = 529270
    print("Loading Galaxy")
    now = datetime.now()
    for row in galaxy_code:
        index = index + 1
        if i/total_file >= perc_now/100:
            print(perc_now, "%", datetime.now() - now)
            perc_now = perc_now + perc
        i = i + 1
print(index)
print(i)