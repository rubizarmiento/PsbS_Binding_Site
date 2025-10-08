
import matplotlib.pyplot as plt
import numpy as np
x = np.linspace(0, 10, 100)
y1 = np.sin(x)
y2 = np.cos(x)
fig, axs = plt.subplots(nrows=2)
axs[0].plot(x, y1)
axs[1].plot(x, y2)

#Save the figure
plt.savefig("/martini/rubiz/Github/PsbS_Binding_Site/4_pairs/analysis/test.png")