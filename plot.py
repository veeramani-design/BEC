import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# Load the data from the text file
data = np.loadtxt('imag1d-th-den.txt')

# Separate the columns into x and y
x = data[:, 0]  # First column (x values)
y = data[:, 1]  # Second column (y values)


# Create the plot
# plt.figure(figsize=(8, 5))
p=PdfPages("results/1dplot.pdf")
plt.rcParams["figure.figsize"]=[16.00,9.00]
plt.rcParams["figure.autolayout"]=True
plt.plot(x, y, label='Data')

# Add labels, title, and grid
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.title('1D')
plt.legend()
plt.grid()
plt.savefig(p,format='pdf')
# Show the plot
# plt.show()
plt.clf()
p.close()