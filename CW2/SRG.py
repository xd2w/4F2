import numpy as np
from matplotlib import pyplot as plt

# x = np.linspace(-2, 2, 100)
# y = np.linspace(-2, 2, 100)

# X, Y = np.meshgrid(x, y)

# s = X + 1j*Y

# z = (X**2 + Y**2) < 1

# G = 1/(s+1)

# plt.contour(np.real(G), np.imag(G), s)
# plt.show()


t = np.linspace(0, 2*np.pi, 100)

x = np.cos(t)/(np.cos(t/2)**2)
y = -np.sin(t)/(np.cos(t/2)**2)

plt.plot(x, y)
plt.show()