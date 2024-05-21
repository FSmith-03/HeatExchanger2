import numpy as np
import matplotlib.pyplot as plt
Y_cold =[]
Y_hot = []
degree = 5
mdot2 = np.arange(0.05, 0.5, 0.01)
for mdoti in mdot2:
    x_hot_data = np.array([0.4826, 0.4340, 0.3924, 0.3507, 0.3021, 0.2535, 0.1979, 0.1493, 0.1111, 0.0694])
    y_hot_data = np.array([0.0944, 0.1662, 0.2297, 0.2820, 0.3294, 0.3856, 0.4447, 0.5006, 0.5311, 0.5615])
    X_hot = np.vander(x_hot_data, degree + 1, increasing=True)
    # Perform least squares polynomial regression
    coefficients, _, _, _ = np.linalg.lstsq(X_hot, y_hot_data, rcond=None)
    Y_hot.append(np.polyval(coefficients[::-1], mdoti))

mdot1 = np.arange(0.15, 0.65, 0.01)
for mdoti in mdot1:
    x_cold_data = np.array([0.6333, 0.6083, 0.5750, 0.5083, 0.4250, 0.3583, 0.3083, 0.2417, 0.1917, 0.1583])
    y_cold_data = np.array([0.1024, 0.1444, 0.1870, 0.2717, 0.3568, 0.4203, 0.4626, 0.5152, 0.5597, 0.5776])
    X_cold = np.vander(x_cold_data, degree + 1, increasing=True)
    # Perform least squares polynomial regression
    coefficients, _, _, _ = np.linalg.lstsq(X_cold, y_cold_data, rcond=None)
    Y_cold.append(np.polyval(coefficients[::-1], mdoti))


plt.scatter(x_hot_data, y_hot_data, s=2, color = "black", label = "Hot Compressor Data")
plt.plot(mdot2, Y_hot, color = "red", label = "Hot Compressor Trends")
plt.scatter(x_cold_data, y_cold_data, s=2, color = "black", label = "Cold Compressor Data")
plt.plot(mdot1, Y_cold, color = "cyan", label = "Cold Compressor Trends")
plt.xlabel('Mdot (kg/s)')
plt.ylabel('Pressure Rise (Bar)')
plt.title('Compressor Pressure Rise Characteristics')
plt.legend(fontsize = 'large')
plt.show()
