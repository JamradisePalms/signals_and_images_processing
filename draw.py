import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

data_sin = pd.read_csv('sin_signal_sampled.csv')
plt.figure()
plt.plot(data_sin.iloc[:, 1])
plt.xlabel('Samples')
plt.ylabel('Amplitude')
plt.show()

data_triangle = pd.read_csv("triangle_signal_sampled.csv")
plt.figure()
plt.plot(data_triangle.iloc[:, 1])
plt.xlabel('Samples')
plt.ylabel('Amplitude')
plt.show()

data_square = pd.read_csv("square_signal_sampled.csv")
plt.figure()
plt.plot(data_square.iloc[:, 1])
plt.xlabel('Samples')
plt.ylabel('Amplitude')
plt.show()

data_dsquare = pd.read_csv("double_square_signal_sampled.csv")
plt.figure()
plt.plot(data_dsquare.iloc[:, 1])
plt.xlabel('Samples')
plt.ylabel('Amplitude')
plt.show()

data_dsquare_complex = pd.read_csv("double_square_signal_complex.csv")
plt.figure()
plt.plot(data_dsquare_complex)
plt.xlabel('freq')
plt.ylabel('Amplitude')
plt.show()