import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv('histogram.csv', sep=',', index_col =0)

data.plot(kind='bar')
plt.ylabel('Frequency')
plt.xlabel('Words')
plt.title('Title')

plt.show()