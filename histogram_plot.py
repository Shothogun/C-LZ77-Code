import pandas as pd
import matplotlib.pyplot as plt

file = input()
data = pd.read_csv(file+'.csv', sep=',', index_col =0)

data.plot(kind='bar')
plt.ylabel('Frequency')
plt.xlabel('Words')
plt.title('Title')

plt.show()