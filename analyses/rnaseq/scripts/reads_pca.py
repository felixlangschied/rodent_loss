import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
pd.set_option('mode.chained_assignment', None)
import sys


path = sys.argv[1]


#organism = path.split('/')[-1].split('_')[0]
organism = sys.argv[2]
mirna = sys.argv[3]
outpath = sys.argv[4]

print(organism, mirna)

# Sample data
x = [1, 2, 3, 4, 5]
y = [10, 20, 25, 30, 40]
sizes = [100, 200, 250, 300, 400]  # Sizes for the scatter points
colors = ['red', 'green', 'blue', 'orange', 'purple']  # Colors for each point

# Create a scatter plot
plt.figure(figsize=(8, 6))
plt.scatter(x, y, s=sizes, c=colors, alpha=0.7, edgecolors="w", linewidth=2)

# Add labels and title
plt.title('Sample Scatter Plot', fontsize=14)
plt.xlabel('X Axis Label')
plt.ylabel('Y Axis Label')

# Show grid
plt.grid(True)

# Display the plot
plt.savefig(outpath)



# #######################
# rawdf = pd.read_csv(path, sep='\t')
# conds = [name for name in rawdf.columns if '_' in name]

# df = rawdf.filter(conds)
# df = df.astype(int)
# tdf = df.transpose()
# #display(tdf)

# # Separating out the features
# x = tdf.loc[:, tdf.columns].values
# y = tdf.index# Standardizing the features
# x = StandardScaler().fit_transform(x)
# pca = PCA(n_components=comps)
# pca.fit(x)

# exp1, exp2 = pca.explained_variance_ratio_

# exP1 = round(exp1*100, 1)
# exP2 = round(exp2*100, 1)

# principalComponents = pca.fit_transform(x)
# principalDf = pd.DataFrame(data = principalComponents, columns = ['Dim 1', 'Dim 2'], index=tdf.index)
# pdf = principalDf.reset_index()
# pdf['probe'] = pdf['index'].apply(lambda x: x.split('_')[0])
# pdf['probe'][pdf['probe'] == 'Neg-Ctl'] = 'CTRL'
# #display(pdf)

# sns.set(rc={'figure.figsize':(6,4), 'ytick.left': True, 'xtick.bottom': True}, font_scale = 1.3, style='whitegrid')
# sns.scatterplot(data=pdf, x='Dim 1', y='Dim 2', hue='probe', s=90)
# plt.xlabel(f'Dim 1 ({exP1}%)')
# plt.ylabel(f'Dim 2 ({exP2}%)')
# plt.ylim([-100, 125])
# plt.xlim([-100, 125])
# plt.tight_layout()
# plt.savefig(snakemake.output[0])
# plt.show()