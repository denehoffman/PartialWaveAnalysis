import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


df = pd.read_csv('output/fit_results.txt', delimiter='\t', index_col=False)
print(df)
df.sort_values(['Bin', 'likelihood', 'total_intensity_err'], ascending=[True, False, True], inplace=True)
print(df)

def mask_first(x):
    result = np.zeros_like(x)
    result[0] = 1
    return result

mask = df.groupby(['Bin'])['Bin'].transform(mask_first).astype(bool)
df_filtered = df.loc[mask]

bin_df = pd.read_csv('output/bin_info.txt', delimiter='\t')
amplitudes = df.columns[2:-3].to_list()[::2]
amperrors = df.columns[2:-3].to_list()[1::2]
colors = ['aqua', 'blue', 'chartreuse', 'coral', 'crimson', 'darkblue', 'darkgreen', 'fuchsia', 'gold', 'indigo', 'lime', 'orangered', 'teal', 'sienna']
for i in range(len(amplitudes)):
    plt.errorbar(bin_df['mass'].iloc[df_filtered['Bin']], df_filtered[amplitudes[i]], yerr=df_filtered[amperrors[i]], fmt='.', label=amplitudes[i].split("::")[-1], color=colors[i])

plt.errorbar(bin_df['mass'].iloc[df_filtered['Bin']], df_filtered['total_intensity'], yerr=df_filtered['total_intensity_err'], fmt='.', color='k', label="total")
    
plt.legend(loc="upper right")
plt.savefig("fig_test.png")
