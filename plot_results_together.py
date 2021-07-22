import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys


df = pd.read_csv(f'{sys.argv[1]}/fit_results.txt', delimiter='\t', index_col=False)
df.sort_values(['Bin', 'likelihood', 'total_intensity_err'], ascending=[True, False, True], inplace=True)

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
    print(f"Plotting {amplitudes[i]}")
    plt.errorbar(bin_df['mass'].iloc[df_filtered['Bin']], df_filtered[amplitudes[i]], yerr=df_filtered[amperrors[i]], linestyle='-', linewidth=1, elinewidth=0.5, marker='.', color=colors[i], label=amplitudes[i].split("::")[-1])

print("Plotting Total")
plt.errorbar(bin_df['mass'].iloc[df_filtered['Bin']], df_filtered['total_intensity'], yerr=df_filtered['total_intensity_err'], linestyle='-', linewidth=1, elinewidth=0.5, marker='.', color='k', label="Total")
plt.xlim(0.9, 2.5)
plt.ylim(bottom=-100)
plt.legend(loc="upper right")
plt.axhline(0, color='k')
    
plt.savefig(f"fig_{sys.argv[1]}_together.png", dpi=300)
