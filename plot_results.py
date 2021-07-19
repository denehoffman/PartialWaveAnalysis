import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('output_twopsangles/fit_results.txt', delimiter='\t', index_col=False)
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

nrows = int(np.ceil(np.sqrt(len(amplitudes) + 1)))
fig, axes = plt.subplots(nrows=nrows, ncols=nrows)
fig.tight_layout()
plt.rcParams["figure.figsize"] = (100, 100)
indexes = [idx for idx in np.ndindex(axes.shape)]

colors = ['aqua', 'blue', 'chartreuse', 'coral', 'crimson', 'darkblue', 'darkgreen', 'fuchsia', 'gold', 'indigo', 'lime', 'orangered', 'teal', 'sienna']
for i in range(len(amplitudes)):
    axes[indexes[i]].errorbar(bin_df['mass'].iloc[df_filtered['Bin']], df_filtered[amplitudes[i]], yerr=df_filtered[amperrors[i]], linestyle='-', linewidth=1, elinewidth=0.5, marker='.', color=colors[i])
    axes[indexes[i]].set_title(amplitudes[i].split("::")[-1])
    axes[indexes[i]].set_xlim(0.9, 2.3)
    axes[indexes[i]].set_ylim(bottom=0)

axes[indexes[len(amplitudes)]].errorbar(bin_df['mass'].iloc[df_filtered['Bin']], df_filtered['total_intensity'], yerr=df_filtered['total_intensity_err'], linestyle='-', linewidth=1, elinewidth=0.5, marker='.', color='k')
axes[indexes[len(amplitudes)]].set_title('Total')
axes[indexes[len(amplitudes)]].set_xlim(0.9, 2.3)
axes[indexes[len(amplitudes)]].set_ylim(bottom=0)
    
plt.savefig("fig_test_twopsangles.png", dpi=300)
