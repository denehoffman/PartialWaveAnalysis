#!/usr/bin/python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
from pathlib import Path

input_folder = Path(sys.argv[1]).resolve()

df = pd.read_csv(input_folder / 'fit_results.txt', delimiter='\t', index_col=False)
df.sort_values(['Bin', 'likelihood', 'total_intensity_err'], ascending=[True, False, True], inplace=True)

def mask_first(x):
    result = np.zeros_like(x)
    result[0] = 1
    return result

mask = df.groupby(['Bin'])['Bin'].transform(mask_first).astype(bool)
df_filtered = df.loc[mask]

bin_df = pd.read_csv(input_folder / 'bin_info.txt', delimiter='\t')
amplitudes = df.columns[3:-3].to_list()[::2]
amplitudes_pos = [a for a in amplitudes if a.endswith("+")]
amplitudes_neg = [a for a in amplitudes if a.endswith("-")]
amperrors = df.columns[3:-3].to_list()[1::2]
amperrors_pos = [a for a in amperrors if a.endswith("+_err")]
amperrors_neg = [a for a in amperrors if a.endswith("-_err")]

n_amps = max(len(amplitudes_pos), len(amplitudes_neg))

nrows = int(np.ceil(np.sqrt(n_amps + 1)))
plt.rcParams["figure.figsize"] = (10, 10)
fig, axes = plt.subplots(nrows=nrows, ncols=nrows)
fig.tight_layout()
indexes = [idx for idx in np.ndindex(axes.shape)]

print("Plotting Separate Amplitudes")
# Positive
for i in range(len(amplitudes_pos)):
    print(amplitudes_pos[i])
    axes[indexes[i]].errorbar(bin_df['mass'].iloc[df_filtered['Bin']], df_filtered[amplitudes_pos[i]], yerr=df_filtered[amperrors_pos[i]], linestyle='-', linewidth=1, elinewidth=0.5, marker='.', color='k')
    axes[indexes[i]].set_title(amplitudes_pos[i].split("::")[-1][:-1])
    axes[indexes[i]].set_xlim(bin_df['mass'].iloc[0] - 0.1, bin_df['mass'].iloc[-1] + 0.1)
    axes[indexes[i]].set_ylim(bottom=0)

# Negative
for i in range(len(amplitudes_neg)):
    print(amplitudes_neg[i])
    axes[indexes[i]].errorbar(bin_df['mass'].iloc[df_filtered['Bin']], df_filtered[amplitudes_neg[i]], yerr=df_filtered[amperrors_neg[i]], linestyle='-', linewidth=1, elinewidth=0.5, marker='.', color='r')
    axes[indexes[i]].set_title(amplitudes_neg[i].split("::")[-1][:-1])
    axes[indexes[i]].set_xlim(bin_df['mass'].iloc[0] - 0.1, bin_df['mass'].iloc[-1] + 0.1)
    axes[indexes[i]].set_ylim(bottom=0)


axes[indexes[n_amps]].errorbar(bin_df['mass'].iloc[df_filtered['Bin']], df_filtered['total_intensity'], yerr=df_filtered['total_intensity_err'], linestyle='-', linewidth=1, elinewidth=0.5, marker='.', color='b')
axes[indexes[n_amps]].set_title('Total')
axes[indexes[n_amps]].set_xlim(bin_df['mass'].iloc[0] - 0.1, bin_df['mass'].iloc[-1] + 0.1)
axes[indexes[n_amps]].set_ylim(bottom=0)

plt.savefig(f"fig_{input_folder.stem}.png", dpi=300)
plt.close()

colors = ['aqua', 'blue', 'chartreuse', 'coral', 'crimson', 'darkblue', 'darkgreen', 'fuchsia', 'gold', 'indigo', 'lime', 'orangered', 'teal', 'sienna']

print("Plotting Combined Amplitudes")
plt.figure()
if len(amplitudes_pos) != 0:
    print("Positive Only")
    for i in range(len(amplitudes_pos)):
        print(amplitudes_pos[i])
        plt.errorbar(bin_df['mass'].iloc[df_filtered['Bin']], df_filtered[amplitudes_pos[i]], yerr=df_filtered[amperrors_pos[i]], linestyle='-', linewidth=1, elinewidth=0.5, marker='.', color=colors[i], label=amplitudes_pos[i].split("::")[-1][:-1])
    plt.errorbar(bin_df['mass'].iloc[df_filtered['Bin']], df_filtered['total_intensity'], yerr=df_filtered['total_intensity_err'], linestyle='-', linewidth=1, elinewidth=0.5, marker='.', color='k', label="Total")
    plt.xlim(bin_df['mass'].iloc[0] - 0.1, bin_df['mass'].iloc[-1] + 0.1)
    plt.ylim(bottom=-100)
    plt.legend(loc="upper right")
    plt.axhline(0, color='k')
    plt.title("Positive Reflectivity")
    plt.savefig(f"fig_{input_folder.stem}_together_pos", dpi=300)
    plt.close()
 
plt.figure()
if len(amplitudes_neg) != 0:
    print("Negative Only")
    for i in range(len(amplitudes_neg)):
        print(amplitudes_neg[i])
        plt.errorbar(bin_df['mass'].iloc[df_filtered['Bin']], df_filtered[amplitudes_neg[i]], yerr=df_filtered[amperrors_neg[i]], linestyle='-', linewidth=1, elinewidth=0.5, marker='.', color=colors[i], label=amplitudes_neg[i].split("::")[-1][:-1])
    plt.errorbar(bin_df['mass'].iloc[df_filtered['Bin']], df_filtered['total_intensity'], yerr=df_filtered['total_intensity_err'], linestyle='-', linewidth=1, elinewidth=0.5, marker='.', color='k', label="Total")
    plt.xlim(bin_df['mass'].iloc[0] - 0.1, bin_df['mass'].iloc[-1] + 0.1)
    plt.ylim(bottom=-100)
    plt.legend(loc="upper right")
    plt.axhline(0, color='k')
    plt.title("Negative Reflectivity")
    plt.savefig(f"fig_{input_folder.stem}_together_neg", dpi=300)
    plt.close()
 
print("Positive and Negative Combined")
plt.figure()
for i in range(len(amplitudes_pos)):
    print(amplitudes_pos[i])
    plt.errorbar(bin_df['mass'].iloc[df_filtered['Bin']], df_filtered[amplitudes_pos[i]], yerr=df_filtered[amperrors_pos[i]], linestyle='-', linewidth=1, elinewidth=0.5, marker='.', color=colors[i], label=amplitudes_pos[i].split("::")[-1])

for i in range(len(amplitudes_neg)):
    print(amplitudes_neg[i])
    plt.errorbar(bin_df['mass'].iloc[df_filtered['Bin']], df_filtered[amplitudes_neg[i]], yerr=df_filtered[amperrors_neg[i]], linestyle='--', linewidth=1, elinewidth=0.5, marker='.', color=colors[i], label=amplitudes_neg[i].split("::")[-1])

plt.errorbar(bin_df['mass'].iloc[df_filtered['Bin']], df_filtered['total_intensity'], yerr=df_filtered['total_intensity_err'], linestyle='-', linewidth=1, elinewidth=0.5, marker='.', color='k', label="Total")
plt.xlim(bin_df['mass'].iloc[0] - 0.1, bin_df['mass'].iloc[-1] + 0.1)
plt.ylim(bottom=-100)
plt.legend(loc="upper right")
plt.axhline(0, color='k')
plt.title("Fit Results")
plt.savefig(f"fig_{input_folder.stem}_together", dpi=300)
plt.close()
