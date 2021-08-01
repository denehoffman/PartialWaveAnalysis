#!/usr/bin/python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import sys
from pathlib import Path

input_folder = Path(sys.argv[1]).resolve()
xlabel = r"$m(K_SK_S)\ GeV/c^2$"
if len(sys.argv) == 3:
    xlabel = rf"${sys.argv[2]}$"

pdf = matplotlib.backends.backend_pdf.PdfPages(f"stats_{input_folder.stem}.pdf")

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
plt.rcParams["figure.figsize"] = (30, 10)
plt.rcParams["font.size"] = 24

print("Plotting Separate Amplitudes")
# Positive
for i in range(len(amplitudes)):
    print(amplitudes[i])
    fig = plt.figure()
    all_runs_by_bin = [df[amplitudes[i]].loc[df['Bin'] == bin_n] for bin_n in bin_df['bin']]
    plt.scatter(bin_df['mass'].iloc[df['Bin']], df[amplitudes[i]], marker='.', color='k', label="Fit Minima")
    plt.violinplot(all_runs_by_bin, bin_df['mass'], widths=bin_df['mass'].iloc[1]-bin_df['mass'].iloc[0], showmeans=True, showextrema=True, showmedians=True)
    plt.scatter(bin_df['mass'].iloc[df_filtered['Bin']], df_filtered[amplitudes[i]], marker='o', color='r', label="Selected Minimum")
    plt.title(amplitudes[i].split("::")[-1])
    plt.xlim(bin_df['mass'].iloc[0] - 0.1, bin_df['mass'].iloc[-1] + 0.1)
    plt.ylim(bottom=0)
    plt.ylabel("Intensity")
    plt.xlabel(xlabel)
    plt.legend()
    plt.tight_layout()
    pdf.savefig(fig, dpi=300)


fig = plt.figure()
all_runs_by_bin = [df['total_intensity'].loc[df['Bin'] == bin_n] for bin_n in bin_df['bin']]
plt.scatter(bin_df['mass'].iloc[df['Bin']], df['total_intensity'], marker='.', color='k', label="Fit Minima")
plt.violinplot(all_runs_by_bin, bin_df['mass'], widths=bin_df['mass'].iloc[1]-bin_df['mass'].iloc[0], showmeans=True, showextrema=True, showmedians=True)
plt.scatter(bin_df['mass'].iloc[df_filtered['Bin']], df_filtered['total_intensity'], marker='o', color='r', label="Selected Minimum")
plt.title('Total Intensity')
plt.xlim(bin_df['mass'].iloc[0] - 0.1, bin_df['mass'].iloc[-1] + 0.1)
plt.ylim(bottom=0)
plt.ylabel("Intensity")
plt.xlabel(xlabel)
plt.legend()
plt.tight_layout()
pdf.savefig(fig, dpi=300)

fig = plt.figure()
all_runs_by_bin = [df['likelihood'].loc[df['Bin'] == bin_n] for bin_n in bin_df['bin']]
plt.scatter(bin_df['mass'].iloc[df['Bin']], df['likelihood'], marker='.', color='k', label="Fit Minima")
plt.violinplot(all_runs_by_bin, bin_df['mass'], widths=bin_df['mass'].iloc[1]-bin_df['mass'].iloc[0], showmeans=True, showextrema=True, showmedians=True)
plt.scatter(bin_df['mass'].iloc[df_filtered['Bin']], df_filtered['likelihood'], marker='o', color='r', label="Selected Minimum")
plt.title('Likelihood')
plt.xlim(bin_df['mass'].iloc[0] - 0.1, bin_df['mass'].iloc[-1] + 0.1)
plt.ylabel("Log-Likelihood")
plt.xlabel(xlabel)
plt.legend()
plt.tight_layout()
pdf.savefig(fig, dpi=300)

fig = plt.figure()
all_runs_by_bin = [df['total_intensity_err'].loc[df['Bin'] == bin_n] for bin_n in bin_df['bin']]
plt.scatter(bin_df['mass'].iloc[df['Bin']], df['total_intensity_err'], marker='.', color='k', label="Fit Minima")
plt.violinplot(all_runs_by_bin, bin_df['mass'], widths=bin_df['mass'].iloc[1]-bin_df['mass'].iloc[0], showmeans=True, showextrema=True, showmedians=True)
plt.scatter(bin_df['mass'].iloc[df_filtered['Bin']], df_filtered['total_intensity_err'], marker='o', color='r', label="Selected Minimum")
plt.title('Error in Total Intensity')
plt.xlim(bin_df['mass'].iloc[0] - 0.1, bin_df['mass'].iloc[-1] + 0.1)
plt.ylabel("Intensity Error")
plt.xlabel(xlabel)
plt.legend()
plt.tight_layout()
pdf.savefig(fig, dpi=300)

pdf.close()
