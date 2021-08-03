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

pdf = matplotlib.backends.backend_pdf.PdfPages(f"bootstrap_{input_folder.stem}.pdf")

df = pd.read_csv(input_folder / 'fit_results.txt', delimiter='\t', index_col=False)
df.sort_values(['Bin', 'likelihood', 'total_intensity_err'], ascending=[True, False, True], inplace=True)

def mask_first(x):
    result = np.zeros_like(x)
    result[0] = 1
    return result

mask = df.groupby(['Bin'])['Bin'].transform(mask_first).astype(bool)

#############
df_filtered = df.loc[mask]
df_bootstrap = pd.read_csv(input_folder / 'bootstrap.txt', delimiter='\t', index_col=False)
bin_df = pd.read_csv(input_folder / 'bin_info.txt', delimiter='\t')
#############

amplitudes = [column for column in df.columns[3:-3].to_list()[::2] if not (column.endswith("_re") or column.endswith("_im"))]
amplitudes_pos = [a for a in amplitudes if a.endswith("+")]
amplitudes_neg = [a for a in amplitudes if a.endswith("-")]

n_amps = max(len(amplitudes_pos), len(amplitudes_neg))

nrows = int(np.ceil(np.sqrt(n_amps)))
plt.rcParams["figure.figsize"] = (20, 10)
plt.rcParams["font.size"] = 24
print(df_bootstrap.head())

############# Histograms
n_bins = 20
for bin_n in range(len(bin_df)):
    print(bin_n)
    fig, axes = plt.subplots(nrows=2, ncols=n_amps)
    for i, amp in enumerate(amplitudes_pos):
        df_bin = df_bootstrap.loc[df_bootstrap['Bin'] == bin_n]
        df_conv = df_bin[df_bin['Convergence'] == 'C']
        axes[0, i].hist(df_conv.loc[:, amp], bins=n_bins, label=f"Var: {df_conv.loc[:, amp].var()}")
        axes[0, i].legend()
    for i, amp in enumerate(amplitudes_neg):
        df_bin = df_bootstrap.loc[df_bootstrap['Bin'] == bin_n]
        df_conv = df_bin[df_bin['Convergence'] == 'C']
        axes[1, i].hist(df_conv.loc[:, amp], bins=n_bins, label=f"Var: {df_conv.loc[:, amp].var()}")
        axes[1, i].legend()
    plt.tight_layout()
    pdf.savefig(fig, dpi=300)
    plt.close()

############# Bootstrap errors
for amp in amplitudes:
    amp_bootstrap_errors = []
    for bin_n in range(len(bin_df)):
        df_bin = df_bootstrap.loc[df_bootstrap['Bin'] == bin_n]
        df_conv = df_bin[df_bin['Convergence'] == 'C']
        amp_bootstrap_errors.append(df_conv.loc[:, amp].var())
    df_filtered.loc[:, f"{amp}_bootstrap_err"] = amp_bootstrap_errors
    print(df_filtered)

amperrors = [column for column in df_filtered.columns.to_list() if "bootstrap" in column]
amperrors_pos = [a for a in amperrors if a.endswith("+_bootstrap_err")]
amperrors_neg = [a for a in amperrors if a.endswith("-_bootstrap_err")]


fig, axes = plt.subplots(nrows=nrows, ncols=nrows)
indexes = [idx for idx in np.ndindex(axes.shape)]


print("Plotting Separate Amplitudes")
# Positive
for i in range(len(amplitudes_pos)):
    print(amplitudes_pos[i])
    axes[indexes[i]].errorbar(bin_df['mass'].iloc[df_filtered['Bin']], df_filtered[amplitudes_pos[i]], yerr=df_filtered[amperrors_pos[i]], elinewidth=0.5, fmt='.', color='r', label=r'$+\epsilon$')
    amp_letter = amplitudes_pos[i].split("::")[-1][0]
    amp_m = amplitudes_pos[i].split("::")[-1][1]
    if int(amp_m) > 0:
        amp_m_sign = amplitudes_pos[i].split("::")[-1][2]
    else:
        amp_m_sign = ""
    axes[indexes[i]].set_title(rf"${amp_letter}_{{{amp_m_sign}{amp_m}}}$")

# Negative
for i in range(len(amplitudes_neg)):
    print(amplitudes_neg[i])
    axes[indexes[i]].errorbar(bin_df['mass'].iloc[df_filtered['Bin']], df_filtered[amplitudes_neg[i]], yerr=df_filtered[amperrors_neg[i]], elinewidth=0.5, fmt='.', color='b', label=r'$-\epsilon$')
    amp_letter = amplitudes_neg[i].split("::")[-1][0]
    amp_m = amplitudes_neg[i].split("::")[-1][1]
    if int(amp_m) > 0:
        amp_m_sign = amplitudes_neg[i].split("::")[-1][2]
    else:
        amp_m_sign = ""
    axes[indexes[i]].set_title(rf"${amp_letter}_{{{amp_m_sign}{amp_m}}}$")

for i in range(n_amps):
    axes[indexes[i]].errorbar(bin_df['mass'].iloc[df_filtered['Bin']], df_filtered['total_intensity'], yerr=df_filtered['total_intensity_err'], elinewidth=0.5, fmt='.', color='k', label='Total')
    axes[indexes[i]].set_ylim(bottom=0)
    axes[indexes[i]].set_xlim(bin_df['mass'].iloc[0] - 0.1, bin_df['mass'].iloc[-1] + 0.1)
    axes[indexes[i]].set_ylabel("Intensity")
    axes[indexes[i]].set_xlabel(xlabel)
    axes[indexes[i]].legend()

plt.tight_layout()
pdf.savefig(fig, dpi=300)
plt.close()

colors = ['aqua', 'blue', 'chartreuse', 'coral', 'crimson', 'darkblue', 'darkgreen', 'fuchsia', 'gold', 'indigo', 'lime', 'orangered', 'teal', 'sienna']

print("Plotting Combined Amplitudes")
fig = plt.figure()
if len(amplitudes_pos) != 0:
    print("Positive Only")
    for i in range(len(amplitudes_pos)):
        print(amplitudes_pos[i])
        amp_letter = amplitudes_pos[i].split("::")[-1][0]
        amp_m = amplitudes_pos[i].split("::")[-1][1]
        if int(amp_m) > 0:
            amp_m_sign = amplitudes_pos[i].split("::")[-1][2]
        else:
            amp_m_sign = ""
        plt.errorbar(bin_df['mass'].iloc[df_filtered['Bin']], df_filtered[amplitudes_pos[i]], yerr=df_filtered[amperrors_pos[i]], linestyle='-', linewidth=1, elinewidth=0.5, marker='.', color=colors[i], label=rf"${amp_letter}_{{{amp_m_sign}{amp_m}}}$")
    plt.errorbar(bin_df['mass'].iloc[df_filtered['Bin']], df_filtered['total_intensity'], yerr=df_filtered['total_intensity_err'], linestyle='-', linewidth=1, elinewidth=0.5, marker='.', color='k', label="Total")
    plt.xlim(bin_df['mass'].iloc[0] - 0.1, bin_df['mass'].iloc[-1] + 0.1)
    plt.ylim(bottom=-100)
    plt.legend(loc="upper right")
    plt.axhline(0, color='k')
    plt.title("Positive Reflectivity")
    plt.ylabel("Intensity")
    plt.xlabel(xlabel)
    plt.tight_layout()
    pdf.savefig(fig, dpi=300)
    plt.close()
 
fig = plt.figure()
if len(amplitudes_neg) != 0:
    print("Negative Only")
    for i in range(len(amplitudes_neg)):
        print(amplitudes_neg[i])
        amp_letter = amplitudes_neg[i].split("::")[-1][0]
        amp_m = amplitudes_neg[i].split("::")[-1][1]
        if int(amp_m) > 0:
            amp_m_sign = amplitudes_neg[i].split("::")[-1][2]
        else:
            amp_m_sign = ""
        plt.errorbar(bin_df['mass'].iloc[df_filtered['Bin']], df_filtered[amplitudes_neg[i]], yerr=df_filtered[amperrors_neg[i]], linestyle='-', linewidth=1, elinewidth=0.5, marker='.', color=colors[i], label=rf"${amp_letter}_{{{amp_m_sign}{amp_m}}}$")
    plt.errorbar(bin_df['mass'].iloc[df_filtered['Bin']], df_filtered['total_intensity'], yerr=df_filtered['total_intensity_err'], linestyle='-', linewidth=1, elinewidth=0.5, marker='.', color='k', label="Total")
    plt.xlim(bin_df['mass'].iloc[0] - 0.1, bin_df['mass'].iloc[-1] + 0.1)
    plt.ylim(bottom=-100)
    plt.legend(loc="upper right")
    plt.axhline(0, color='k')
    plt.title("Negative Reflectivity")
    plt.ylabel("Intensity")
    plt.xlabel(xlabel)
    plt.tight_layout()
    pdf.savefig(fig, dpi=300)
    plt.close()
 
print("Positive and Negative Combined")
fig = plt.figure()
for i in range(len(amplitudes_pos)):
    print(amplitudes_pos[i])
    amp_letter = amplitudes_pos[i].split("::")[-1][0]
    amp_m = amplitudes_pos[i].split("::")[-1][1]
    if int(amp_m) > 0:
        amp_m_sign = amplitudes_pos[i].split("::")[-1][2]
    else:
        amp_m_sign = ""
    plt.errorbar(bin_df['mass'].iloc[df_filtered['Bin']], df_filtered[amplitudes_pos[i]], yerr=df_filtered[amperrors_pos[i]], linestyle='-', linewidth=1, elinewidth=0.5, marker='.', color=colors[i], label=rf"${amp_letter}^+_{{{amp_m_sign}{amp_m}}}$")

for i in range(len(amplitudes_neg)):
    print(amplitudes_neg[i])
    amp_letter = amplitudes_neg[i].split("::")[-1][0]
    amp_m = amplitudes_neg[i].split("::")[-1][1]
    if int(amp_m) > 0:
        amp_m_sign = amplitudes_neg[i].split("::")[-1][2]
    else:
        amp_m_sign = ""
    plt.errorbar(bin_df['mass'].iloc[df_filtered['Bin']], df_filtered[amplitudes_neg[i]], yerr=df_filtered[amperrors_neg[i]], linestyle='--', linewidth=1, elinewidth=0.5, marker='.', color=colors[i], label=rf"${amp_letter}^-_{{{amp_m_sign}{amp_m}}}$")

plt.errorbar(bin_df['mass'].iloc[df_filtered['Bin']], df_filtered['total_intensity'], yerr=df_filtered['total_intensity_err'], linestyle='-', linewidth=1, elinewidth=0.5, marker='.', color='k', label="Total")
plt.xlim(bin_df['mass'].iloc[0] - 0.1, bin_df['mass'].iloc[-1] + 0.1)
plt.ylim(bottom=-100)
plt.legend(loc="upper right")
plt.axhline(0, color='k')
plt.title("Fit Results")
plt.ylabel("Intensity")
plt.xlabel(xlabel)
plt.tight_layout()
pdf.savefig(fig, dpi=300)
plt.close()

pdf.close()
