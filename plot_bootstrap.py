#!/usr/bin/python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import sys
from scipy.optimize import curve_fit
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
df_filtered.set_index('Bin', inplace=True)
df_bootstrap = pd.read_csv(input_folder / 'bootstrap.txt', delimiter='\t', index_col=False)
bin_df = pd.read_csv(input_folder / 'bin_info.txt', delimiter='\t')
#############

amplitudes = [column[:-4] for column in df.columns[3:-3].to_list()[::2] if (column.endswith("_INT") and not column.endswith("_AC_INT"))]

wave_set = set([amp[:-1] for amp in amplitudes])
wave_dict = {'S': 0, 'P': 1, 'D': 2, 'F': 3, 'G': 4}
waves_sorted = sorted(list(wave_set), key=lambda wave: 100 * wave_dict[wave[0]] + (-1 if wave[-1] == '-' else 1) * int(wave[1]))
# kind of sneaky wave of sorting waves by L-letter and M-number without creating a whole new class

n_amps = len(waves_sorted) + 1 # plus 1 for the total plot

plt.rcParams["figure.figsize"] = (20, 10)
plt.rcParams["font.size"] = 24

############# Histograms
n_bins = 30
for bin_n in range(len(bin_df)):
    print(bin_n)
    fig, axes = plt.subplots(nrows=2, ncols=n_amps)
    for i, wave in enumerate(waves_sorted):
        print(wave, end='\t')
        df_bin = df_bootstrap.loc[df_bootstrap['Bin'] == bin_n]
        df_conv = df_bin[df_bin['Convergence'] == 'C']
        wave_pos = wave + "+"
        if wave_pos in amplitudes:
            print("+e", end='\t')
            amp_max = df_conv.loc[:, wave_pos].max()
            amp_min = df_conv.loc[:, wave_pos].min()
            correct_value = df_filtered[wave_pos][bin_n]
            entries, bins = np.histogram(df_conv.loc[:, wave_pos], bins=np.linspace(amp_min, amp_max, n_bins))
            bin_centers = np.array([0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)])
            axes[0, i].bar(bin_centers, entries, width=bins[1] - bins[0], color='navy', label=rf"$\mu$: {int(df_conv.loc[:, wave_pos].mean())}, $\sigma$: {int(df_conv.loc[:, wave_pos].std())}")
            axes[0, i].legend()
            amp_letter = wave[0] 
            amp_m = wave[1]
            if int(amp_m) > 0:
                amp_m_sign = wave[2]
            else:
                amp_m_sign = ""
            axes[0, i].set_title(rf"${amp_letter}^+_{{{amp_m_sign}{amp_m}}}$")

        wave_neg = wave + "-"
        if wave_neg in amplitudes:
            print("-e", end='\t')
            amp_max = df_conv.loc[:, wave_neg].max()
            amp_min = df_conv.loc[:, wave_neg].min()
            correct_value = df_filtered[wave_neg][bin_n]
            entries, bins = np.histogram(df_conv.loc[:, wave_neg], bins=np.linspace(amp_min, amp_max, n_bins))
            bin_centers = np.array([0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)])
            axes[1, i].bar(bin_centers, entries, width=bins[1] - bins[0], color='navy', label=rf"$\mu$: {int(df_conv.loc[:, wave_neg].mean())}, $\sigma$: {int(df_conv.loc[:, wave_neg].std())}")
            axes[1, i].legend()
            amp_letter = wave[0] 
            amp_m = wave[1]
            if int(amp_m) > 0:
                amp_m_sign = wave[2]
            else:
                amp_m_sign = ""
            axes[1, i].set_title(rf"${amp_letter}^-_{{{amp_m_sign}{amp_m}}}$")
        print()
    print("tot")
    amp_max = df_conv.loc[:, 'total_intensity'].max()
    amp_min = df_conv.loc[:, 'total_intensity'].min()
    correct_value = df_filtered['total_intensity'][bin_n]
    entries, bins = np.histogram(df_conv.loc[:, 'total_intensity'], bins=np.linspace(amp_min, amp_max, n_bins))
    bin_centers = np.array([0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)])
    axes[1, -1].bar(bin_centers, entries, width=bins[1] - bins[0], color='navy', label=rf"$\mu$: {int(df_conv.loc[:, 'total_intensity'].mean())}, $\sigma$: {int(df_conv.loc[:, 'total_intensity'].std())}")
    axes[1, -1].legend()
    axes[1, -1].set_title("Total Intensity")
        
    fig.suptitle(f"Bin {bin_n} Bootstrapped Distributions")
    plt.tight_layout()
    pdf.savefig(fig, dpi=300)
    plt.close()

############# Bootstrap errors
for amp in amplitudes:
    alpha = 0.05
    amp_bootstrap_errors = []
    amp_bootstrap_CIL = []
    amp_bootstrap_CIU = []
    for bin_n in range(len(bin_df)):
        fit_value = df_filtered[amp][bin_n]
        df_bin = df_bootstrap.loc[df_bootstrap['Bin'] == bin_n]
        df_conv = df_bin[df_bin['Convergence'] == 'C']
        amp_bootstrap_errors.append(df_conv.loc[:, amp].std()) # Using mean of distribution for now
        amp_bootstrap_CIL.append(2 * fit_value - np.quantile(df_conv.loc[:, amp], 1 - alpha / 2))
        amp_bootstrap_CIU.append(2 * fit_value - np.quantile(df_conv.loc[:, amp], alpha / 2))

    df_filtered.loc[:, f"{amp}_bootstrap_err"] = amp_bootstrap_errors
    df_filtered.loc[:, f"{amp}_bootstrap_CIL"] = amp_bootstrap_CIL
    df_filtered.loc[:, f"{amp}_bootstrap_CIU"] = amp_bootstrap_CIU

amp_bootstrap_errors = []
amp_bootstrap_CIL = []
amp_bootstrap_CIU = []
alpha = 0.05
for bin_n in range(len(bin_df)):
    fit_value = df_filtered['total_intensity'][bin_n]
    df_bin = df_bootstrap.loc[df_bootstrap['Bin'] == bin_n]
    df_conv = df_bin[df_bin['Convergence'] == 'C']
    amp_bootstrap_errors.append(df_conv.loc[:, 'total_intensity'].std())
    amp_bootstrap_CIL.append(2 * fit_value - np.quantile(df_conv.loc[:, 'total_intensity'], 1 - alpha / 2))
    amp_bootstrap_CIU.append(2 * fit_value - np.quantile(df_conv.loc[:, 'total_intensity'], alpha / 2))
df_filtered.loc[:, f"total_bootstrap_err"] = amp_bootstrap_errors
df_filtered.loc[:, f"total_bootstrap_CIL"] = amp_bootstrap_CIL
df_filtered.loc[:, f"total_bootstrap_CIU"] = amp_bootstrap_CIU

############ Plot new error bars on fits
print("Plotting Separate Amplitudes")
for wave in waves_sorted:
    fig = plt.figure()
    print(wave, end='\t')
    # wave = S0, P1+, D2-, etc. -- no reflectivity info
    amp_letter = wave[0]
    amp_m = wave[1]
    if int(amp_m) > 0:
        amp_m_sign = wave[2]
    else:
        amp_m_sign = ""

    wave_pos = wave + "+"
    if wave_pos in amplitudes:
        print("+e", end='\t')
        plt.errorbar(bin_df['mass'], df_filtered[wave_pos], yerr=df_filtered[wave_pos + "_bootstrap_err"], elinewidth=0.5, fmt='o', color='r', label=r'$+\epsilon$')
        plt.fill_between(bin_df['mass'], df_filtered[wave_pos + "_bootstrap_CIL"], df_filtered[wave_pos + "_bootstrap_CIU"], color='r', alpha=0.1)
    else:
        print("", end='\t')

    wave_neg = wave + "-"
    if wave_neg in amplitudes:
        print("-e", end='\t')
        plt.errorbar(bin_df['mass'], df_filtered[wave_neg], yerr=df_filtered[wave_neg + "_bootstrap_err"], elinewidth=0.5, fmt='o', color='b', label=r'$-\epsilon$')
        plt.fill_between(bin_df['mass'], df_filtered[wave_neg + "_bootstrap_CIL"], df_filtered[wave_neg + "_bootstrap_CIU"], color='b', alpha=0.1)
    else:
        print("", end='\t')
    # Plot total
    print("tot")
    plt.errorbar(bin_df['mass'], df_filtered['total_intensity'], yerr=df_filtered['total_bootstrap_err'], elinewidth=0.5, fmt='o', color='k', label='Total')
    plt.fill_between(bin_df['mass'], df_filtered["total_bootstrap_CIL"], df_filtered["total_bootstrap_CIU"], color='k', alpha=0.1)
    plt.title(rf"${amp_letter}_{{{amp_m_sign}{amp_m}}}$")
    plt.ylim(bottom=0)
    plt.xlim(bin_df['mass'].iloc[0] - 0.1, bin_df['mass'].iloc[-1] + 0.1)
    plt.ylabel("Intensity")
    plt.xlabel(xlabel)
    plt.legend()
    plt.tight_layout()
    pdf.savefig(fig, dpi=300)

colors = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple", "tab:brown", "tab:pink", "tab:gray", "tab:olive", "tab:cyan"]

print("Plotting Combined Amplitudes")
print("Positive Reflectivity")
fig = plt.figure()
for i, wave in enumerate(waves_sorted):
    amp_letter = wave[0]
    amp_m = wave[1]
    if int(amp_m) > 0:
        amp_m_sign = wave[2]
    else:
        amp_m_sign = ""

    wave_pos = wave + "+"
    if wave_pos in amplitudes:
        print(wave + "\t+e")
        plt.errorbar(bin_df['mass'], df_filtered[wave_pos], yerr=df_filtered[wave_pos + "_bootstrap_err"], elinewidth=0.5, fmt='o', color=colors[i], label=rf"${amp_letter}_{{{amp_m_sign}{amp_m}}}$")
print("tot")
plt.errorbar(bin_df['mass'], df_filtered['total_intensity'], yerr=df_filtered['total_bootstrap_err'], elinewidth=0.5, fmt='o', color='k', label="Total")
plt.xlim(bin_df['mass'].iloc[0] - 0.1, bin_df['mass'].iloc[-1] + 0.1)
plt.ylim(bottom=0)
plt.legend(loc="upper right")
plt.title("Positive Reflectivity")
plt.ylabel("Intensity")
plt.xlabel(xlabel)
plt.tight_layout()
pdf.savefig(fig, dpi=300)

print("Negative Reflectivity")
fig = plt.figure()
for i, wave in enumerate(waves_sorted):
    amp_letter = wave[0]
    amp_m = wave[1]
    if int(amp_m) > 0:
        amp_m_sign = wave[2]
    else:
        amp_m_sign = ""

    wave_neg = wave + "-"
    if wave_neg in amplitudes:
        print(wave + "\t\t-e")
        plt.errorbar(bin_df['mass'], df_filtered[wave_neg], yerr=df_filtered[wave_neg + "_bootstrap_err"], elinewidth=0.5, fmt='o', color=colors[i], label=rf"${amp_letter}_{{{amp_m_sign}{amp_m}}}$")
print("tot")
plt.errorbar(bin_df['mass'], df_filtered['total_intensity'], yerr=df_filtered['total_bootstrap_err'], elinewidth=0.5, fmt='o', color='k', label="Total")
plt.xlim(bin_df['mass'].iloc[0] - 0.1, bin_df['mass'].iloc[-1] + 0.1)
plt.ylim(bottom=0)
plt.legend(loc="upper right")
plt.title("Negative Reflectivity")
plt.ylabel("Intensity")
plt.xlabel(xlabel)
plt.tight_layout()
pdf.savefig(fig, dpi=300)

print("Positive and Negative Reflectivity")
fig = plt.figure()
for i, wave in enumerate(waves_sorted):
    amp_letter = wave[0]
    amp_m = wave[1]
    if int(amp_m) > 0:
        amp_m_sign = wave[2]
    else:
        amp_m_sign = ""

    wave_pos = wave + "+"
    print(wave, end='\t')
    if wave_pos in amplitudes:
        print("+e", end='\t')
        plt.errorbar(bin_df['mass'], df_filtered[wave_pos], yerr=df_filtered[wave_pos + "_bootstrap_err"], elinewidth=0.5, fmt='o', color=colors[i], label=rf"${amp_letter}^+_{{{amp_m_sign}{amp_m}}}$")
    else:
        print("", end='\t')

    wave_neg = wave + "-"
    if wave_neg in amplitudes:
        print("-e", end='\t')
        plt.errorbar(bin_df['mass'], df_filtered[wave_neg], yerr=df_filtered[wave_neg + "_bootstrap_err"], elinewidth=0.5, fmt='s', color=colors[i], label=rf"${amp_letter}^-_{{{amp_m_sign}{amp_m}}}$")
    else:
        print("", end='\t')
    print()
print("tot")
plt.errorbar(bin_df['mass'], df_filtered['total_intensity'], yerr=df_filtered['total_bootstrap_err'], elinewidth=0.5, fmt='o', color='k', label="Total")
plt.xlim(bin_df['mass'].iloc[0] - 0.1, bin_df['mass'].iloc[-1] + 0.1)
plt.ylim(bottom=0)
plt.legend(loc="upper right")
plt.title("All Waves")
plt.ylabel("Intensity")
plt.xlabel(xlabel)
plt.tight_layout()
pdf.savefig(fig, dpi=300)

pdf.close()
