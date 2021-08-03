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

amplitudes = [column for column in df.columns[3:-3].to_list()[::2] if not (column.endswith("_re") or column.endswith("_im"))]
amplitudes_pos = [a for a in amplitudes if a.endswith("+")]
amplitudes_neg = [a for a in amplitudes if a.endswith("-")]

n_amps = max(len(amplitudes_pos), len(amplitudes_neg))

nrows = int(np.ceil(np.sqrt(n_amps)))
plt.rcParams["figure.figsize"] = (20, 10)
plt.rcParams["font.size"] = 24

def gaussian(x, A, mu, sigma):
    return A * np.exp(-np.power((x - mu), 2) / (2 * np.power(sigma, 2)))

############# Histograms
n_bins = 30
for bin_n in range(len(bin_df)):
    print(bin_n)
    fig, axes = plt.subplots(nrows=2, ncols=n_amps)
    for i, amp in enumerate(amplitudes_pos):
        df_bin = df_bootstrap.loc[df_bootstrap['Bin'] == bin_n]
        df_conv = df_bin[df_bin['Convergence'] == 'C']
        amp_max = df_conv.loc[:, amp].max()
        amp_min = df_conv.loc[:, amp].min()
        correct_value = df_filtered[amp][bin_n]
        correct_gaussian = lambda x, A, sigma: gaussian(x, A, correct_value, sigma)
        entries, bins = np.histogram(df_conv.loc[:, amp], bins=np.linspace(amp_min, amp_max, n_bins))
        bin_centers = np.array([0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)])
        popt, _ = curve_fit(correct_gaussian, xdata=bin_centers, ydata=entries, p0=[1, 1])
        #popt2, _ = curve_fit(gaussian, xdata=bin_centers, ydata=entries, p0=[1, correct_value, 1])
        axes[0, i].bar(bin_centers, entries, width=bins[1] - bins[0], color='navy', label=rf"$\mu$: {int(df_conv.loc[:, amp].mean())}, $\sigma$: {int(df_conv.loc[:, amp].std())}")
        xs = np.linspace(amp_min, amp_max, 2000)
        axes[0, i].plot(xs, correct_gaussian(xs, *popt), color='darkorange', linewidth=2.5, label=rf"$\mu_G$: {int(correct_value)}, $\sigma_G$: {abs(int(popt[1]))}")
        #axes[0, i].plot(xs, gaussian(xs, *popt2), color='darkorange', linestyle='--', linewidth=2.5, label=rf"$\mu_F$: {int(popt2[1])}, $\sigma_F$: {abs(int(popt2[2]))}")
        axes[0, i].legend()
        amp_letter = amplitudes_neg[i].split("::")[-1][0]
        amp_m = amplitudes_neg[i].split("::")[-1][1]
        if int(amp_m) > 0:
            amp_m_sign = amplitudes_neg[i].split("::")[-1][2]
        else:
            amp_m_sign = ""
        axes[0, i].set_title(rf"${amp_letter}_{{{amp_m_sign}{amp_m}}}$")
    for i, amp in enumerate(amplitudes_neg):
        df_bin = df_bootstrap.loc[df_bootstrap['Bin'] == bin_n]
        df_conv = df_bin[df_bin['Convergence'] == 'C']
        amp_max = df_conv.loc[:, amp].max()
        amp_min = df_conv.loc[:, amp].min()
        correct_value = df_filtered[amp][bin_n]
        correct_gaussian = lambda x, A, sigma: gaussian(x, A, correct_value, sigma)
        entries, bins = np.histogram(df_conv.loc[:, amp], bins=np.linspace(amp_min, amp_max, n_bins))
        bin_centers = np.array([0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)])
        popt, _ = curve_fit(correct_gaussian, xdata=bin_centers, ydata=entries, p0=[1, 1])
        #popt2, _ = curve_fit(gaussian, xdata=bin_centers, ydata=entries, p0=[1, correct_value, 1])
        axes[1, i].bar(bin_centers, entries, width=bins[1] - bins[0], color='navy', label=rf"$\mu$: {int(df_conv.loc[:, amp].mean())}, $\sigma$: {int(df_conv.loc[:, amp].std())}")
        xs = np.linspace(amp_min, amp_max, 2000)
        axes[1, i].plot(xs, correct_gaussian(xs, *popt), color='darkorange', linewidth=2.5, label=rf"$\mu_G$: {int(correct_value)}, $\sigma_G$: {abs(int(popt[1]))}")
        #axes[1, i].plot(xs, gaussian(xs, *popt2), color='darkorange', linestyle='--', linewidth=2.5, label=rf"$\mu_F$: {int(popt2[1])}, $\sigma_F$: {abs(int(popt2[2]))}")
        axes[1, i].legend()
        amp_letter = amplitudes_neg[i].split("::")[-1][0]
        amp_m = amplitudes_neg[i].split("::")[-1][1]
        if int(amp_m) > 0:
            amp_m_sign = amplitudes_neg[i].split("::")[-1][2]
        else:
            amp_m_sign = ""
        axes[1, i].set_title(rf"${amp_letter}_{{{amp_m_sign}{amp_m}}}$")
    fig.suptitle(f"Bin {bin_n} Bootstrapped Distributions")
    plt.tight_layout()
    pdf.savefig(fig, dpi=300)
    plt.close()

############# Bootstrap errors
for amp in amplitudes:
    amp_bootstrap_errors = []
    for bin_n in range(len(bin_df)):
        df_bin = df_bootstrap.loc[df_bootstrap['Bin'] == bin_n]
        df_conv = df_bin[df_bin['Convergence'] == 'C']
        amp_bootstrap_errors.append(df_conv.loc[:, amp].std()) # Using mean of distribution for now, no fit
    df_filtered.loc[:, f"{amp}_bootstrap_err"] = amp_bootstrap_errors

amp_bootstrap_errors = []
for bin_n in range(len(bin_df)):
    df_bin = df_bootstrap.loc[df_bootstrap['Bin'] == bin_n]
    df_conv = df_bin[df_bin['Convergence'] == 'C']
    amp_bootstrap_errors.append(df_conv.loc[:, amp].std())
df_filtered.loc[:, f"total_bootstrap_err"] = amp_bootstrap_errors

amperrors = [column for column in df_filtered.columns.to_list() if "bootstrap" in column]
amperrors_pos = [a for a in amperrors if a.endswith("+_bootstrap_err")]
amperrors_neg = [a for a in amperrors if a.endswith("-_bootstrap_err")]


fig, axes = plt.subplots(nrows=nrows, ncols=nrows)
indexes = [idx for idx in np.ndindex(axes.shape)]


print("Plotting Separate Amplitudes")
# Positive
for i in range(len(amplitudes_pos)):
    print(amplitudes_pos[i])
    axes[indexes[i]].errorbar(bin_df['mass'], df_filtered[amplitudes_pos[i]], yerr=df_filtered[amperrors_pos[i]], elinewidth=0.5, fmt='.', color='r', label=r'$+\epsilon$')
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
    axes[indexes[i]].errorbar(bin_df['mass'], df_filtered[amplitudes_neg[i]], yerr=df_filtered[amperrors_neg[i]], elinewidth=0.5, fmt='.', color='b', label=r'$-\epsilon$')
    amp_letter = amplitudes_neg[i].split("::")[-1][0]
    amp_m = amplitudes_neg[i].split("::")[-1][1]
    if int(amp_m) > 0:
        amp_m_sign = amplitudes_neg[i].split("::")[-1][2]
    else:
        amp_m_sign = ""
    axes[indexes[i]].set_title(rf"${amp_letter}_{{{amp_m_sign}{amp_m}}}$")

for i in range(n_amps):
    axes[indexes[i]].errorbar(bin_df['mass'], df_filtered['total_intensity'], yerr=df_filtered['total_bootstrap_err'], elinewidth=0.5, fmt='.', color='k', label='Total')
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
        plt.errorbar(bin_df['mass'], df_filtered[amplitudes_pos[i]], yerr=df_filtered[amperrors_pos[i]], elinewidth=0.5, fmt='.', color=colors[i], label=rf"${amp_letter}_{{{amp_m_sign}{amp_m}}}$")
    plt.errorbar(bin_df['mass'], df_filtered['total_intensity'], yerr=df_filtered['total_bootstrap_err'], elinewidth=0.5, fmt='.', color='k', label="Total")
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
        plt.errorbar(bin_df['mass'], df_filtered[amplitudes_neg[i]], yerr=df_filtered[amperrors_neg[i]], elinewidth=0.5, fmt='.', color=colors[i], label=rf"${amp_letter}_{{{amp_m_sign}{amp_m}}}$")
    plt.errorbar(bin_df['mass'], df_filtered['total_intensity'], yerr=df_filtered['total_bootstrap_err'], elinewidth=0.5, fmt='.', color='k', label="Total")
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
    plt.errorbar(bin_df['mass'], df_filtered[amplitudes_pos[i]], yerr=df_filtered[amperrors_pos[i]], elinewidth=0.5, fmt='.', color=colors[i], label=rf"${amp_letter}^+_{{{amp_m_sign}{amp_m}}}$")

for i in range(len(amplitudes_neg)):
    print(amplitudes_neg[i])
    amp_letter = amplitudes_neg[i].split("::")[-1][0]
    amp_m = amplitudes_neg[i].split("::")[-1][1]
    if int(amp_m) > 0:
        amp_m_sign = amplitudes_neg[i].split("::")[-1][2]
    else:
        amp_m_sign = ""
    plt.errorbar(bin_df['mass'], df_filtered[amplitudes_neg[i]], yerr=df_filtered[amperrors_neg[i]], elinewidth=0.5, fmt='d', color=colors[i], label=rf"${amp_letter}^-_{{{amp_m_sign}{amp_m}}}$")

plt.errorbar(bin_df['mass'], df_filtered['total_intensity'], yerr=df_filtered['total_bootstrap_err'], elinewidth=0.5, fmt='.', color='k', label="Total")
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