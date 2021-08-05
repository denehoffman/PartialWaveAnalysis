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

pdf = matplotlib.backends.backend_pdf.PdfPages(f"figs_{input_folder.stem}.pdf")

df = pd.read_csv(input_folder / 'fit_results.txt', delimiter='\t', index_col=False)
df.sort_values(['Bin', 'likelihood', 'total_intensity_err'], ascending=[True, False, True], inplace=True)

def mask_first(x):
    result = np.zeros_like(x)
    result[0] = 1
    return result

mask = df.groupby(['Bin'])['Bin'].transform(mask_first).astype(bool)
df_filtered = df.loc[mask]
print(df_filtered.head())

bin_df = pd.read_csv(input_folder / 'bin_info.txt', delimiter='\t')
amplitudes = [column[:-4] for column in df.columns[3:-3].to_list()[::2] if (column.endswith("_INT") and not column.endswith("_AC_INT"))]
amperrors = [column[:-4] for column in df.columns[3:-3].to_list()[1::2] if (column.endswith("_err_INT") and not column.endswith("_err_AC_INT"))]

wave_set = set([amp[:-1] for amp in amplitudes])
wave_dict = {'S': 0, 'P': 1, 'D': 2, 'F': 3, 'G': 4}
waves_sorted = sorted(list(wave_set), key=lambda wave: 100 * wave_dict[wave[0]] + (-1 if wave[-1] == '-' else 1) * int(wave[1]))
# kind of sneaky wave of sorting waves by L-letter and M-number without creating a whole new class


plt.rcParams["figure.figsize"] = (20, 10)
plt.rcParams["font.size"] = 24

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
        plt.errorbar(bin_df['mass'], df_filtered[wave_pos], yerr=df_filtered[wave_pos + "_err"], elinewidth=0.5, fmt='o', color='r', label=r'$+\epsilon$')
    else:
        print("", end='\t')

    wave_neg = wave + "-"
    if wave_neg in amplitudes:
        print("-e", end='\t')
        plt.errorbar(bin_df['mass'], df_filtered[wave_neg], yerr=df_filtered[wave_neg + "_err"], elinewidth=0.5, fmt='o', color='b', label=r'$-\epsilon$')
    else:
        print("", end='\t')
    # Plot total
    print("tot")
    plt.errorbar(bin_df['mass'], df_filtered['total_intensity'], yerr=df_filtered['total_intensity_err'], elinewidth=0.5, fmt='o', color='k', label='Total')
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
        plt.errorbar(bin_df['mass'], df_filtered[wave_pos], yerr=df_filtered[wave_pos + "_err"], elinewidth=0.5, fmt='o', color=colors[i], label=rf"${amp_letter}_{{{amp_m_sign}{amp_m}}}$")
print("tot")
plt.errorbar(bin_df['mass'], df_filtered['total_intensity'], yerr=df_filtered['total_intensity_err'], elinewidth=0.5, fmt='o', color='k', label="Total")
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
        plt.errorbar(bin_df['mass'], df_filtered[wave_neg], yerr=df_filtered[wave_neg + "_err"], elinewidth=0.5, fmt='o', color=colors[i], label=rf"${amp_letter}_{{{amp_m_sign}{amp_m}}}$")
print("tot")
plt.errorbar(bin_df['mass'], df_filtered['total_intensity'], yerr=df_filtered['total_intensity_err'], elinewidth=0.5, fmt='o', color='k', label="Total")
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
        plt.errorbar(bin_df['mass'], df_filtered[wave_pos], yerr=df_filtered[wave_pos + "_err"], elinewidth=0.5, fmt='o', color=colors[i], label=rf"${amp_letter}^+_{{{amp_m_sign}{amp_m}}}$")
    else:
        print("", end='\t')

    wave_neg = wave + "-"
    if wave_neg in amplitudes:
        print("-e", end='\t')
        plt.errorbar(bin_df['mass'], df_filtered[wave_neg], yerr=df_filtered[wave_neg + "_err"], elinewidth=0.5, fmt='s', color=colors[i], label=rf"${amp_letter}^-_{{{amp_m_sign}{amp_m}}}$")
    else:
        print("", end='\t')
    print()
print("tot")
plt.errorbar(bin_df['mass'], df_filtered['total_intensity'], yerr=df_filtered['total_intensity_err'], elinewidth=0.5, fmt='o', color='k', label="Total")
plt.xlim(bin_df['mass'].iloc[0] - 0.1, bin_df['mass'].iloc[-1] + 0.1)
plt.ylim(bottom=0)
plt.legend(loc="upper right")
plt.title("All Waves")
plt.ylabel("Intensity")
plt.xlabel(xlabel)
plt.tight_layout()
pdf.savefig(fig, dpi=300)

pdf.close()
