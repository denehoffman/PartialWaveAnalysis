#!/usr/bin/python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import sys
from pathlib import Path
import uproot
import awkward as ak
import argparse
from tqdm import tqdm
from scipy.special import gammaln
import scipy.stats as st
from simple_term_menu import TerminalMenu


def plot_menu(root):
    plot_menu_title = "Select a fit_results.txt file:"
    fits = [p for p in root.glob("*::fit_results.txt")]
    plot_menu_items = ["Cancel"] + [fit.name[:-17] for fit in fits]
    plot_menu_cursor = "> "
    plot_menu_cursor_style = ("fg_red", "bold")
    plot_menu_style = ("bg_black", "fg_green")
    plot_menu = TerminalMenu(
        menu_entries=plot_menu_items,
        title=plot_menu_title,
        menu_cursor=plot_menu_cursor,
        menu_cursor_style=plot_menu_cursor_style,
        menu_highlight_style=plot_menu_style
    )
    selection_index = plot_menu.show()
    if selection_index == 0:
        print("No plot file chosen")
        sys.exit(1)
    else:
        return fits[selection_index - 1]


parser = argparse.ArgumentParser(description="Plotting tools for PartialWaveAnalysis")
parser.add_argument("path", help="path to fit directory")
parser.add_argument("-l", "--label", default="K_SK_S", help="LaTeX formated string of particles for x-axis label (default: \"K_SK_S\")")
parser.add_argument("-c", "--contours", action='store_true', help="plot contour lines")
parser.add_argument("-d", "--density", action='store_true', help="plot density")

if len(sys.argv) == 1: # if the user doesn't supply any arguments, print the help string and exit
    parser.print_help(sys.stderr)
    sys.exit(1)
args = parser.parse_args()

fit_results = plot_menu(Path(args.path)).resolve()
xlabel = rf"m$({args.label})$ GeV/$c^2$"

pdf = matplotlib.backends.backend_pdf.PdfPages(f"amps_{fit_results.stem.split('::')[0]}.pdf")

df = pd.read_csv(fit_results, delimiter='\t', index_col=False)
df.sort_values(['Bin', 'likelihood', 'total_intensity_err_AC'], ascending=[True, False, True], inplace=True)

def mask_first(x):
    result = np.zeros_like(x)
    result[0] = 1
    return result

mask = df.groupby(['Bin'])['Bin'].transform(mask_first).astype(bool)
df_filtered = df.loc[mask]
# print(df_filtered.columns)
# print(df_filtered.head())

bin_df = pd.read_csv(fit_results.parent / 'bin_info.txt', delimiter='\t')

amplitudes = [column[:-3] for column in df.columns.to_list() if column.endswith('_re')]
print(amplitudes)

wave_set = set([amp[:-1] for amp in amplitudes])
wave_dict = {'S': 0, 'P': 1, 'D': 2, 'F': 3, 'G': 4}
waves_sorted = sorted(list(wave_set), key=lambda wave: 100 * wave_dict[wave[0]] + (-1 if wave[-1] == '-' else 1) * int(wave[1]))
# kind of sneaky wave of sorting waves by L-letter and M-number without creating a whole new class
print(waves_sorted)

plt.rcParams["figure.figsize"] = (5 * len(waves_sorted), 10)
plt.rcParams["font.size"] = 24

for bin_n in range(len(bin_df)):
    print(f"Bin: {bin_n}    ", end='\r')
    fig, axes = plt.subplots(nrows=2, ncols=len(waves_sorted))
    for i, wave in enumerate(waves_sorted):
        amp_letter = wave[0] 
        amp_m = wave[1]
        if int(amp_m) > 0:
            amp_m_sign = wave[2]
        else:
            amp_m_sign = ""
        df_bin = df.loc[df['Bin'] == bin_n]
        df_conv = df_bin[df_bin['Convergence'] == 'C']
        wave_pos = wave + "+"
        if wave_pos in amplitudes:
            xmin = df_conv.loc[:, wave_pos + '_re'].min()
            xmax = df_conv.loc[:, wave_pos + '_re'].max()
            ymin = df_conv.loc[:, wave_pos + '_im'].min()
            ymax = df_conv.loc[:, wave_pos + '_im'].max()
            x_range = max(max(abs(xmax), abs(xmin)), 10) # 10 for real-wave plots
            y_range = max(max(abs(ymax), abs(ymin)), 10)
            try:
                # Credit to Flabetvibes at https://stackoverflow.com/a/30146280 for this code
                xx, yy = np.mgrid[-x_range:x_range:100j, -y_range:y_range:100j]
                positions = np.vstack([xx.ravel(), yy.ravel()])
                values = np.vstack([df_conv.loc[:, wave_pos + '_re'], df_conv.loc[:, wave_pos + '_im']])
                kernel = st.gaussian_kde(values)
                f = np.reshape(kernel(positions).T, xx.shape)
                if args.contours:
                    axes[0, i].contour(xx, yy, f, levels=5, colors='k')
                if args.density:
                    axes[0, i].contourf(xx, yy, f, levels=5, cmap='Blues')
            except:
                pass
            axes[0, i].scatter(df_conv.loc[:, wave_pos + '_re'], df_conv.loc[:, wave_pos + '_im'], color='k', marker=',')
            axes[0, i].scatter(df_filtered[wave_pos + '_re'].iloc[bin_n], df_filtered[wave_pos + '_im'].iloc[bin_n], color='r', marker=',')
            axes[0, i].set_xlabel("Re")
            axes[0, i].set_ylabel("Im")
            axes[0, i].set_title(rf"${amp_letter}^+_{{{amp_m_sign}{amp_m}}}$")
            axes[0, i].set_xlim(-x_range, x_range)
            axes[0, i].set_ylim(-y_range, y_range)
        wave_neg = wave + "-"
        if wave_neg in amplitudes:
            xmin = df_conv.loc[:, wave_neg + '_re'].min()
            xmax = df_conv.loc[:, wave_neg + '_re'].max()
            ymin = df_conv.loc[:, wave_neg + '_im'].min()
            ymax = df_conv.loc[:, wave_neg + '_im'].max()
            x_range = max(max(abs(xmax), abs(xmin)), 10) # 10 for real-wave plots
            y_range = max(max(abs(ymax), abs(ymin)), 10)
            try:
                xx, yy = np.mgrid[-x_range:x_range:100j, -y_range:y_range:100j]
                positions = np.vstack([xx.ravel(), yy.ravel()])
                values = np.vstack([df_conv.loc[:, wave_neg + '_re'], df_conv.loc[:, wave_neg + '_im']])
                kernel = st.gaussian_kde(values)
                f = np.reshape(kernel(positions).T, xx.shape)
                if args.contours:
                    axes[1, i].contour(xx, yy, f, levels=5, colors='k')
                if args.density:
                    axes[1, i].contourf(xx, yy, f, levels=5, cmap='Blues')
            except:
                pass
            axes[1, i].scatter(df_conv.loc[:, wave_neg + '_re'], df_conv.loc[:, wave_neg + '_im'], color='k', marker=',')
            axes[1, i].scatter(df_filtered[wave_neg + '_re'].iloc[bin_n], df_filtered[wave_neg + '_im'].iloc[bin_n], color='r', marker=',')
            axes[1, i].set_xlabel("Re")
            axes[1, i].set_ylabel("Im")
            axes[1, i].set_xlim(-x_range, x_range)
            axes[1, i].set_ylim(-y_range, y_range)
            axes[1, i].set_title(rf"${amp_letter}^-_{{{amp_m_sign}{amp_m}}}$")
    fig.suptitle(f"Bin {bin_n} Fit Amplitude Distributions")
    plt.tight_layout()
    pdf.savefig(fig, dpi=300)
    plt.close()

pdf.close()
print("Done!")
