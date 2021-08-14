#!/usr/bin/env python3

import argparse
import sys
from pathlib import Path

import matplotlib.backends.backend_pdf
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from simple_term_menu import TerminalMenu


def plot_menu(root):
    """

    :param root: 

    """
    plot_menu_title = "Select a fit_results.txt file:"
    fits = [p for p in root.glob("*::fit_results.txt")]
    plot_menu_items = ["Cancel"] + [fit.name[:-17] for fit in fits]
    plot_menu_cursor = "> "
    plot_menu_cursor_style = ("fg_red", "bold")
    plot_menu_style = ("bg_black", "fg_green")
    plot_menu = TerminalMenu(menu_entries=plot_menu_items,
                             title=plot_menu_title,
                             menu_cursor=plot_menu_cursor,
                             menu_cursor_style=plot_menu_cursor_style,
                             menu_highlight_style=plot_menu_style)
    selection_index = plot_menu.show()
    if selection_index == 0:
        print("No plot file chosen")
        sys.exit(1)
    else:
        return fits[selection_index - 1]


parser = argparse.ArgumentParser(
    description="Plotting tools for PartialWaveAnalysis")
parser.add_argument("path", help="path to fit directory")
parser.add_argument("-u",
                    "--uncorrected",
                    action='store_true',
                    help="plot intensities without acceptance correction")
parser.add_argument(
    "-l",
    "--label",
    default="K_SK_S",
    help=
    "LaTeX formated string of particles for x-axis label (default: \"K_SK_S\")")

if len(
        sys.argv
) == 1:  # if the user doesn't supply any arguments, print the help string and exit
    parser.print_help(sys.stderr)
    sys.exit(1)
args = parser.parse_args()

fit_results = plot_menu(Path(args.path)).resolve()
xlabel = f"m$({args.label})$ GeV/$c^2$"

pdf = matplotlib.backends.backend_pdf.PdfPages(
    f"stats_{fit_results.stem.split('::')[0]}.pdf")

df = pd.read_csv(fit_results, delimiter='\t', index_col=False)
df.sort_values(['Bin', 'likelihood', 'total_intensity_err_AC'],
               ascending=[True, False, True],
               inplace=True)


def mask_first(x):
    """

    :param x: 

    """
    result = np.zeros_like(x)
    result[0] = 1
    return result


mask = df.groupby(['Bin'])['Bin'].transform(mask_first).astype(bool)
df_filtered = df.loc[mask]

bin_df = pd.read_csv(fit_results.parent / 'bin_info.txt', delimiter='\t')
if not args.uncorrected:
    ac_tag_total = "_AC"
    ac_tag = "_AC_INT"
    amplitudes = [
        column[:-len(ac_tag)]
        for column in df.columns.to_list()
        if column.endswith(ac_tag) and not "_err" in column
    ]
else:
    ac_tag_total = ""
    ac_tag = "_INT"
    amplitudes = [
        column[:-len(ac_tag)]
        for column in df.columns.to_list()
        if column.endswith(ac_tag) and not "_err" in column and
        not "_AC_" in column
    ]

plt.rcParams["figure.figsize"] = (30, 10)
plt.rcParams["font.size"] = 24

print("Plotting Separate Amplitudes")
# Positive
for i in range(len(amplitudes)):
    fig = plt.figure()
    all_runs_by_bin = [
        df[amplitudes[i] + ac_tag].loc[df['Bin'] == bin_n]
        for bin_n in bin_df['bin']
    ]
    plt.scatter(bin_df['mass'].iloc[df['Bin']],
                df[amplitudes[i] + ac_tag],
                marker='.',
                color='k',
                label="Fit Minima")
    plt.violinplot(all_runs_by_bin,
                   bin_df['mass'],
                   widths=bin_df['mass'].iloc[1] - bin_df['mass'].iloc[0],
                   showmeans=True,
                   showextrema=True,
                   showmedians=True)
    plt.scatter(bin_df['mass'].iloc[df_filtered['Bin']],
                df_filtered[amplitudes[i] + ac_tag],
                marker='o',
                color='r',
                label="Selected Minimum")
    plt.title(amplitudes[i].split("::")[-1])
    plt.xlim(bin_df['mass'].iloc[0] - 0.1, bin_df['mass'].iloc[-1] + 0.1)
    plt.ylim(bottom=0)
    plt.ylabel("Intensity")
    plt.xlabel(xlabel)
    plt.legend()
    plt.tight_layout()
    pdf.savefig(fig, dpi=300)

print("Plotting Total Intensity")
fig = plt.figure()
all_runs_by_bin = [
    df['total_intensity' + ac_tag_total].loc[df['Bin'] == bin_n]
    for bin_n in bin_df['bin']
]
plt.scatter(bin_df['mass'].iloc[df['Bin']],
            df['total_intensity' + ac_tag_total],
            marker='.',
            color='k',
            label="Fit Minima")
plt.violinplot(all_runs_by_bin,
               bin_df['mass'],
               widths=bin_df['mass'].iloc[1] - bin_df['mass'].iloc[0],
               showmeans=True,
               showextrema=True,
               showmedians=True)
plt.scatter(bin_df['mass'].iloc[df_filtered['Bin']],
            df_filtered['total_intensity' + ac_tag_total],
            marker='o',
            color='r',
            label="Selected Minimum")
plt.title('Total Intensity')
plt.xlim(bin_df['mass'].iloc[0] - 0.1, bin_df['mass'].iloc[-1] + 0.1)
plt.ylim(bottom=0)
plt.ylabel("Intensity")
plt.xlabel(xlabel)
plt.legend()
plt.tight_layout()
pdf.savefig(fig, dpi=300)

print("Plotting Normalized Log(Likelihood)")
fig = plt.figure()
all_runs_by_bin = [
    df['likelihood'].loc[df['Bin'] == bin_n].to_numpy() /
    df['total_intensity' + ac_tag_total].loc[df['Bin'] == bin_n].to_numpy()
    for bin_n in bin_df['bin']
]
plt.scatter(bin_df['mass'].iloc[df['Bin']],
            df['likelihood'].to_numpy() /
            df['total_intensity' + ac_tag_total].to_numpy(),
            marker='.',
            color='k',
            label="Fit Minima")
plt.violinplot(all_runs_by_bin,
               bin_df['mass'],
               widths=bin_df['mass'].iloc[1] - bin_df['mass'].iloc[0],
               showmeans=True,
               showextrema=True,
               showmedians=True)
plt.scatter(bin_df['mass'].iloc[df_filtered['Bin']],
            df_filtered['likelihood'].to_numpy() /
            df_filtered['total_intensity' + ac_tag_total].to_numpy(),
            marker='o',
            color='r',
            label="Selected Minimum")
plt.title('Likelihood')
plt.xlim(bin_df['mass'].iloc[0] - 0.1, bin_df['mass'].iloc[-1] + 0.1)
plt.ylabel("Log-Likelihood / Total Intensity")
plt.xlabel(xlabel)
plt.legend()
plt.tight_layout()
pdf.savefig(fig, dpi=300)

print("Plotting Total Intensity Error")
fig = plt.figure()
all_runs_by_bin = [
    df['total_intensity_err' + ac_tag_total].loc[df['Bin'] == bin_n]
    for bin_n in bin_df['bin']
]
plt.scatter(bin_df['mass'].iloc[df['Bin']],
            df['total_intensity_err' + ac_tag_total],
            marker='.',
            color='k',
            label="Fit Minima")
plt.violinplot(all_runs_by_bin,
               bin_df['mass'],
               widths=bin_df['mass'].iloc[1] - bin_df['mass'].iloc[0],
               showmeans=True,
               showextrema=True,
               showmedians=True)
plt.scatter(bin_df['mass'].iloc[df_filtered['Bin']],
            df_filtered['total_intensity_err' + ac_tag_total],
            marker='o',
            color='r',
            label="Selected Minimum")
plt.title('Error in Total Intensity')
plt.xlim(bin_df['mass'].iloc[0] - 0.1, bin_df['mass'].iloc[-1] + 0.1)
plt.ylabel("Intensity Error")
plt.xlabel(xlabel)
plt.legend()
plt.tight_layout()
pdf.savefig(fig, dpi=300)

pdf.close()
