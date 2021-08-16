#!/usr/bin/env python3

import argparse
import sys
from pathlib import Path

import matplotlib.backends.backend_pdf
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from simple_term_menu import TerminalMenu


def plot_menu(root):
    """

    :param root: 

    """
    plot_menu_title = "Select a fit_results.txt file:"
    fits = [p for p in root.glob("*::bootstrap.txt")]
    plot_menu_items = ["Cancel"] + [fit.name[:-15] for fit in fits]
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
        return root / Path(fits[selection_index - 1].name[:-15] + "::fit_results.txt")


parser = argparse.ArgumentParser(
    description="Plotting tools for PartialWaveAnalysis")
parser.add_argument("path", help="path to fit directory")
parser.add_argument("-u",
                    "--uncorrected",
                    action='store_true',
                    help="plot intensities without acceptance correction")
parser.add_argument("-c",
                    "--confidence",
                    action='store_true',
                    help="plot confidence intervals")
parser.add_argument("-p",
                    "--probability",
                    default=0.95,
                    help="set confidence level for intervals (default: 0.95)")
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

colors = [
    "tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple", "tab:brown",
    "tab:pink", "tab:olive", "tab:cyan"
]

pdf = matplotlib.backends.backend_pdf.PdfPages(
    f"bootstrap_{fit_results.stem.split('::')[0]}.pdf")

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

#############
df_filtered = df.loc[mask]
df_filtered['index'] = df_filtered[
    'Bin']  # do this so we don't lose the 'Bin' column
df_filtered.set_index('index', inplace=True)
df_bootstrap = pd.read_csv(fit_results.parent /
                           f"{fit_results.stem.split('::')[0]}::bootstrap.txt",
                           delimiter='\t',
                           index_col=False)
bin_df = pd.read_csv(fit_results.parent / 'bin_info.txt', delimiter='\t')
#############

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

wave_set = set([amp[:-1] for amp in amplitudes])
wave_dict = {'S': 0, 'P': 1, 'D': 2, 'F': 3, 'G': 4}
waves_sorted = sorted(list(wave_set),
                      key=lambda wave: 100 * wave_dict[wave[0]] +
                      (-1 if wave[-1] == '-' else 1) * int(wave[1]))
# kind of sneaky wave of sorting waves by L-letter and M-number without creating a whole new class

n_amps = len(waves_sorted) + 1  # plus 1 for the total plot

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
            amp_max = df_conv.loc[:, wave_pos + ac_tag].max()
            amp_min = df_conv.loc[:, wave_pos + ac_tag].min()
            correct_value = df_filtered[wave_pos + ac_tag][bin_n]
            entries, bins = np.histogram(df_conv.loc[:, wave_pos + ac_tag],
                                         bins=np.linspace(
                                             amp_min, amp_max, n_bins))
            bin_centers = np.array(
                [0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)])
            axes[0, i].bar(
                bin_centers,
                entries,
                width=bins[1] - bins[0],
                color='navy',
                label=
                f"$\mu$: {int(df_conv.loc[:, wave_pos + ac_tag].mean())}, $\sigma$: {int(df_conv.loc[:, wave_pos + ac_tag].std())}"
            )
            axes[0, i].legend()
            amp_letter = wave[0]
            amp_m = wave[1]
            if int(amp_m) > 0:
                amp_m_sign = wave[2]
            else:
                amp_m_sign = ""
            axes[0, i].set_title(f"${amp_letter}^+_{{{amp_m_sign}{amp_m}}}$")

        wave_neg = wave + "-"
        if wave_neg in amplitudes:
            print("-e", end='\t')
            amp_max = df_conv.loc[:, wave_neg + ac_tag].max()
            amp_min = df_conv.loc[:, wave_neg + ac_tag].min()
            correct_value = df_filtered[wave_neg + ac_tag][bin_n]
            entries, bins = np.histogram(df_conv.loc[:, wave_neg + ac_tag],
                                         bins=np.linspace(
                                             amp_min, amp_max, n_bins))
            bin_centers = np.array(
                [0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)])
            axes[1, i].bar(
                bin_centers,
                entries,
                width=bins[1] - bins[0],
                color='navy',
                label=
                f"$\mu$: {int(df_conv.loc[:, wave_neg + ac_tag].mean())}, $\sigma$: {int(df_conv.loc[:, wave_neg + ac_tag].std())}"
            )
            axes[1, i].legend()
            amp_letter = wave[0]
            amp_m = wave[1]
            if int(amp_m) > 0:
                amp_m_sign = wave[2]
            else:
                amp_m_sign = ""
            axes[1, i].set_title(f"${amp_letter}^-_{{{amp_m_sign}{amp_m}}}$")
        print()
    print("tot")
    amp_max = df_conv.loc[:, 'total_intensity' + ac_tag_total].max()
    amp_min = df_conv.loc[:, 'total_intensity' + ac_tag_total].min()
    correct_value = df_filtered['total_intensity' + ac_tag_total][bin_n]
    entries, bins = np.histogram(df_conv.loc[:,
                                             'total_intensity' + ac_tag_total],
                                 bins=np.linspace(amp_min, amp_max, n_bins))
    bin_centers = np.array(
        [0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)])
    axes[1, -1].bar(
        bin_centers,
        entries,
        width=bins[1] - bins[0],
        color='navy',
        label=
        f"$\mu$: {int(df_conv.loc[:, 'total_intensity' + ac_tag_total].mean())}, $\sigma$: {int(df_conv.loc[:, 'total_intensity' + ac_tag_total].std())}"
    )
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
        fit_value = df_filtered[amp + ac_tag][bin_n]
        df_bin = df_bootstrap.loc[df_bootstrap['Bin'] == bin_n]
        df_conv = df_bin[df_bin['Convergence'] == 'C']
        amp_bootstrap_errors.append(
            df_conv.loc[:, amp +
                        ac_tag].std())  # Using mean of distribution for now
        amp_bootstrap_CIL.append(2 * fit_value -
                                 np.quantile(df_conv.loc[:, amp + ac_tag], 1 -
                                             alpha / 2))
        amp_bootstrap_CIU.append(2 * fit_value -
                                 np.quantile(df_conv.loc[:, amp +
                                                         ac_tag], alpha / 2))

    df_filtered.loc[:, f"{amp}_bootstrap_err" + ac_tag] = amp_bootstrap_errors
    df_filtered.loc[:, f"{amp}_bootstrap_CIL" + ac_tag] = amp_bootstrap_CIL
    df_filtered.loc[:, f"{amp}_bootstrap_CIU" + ac_tag] = amp_bootstrap_CIU

amp_bootstrap_errors = []
amp_bootstrap_CIL = []
amp_bootstrap_CIU = []
alpha = 1 - args.probability
for bin_n in range(len(bin_df)):
    fit_value = df_filtered['total_intensity' + ac_tag_total][bin_n]
    df_bin = df_bootstrap.loc[df_bootstrap['Bin'] == bin_n]
    df_conv = df_bin[df_bin['Convergence'] == 'C']
    amp_bootstrap_errors.append(df_conv.loc[:, 'total_intensity' +
                                            ac_tag_total].std())
    amp_bootstrap_CIL.append(2 * fit_value - np.quantile(
        df_conv.loc[:, 'total_intensity' + ac_tag_total], 1 - alpha / 2))
    amp_bootstrap_CIU.append(2 * fit_value - np.quantile(
        df_conv.loc[:, 'total_intensity' + ac_tag_total], alpha / 2))
df_filtered.loc[:, f"total_bootstrap_err" + ac_tag_total] = amp_bootstrap_errors
df_filtered.loc[:, f"total_bootstrap_CIL" + ac_tag_total] = amp_bootstrap_CIL
df_filtered.loc[:, f"total_bootstrap_CIU" + ac_tag_total] = amp_bootstrap_CIU

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
        plt.errorbar(bin_df['mass'].iloc[df_filtered['Bin']],
                     df_filtered[wave_pos + ac_tag],
                     yerr=df_filtered[wave_pos + "_bootstrap_err" + ac_tag],
                     elinewidth=0.5,
                     fmt='o',
                     color='r',
                     label=r'$+\epsilon$')
        if args.confidence:
            plt.fill_between(bin_df['mass'].iloc[df_filtered['Bin']],
                             df_filtered[wave_pos + "_bootstrap_CIL" + ac_tag],
                             df_filtered[wave_pos + "_bootstrap_CIU" + ac_tag],
                             color='r',
                             alpha=0.1)
    else:
        print("", end='\t')

    wave_neg = wave + "-"
    if wave_neg in amplitudes:
        print("-e", end='\t')
        plt.errorbar(bin_df['mass'].iloc[df_filtered['Bin']],
                     df_filtered[wave_neg + ac_tag],
                     yerr=df_filtered[wave_neg + "_bootstrap_err" + ac_tag],
                     elinewidth=0.5,
                     fmt='o',
                     color='b',
                     label=r'$-\epsilon$')
        if args.confidence:
            plt.fill_between(bin_df['mass'],
                             df_filtered[wave_neg + "_bootstrap_CIL" + ac_tag],
                             df_filtered[wave_neg + "_bootstrap_CIU" + ac_tag],
                             color='b',
                             alpha=0.1)
    else:
        print("", end='\t')
    # Plot total
    print("tot")
    plt.errorbar(bin_df['mass'].iloc[df_filtered['Bin']],
                 df_filtered['total_intensity' + ac_tag_total],
                 yerr=df_filtered['total_bootstrap_err' + ac_tag_total],
                 elinewidth=0.5,
                 fmt='o',
                 color='k',
                 label='Total')
    if args.confidence:
        plt.fill_between(bin_df['mass'],
                         df_filtered["total_bootstrap_CIL" + ac_tag_total],
                         df_filtered["total_bootstrap_CIU" + ac_tag_total],
                         color='k',
                         alpha=0.1)
    plt.title(f"${amp_letter}_{{{amp_m_sign}{amp_m}}}$")
    plt.ylim(bottom=0)
    plt.xlim(bin_df['mass'].iloc[0] - 0.1, bin_df['mass'].iloc[-1] + 0.1)
    plt.ylabel("Intensity")
    plt.xlabel(xlabel)
    plt.legend()
    plt.tight_layout()
    pdf.savefig(fig, dpi=300)

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
        plt.errorbar(bin_df['mass'].iloc[df_filtered['Bin']],
                     df_filtered[wave_pos + ac_tag],
                     yerr=df_filtered[wave_pos + "_bootstrap_err" + ac_tag],
                     elinewidth=0.5,
                     fmt='o',
                     color=colors[i],
                     label=f"${amp_letter}_{{{amp_m_sign}{amp_m}}}$")
print("tot")
plt.errorbar(bin_df['mass'].iloc[df_filtered['Bin']],
             df_filtered['total_intensity' + ac_tag_total],
             yerr=df_filtered['total_bootstrap_err' + ac_tag_total],
             elinewidth=0.5,
             fmt='o',
             color='k',
             label="Total")
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
        plt.errorbar(bin_df['mass'].iloc[df_filtered['Bin']],
                     df_filtered[wave_neg + ac_tag],
                     yerr=df_filtered[wave_neg + "_bootstrap_err" + ac_tag],
                     elinewidth=0.5,
                     fmt='o',
                     color=colors[i],
                     label=f"${amp_letter}_{{{amp_m_sign}{amp_m}}}$")
print("tot")
plt.errorbar(bin_df['mass'].iloc[df_filtered['Bin']],
             df_filtered['total_intensity' + ac_tag_total],
             yerr=df_filtered['total_bootstrap_err' + ac_tag_total],
             elinewidth=0.5,
             fmt='o',
             color='k',
             label="Total")
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
        plt.errorbar(bin_df['mass'].iloc[df_filtered['Bin']],
                     df_filtered[wave_pos + ac_tag],
                     yerr=df_filtered[wave_pos + "_bootstrap_err" + ac_tag],
                     elinewidth=0.5,
                     fmt='o',
                     color=colors[i],
                     label=f"${amp_letter}^+_{{{amp_m_sign}{amp_m}}}$")
    else:
        print("", end='\t')

    wave_neg = wave + "-"
    if wave_neg in amplitudes:
        print("-e", end='\t')
        plt.errorbar(bin_df['mass'].iloc[df_filtered['Bin']],
                     df_filtered[wave_neg + ac_tag],
                     yerr=df_filtered[wave_neg + "_bootstrap_err" + ac_tag],
                     elinewidth=0.5,
                     fmt='s',
                     color=colors[i],
                     label=f"${amp_letter}^-_{{{amp_m_sign}{amp_m}}}$")
    else:
        print("", end='\t')
    print()
print("tot")
plt.errorbar(bin_df['mass'].iloc[df_filtered['Bin']],
             df_filtered['total_intensity' + ac_tag_total],
             yerr=df_filtered['total_bootstrap_err' + ac_tag_total],
             elinewidth=0.5,
             fmt='o',
             color='k',
             label="Total")
plt.xlim(bin_df['mass'].iloc[0] - 0.1, bin_df['mass'].iloc[-1] + 0.1)
plt.ylim(bottom=0)
plt.legend(loc="upper right")
plt.title("All Waves")
plt.ylabel("Intensity")
plt.xlabel(xlabel)
plt.tight_layout()
pdf.savefig(fig, dpi=300)

pdf.close()
