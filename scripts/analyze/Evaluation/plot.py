

from matplotlib.legend_handler import HandlerTuple
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt 
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
from matplotlib.patches import Patch
from matplotlib.legend import Legend

import numpy as np
import pandas as pd
import seaborn as sns
from statsmodels.stats.proportion import proportion_confint
import ast, re
import matplotlib.cm as cm
from matplotlib.ticker import MaxNLocator, MultipleLocator, FuncFormatter

import warnings
warnings.filterwarnings("ignore", category=UserWarning) # ignore all UserWarning
warnings.filterwarnings("ignore", category=FutureWarning) # ignore FutureWarning




def plotCumulative(x, y_list, label_list, color_list, fontsize=20, filename='Cumulative.png'):
    plt.figure(figsize=(8, 6))
    X_ = np.linspace(x.min(), x.max(), 500)
    for (y, label, color) in zip(y_list, label_list, color_list):
        cubic_interpolation = interp1d(x, y, kind="linear")
        Y_ = cubic_interpolation(X_)
        plt.plot(X_, Y_, label=label, color=color)    
    plt.ylim(0, 1)
    plt.xlim(0, 10)
    plt.xticks(np.arange(0,11,1))
    plt.yticks(np.arange(0,1.1,0.1))
    plt.tick_params(axis='both', labelsize=fontsize)
    plt.xlabel("RMSD cutoff (Å)", fontsize=fontsize)
    plt.ylabel("Cumulative Frequency", fontsize=fontsize)
    #plt.legend(loc='lower right', fontsize=fontsize)  # Show labels
    plt.grid(True)
    #plt.show()
    plt.tight_layout()
    plt.savefig(filename,dpi=300)
    plt.clf()


def plot_bar(categories, values, colors, figname, lipidclass, yerrors, fontsize=12):
    """
    Function to generate a bar plot.
    
    Parameters:
    categories (list): List of category names.
    values (list): List of corresponding values.
    colors (list): List of colors for each bar.
    """
    lower_errors = [val_list[0] for val_list in yerrors]
    upper_errors = [val_list[1] for val_list in yerrors]
    # Create a bar plot
    yerr_list = [lower_errors, upper_errors]
    #print(f'Plotting test: yerr_list: {yerr_list}')

    lower_errors = np.array(values) - np.array(lower_errors)
    upper_errors = np.array(upper_errors) - np.array(values)
    yerr = [lower_errors, upper_errors]

    plt.bar(categories, values, color=colors, yerr=yerr, error_kw=dict(ecolor='black', lw=0.5, capsize=3, capthick=0.5))
    plt.yticks(np.arange(0,1.1,0.1))
    
    # Add labels and title
    plt.xlabel('Methods', fontsize=fontsize)
    plt.ylabel('Frequency of predictions below 2 Å', fontsize=fontsize)
    plt.title(lipidclass)
    
    # Display the plot
    #plt.show()
    plt.xticks(rotation=-45)
    plt.subplots_adjust(bottom=0.2)
    plt.tight_layout()

    plt.savefig(figname,dpi=300)
    plt.clf()




def _shades_from_base(base_color, n):
    """Generate n shades from base (dark) to light."""
    rgb = np.array(mcolors.to_rgb(base_color))
    white = np.ones(3)
    ts = np.linspace(1.0, 0.3, n)   # 1.0 ~ base (darkest), 0.3 ~ light
    return [tuple((1 - t) * white + t * rgb) for t in ts]


def plot_bar_new(categories, values, colors, figname, lipidclass="", fontsize=20, series_labels=None, percent=True):
    """
    Stacked bar plot from cumulative values.
    
    categories: list[str] of category names (one bar per category)
    values: list[list[float]] cumulative values per category, len == n_segments
    colors: list[str] base color per category
    series_labels: optional list[str] labels for the stacked segments
    """
    n_groups = len(categories)
    assert len(values) == n_groups == len(colors), "categories, values, colors length mismatch"
    n_segments = len(values[0])    
    plt.figure(figsize=(8, 6))

    cutoffs = [2, 2.5, 3]  # Example cutoffs for segments
    series_labels = [f"cutoff={c}" for c in cutoffs]

    if series_labels is None:
        series_labels = [f"Level {i+1}" for i in range(n_segments)]

    # Validate all categories have the same number of segments
    for v in values:
        assert len(v) == n_segments, "All inner lists in values must have equal length"

    x = np.arange(n_groups)

    # Patterns for the bar
    #patterns = ["*****", ".....", "/////"]   # match cutoffs
    patterns = ["***", "..", "//"]

    for i, (cat, cum_vals, base_color) in enumerate(zip(categories, values, colors)):
        cum_vals = np.asarray(cum_vals, dtype=float)

        # Ensure non-decreasing cumulative input
        if np.any(np.diff(cum_vals) < -1e-12):
            raise ValueError(f"Values for category '{cat}' must be non-decreasing (cumulative).")

        # If percent option is set, scale cumulative values to percent
        if percent:
            cum_vals = cum_vals * 100.0

        # Convert cumulative → segment sizes, total height = last cumulative value
        segs = np.diff(np.concatenate(([0.0], cum_vals)))

        shades = _shades_from_base(base_color, n_segments)
        bottom = 0.0
        for j, (seg, shade) in enumerate(zip(segs, shades)):
            if j == 0:
                plt.bar(
                    x[i], seg, bottom=bottom, width=0.6,
                    color=base_color,
                    #edgecolor=shade,  # edgecolor important for hatch visibility
                    edgecolor='black',
                    label=series_labels[j] if i == 0 else None,
                    #hatch=patterns[j]
                )
            else:
                plt.bar(
                    x[i], seg, bottom=bottom, width=0.6,
                    #color='white',
                    color=base_color,
                    #edgecolor=shade,  # edgecolor important for hatch visibility
                    edgecolor='black',
                    label=series_labels[j] if i == 0 else None,
                    hatch=patterns[j]
                )
            bottom += seg

    # Axes labels
    if percent:
        plt.ylim(0, 100.0)
        plt.yticks(np.arange(0, 110, 10), fontsize=fontsize)
        plt.ylabel('Success Rate (%)', fontsize=fontsize)
    else:
        plt.ylim(0, 1.0)
        plt.yticks(np.arange(0, 1.1, 0.1), fontsize=fontsize)
        plt.ylabel('Frequency of predictions below cutoff', fontsize=fontsize)
    plt.xticks(x, categories, rotation=-45, fontsize=fontsize, ha='left')
    plt.title(lipidclass, fontsize=fontsize)
    plt.subplots_adjust(bottom=0.2)

    legend_handles = []
    for i in range(len(cutoffs)-1, -1, -1):
        if i == 0:
            legend_handles.append(Patch(facecolor='lightgray', edgecolor="gray", label=f"cutoff={cutoffs[i]} Å"))
        else:
            legend_handles.append(Patch(facecolor='white', edgecolor="gray", hatch=patterns[i], label=f"cutoff={cutoffs[i]} Å"))

    plt.legend(handles=legend_handles, frameon=False, fontsize=fontsize)
    #plt.legend(legend_handles, legend_labels, handler_map={tuple: HandlerTuple(ndivide=len(colors))}, frameon=False, borderpad=0)
    plt.tight_layout()  # automatically adjusts padding
    plt.savefig(figname, dpi=300)
    plt.clf()





def plot_bar_MW(sucessrates_MW_dict, ci_intervals_MW_dict, num_datapoints_list, colors, fontsize=16, figname='Bar_MW.png', percent=True):

    """
    Function to generate a grouped bar plot for success rates by lipid molecular weight.
    """
    categories = list(sucessrates_MW_dict.keys())
    categories_num = [f'{cat}\n(n={num})' for cat, num in zip(categories, num_datapoints_list)]
    values_af = [sucessrates_MW_dict[cat][0] for cat in categories]
    values_chai = [sucessrates_MW_dict[cat][1] for cat in categories]
    values_vina = [sucessrates_MW_dict[cat][2] for cat in categories]
    values_rs = [sucessrates_MW_dict[cat][3] for cat in categories]
    values_dd = [sucessrates_MW_dict[cat][4] for cat in categories]

    if percent:
        values_af = [v * 100.0 for v in values_af]
        values_chai = [v * 100.0 for v in values_chai]
        values_vina = [v * 100.0 for v in values_vina]
        values_rs = [v * 100.0 for v in values_rs]
        values_dd = [v * 100.0 for v in values_dd]
    

    # plot
    if ci_intervals_MW_dict is not None:
        values = [values_af, values_chai, values_vina, values_rs, values_dd]
        ci_intervals_af = [ci_intervals_MW_dict[cat][0] for cat in categories]
        ci_intervals_chai = [ci_intervals_MW_dict[cat][1] for cat in categories]
        ci_intervals_vina = [ci_intervals_MW_dict[cat][2] for cat in categories]
        ci_intervals_rs = [ci_intervals_MW_dict[cat][3] for cat in categories]
        ci_intervals_dd = [ci_intervals_MW_dict[cat][4] for cat in categories]
        if percent:
            ci_intervals_af = [v * 100.0 for v in ci_intervals_af]
            ci_intervals_chai = [v * 100.0 for v in ci_intervals_chai]
            ci_intervals_vina = [v * 100.0 for v in ci_intervals_vina]
            ci_intervals_rs = [v * 100.0 for v in ci_intervals_rs]
            ci_intervals_dd = [v * 100.0 for v in ci_intervals_dd]
        ci_intervals = [ci_intervals_af, ci_intervals_chai, ci_intervals_vina, ci_intervals_rs, ci_intervals_dd]
        x = np.arange(len(categories))  # the label locations
        width = 0.15  # the width of the bars
        fig, ax = plt.subplots(figsize=(12, 6))
        labels = ['AlphaFold 3', 'Chai-1', 'Autodock Vina', 'RoseTTAFoldAA', 'DiffDock-L']
        locations = [x - 2*width, x - width, x, x + width, x + 2*width]
        for i, (values, ci_intervals, label, location) in enumerate(zip(values, ci_intervals, labels, locations)):
            lower_errors = [val_list[0] for val_list in ci_intervals]
            upper_errors = [val_list[1] for val_list in ci_intervals]
            lows = np.array(values) - np.array(lower_errors)
            uppers = np.array(upper_errors) - np.array(values)
            yerr = [lows, uppers]
            rects = ax.bar(location, values, width, label=label, color=colors[i], yerr=yerr, error_kw=dict(ecolor='black', lw=0.5, capsize=3, capthick=0.5))

    else:
        x = np.arange(len(categories))  # the label locations
        width = 0.15  # the width of the bars

        fig, ax = plt.subplots(figsize=(12, 6))
        rects1 = ax.bar(x - 2*width, values_af, width, label='AlphaFold 3', color=colors[0])
        rects2 = ax.bar(x - width, values_chai, width, label='Chai-1', color=colors[1])
        rects3 = ax.bar(x, values_vina, width, label='Autodock Vina', color=colors[2])
        rects4 = ax.bar(x + width, values_rs, width, label='RoseTTAFoldAA', color=colors[3])
        rects5 = ax.bar(x + 2*width, values_dd, width, label='DiffDock-L', color=colors[4])

    # Add  text for labels, title and custom x-axis tick labels, etc.
    ax.tick_params(axis='both', labelsize=fontsize)
    ax.set_xlabel('Lipid Molecular Weight (Da)', fontsize=fontsize)
    #ax.set_title('Success Rates by Lipid Molecular Weight', fontsize=fontsize)
    ax.set_xticks(x) 
    ax.set_xticklabels(categories_num, rotation=-45, fontsize=fontsize, ha='left')
    if percent:
        ax.set_yticks(np.arange(0, 110, 10))
        ax.set_ylabel('Success Rate (%)', fontsize=fontsize)
        ax.set_yticklabels([f"{tick:.0f}" for tick in np.arange(0, 110, 10)], fontsize=fontsize)
    else:
        ax.set_yticks(np.arange(0, 1.1, 0.1))
        ax.set_ylabel('Frequency of predictions below 2 Å', fontsize=fontsize)
        ax.set_yticklabels([f"{tick:.1f}" for tick in np.arange(0, 1.1, 0.1)], fontsize=fontsize)
    #ax.legend(loc='upper left', fontsize=fontsize)

    fig.tight_layout()  
    plt.savefig(figname,dpi=300)
    plt.clf()




def plot_bar_lipid(sucessrates_lipid_dict, ci_intervals_lipid_dict, num_datapoints_list, colors, fontsize=16, figname='Bar_LipidClass.png', percent=True):
    """
    Function to generate a bar plot for success rates by lipid class.
    """
    categories = list(sucessrates_lipid_dict.keys())
    categories_num = [f'{cat}\n(n={num})' for cat, num in zip(categories, num_datapoints_list)]
    def scale_percent(seq):
        # Handles both flat lists and numpy arrays, and ignores non-numeric types
        return [float(v) * 100.0 if isinstance(v, (int, float, np.integer, np.floating)) else v for v in seq]

    values_af = scale_percent([sucessrates_lipid_dict[cat][0] for cat in categories]) if percent else [sucessrates_lipid_dict[cat][0] for cat in categories]
    values_chai = scale_percent([sucessrates_lipid_dict[cat][1] for cat in categories]) if percent else [sucessrates_lipid_dict[cat][1] for cat in categories]
    values_vina = scale_percent([sucessrates_lipid_dict[cat][2] for cat in categories]) if percent else [sucessrates_lipid_dict[cat][2] for cat in categories]
    values_rs = scale_percent([sucessrates_lipid_dict[cat][3] for cat in categories]) if percent else [sucessrates_lipid_dict[cat][3] for cat in categories]
    values_dd = scale_percent([sucessrates_lipid_dict[cat][4] for cat in categories]) if percent else [sucessrates_lipid_dict[cat][4] for cat in categories]
    values = [values_af, values_chai, values_vina, values_rs, values_dd]

    ci_intervals_af = scale_percent([ci_intervals_lipid_dict[cat][0] for cat in categories]) if percent else [ci_intervals_lipid_dict[cat][0] for cat in categories]
    ci_intervals_chai = scale_percent([ci_intervals_lipid_dict[cat][1] for cat in categories]) if percent else [ci_intervals_lipid_dict[cat][1] for cat in categories]
    ci_intervals_vina = scale_percent([ci_intervals_lipid_dict[cat][2] for cat in categories]) if percent else [ci_intervals_lipid_dict[cat][2] for cat in categories]
    ci_intervals_rs = scale_percent([ci_intervals_lipid_dict[cat][3] for cat in categories]) if percent else [ci_intervals_lipid_dict[cat][3] for cat in categories]
    ci_intervals_dd = scale_percent([ci_intervals_lipid_dict[cat][4] for cat in categories]) if percent else [ci_intervals_lipid_dict[cat][4] for cat in categories]
    ci_intervals = [ci_intervals_af, ci_intervals_chai, ci_intervals_vina, ci_intervals_rs, ci_intervals_dd]

    # plot
    x = np.arange(len(categories))  # the label locations
    width = 0.15  # the width of the bars
    fig, ax = plt.subplots(figsize=(12, 6))
    labels = ['AlphaFold 3', 'Chai-1', 'Autodock Vina', 'RoseTTAFoldAA', 'DiffDock-L']
    locations = [x - 2*width, x - width, x, x + width, x + 2*width]
    for i, (values, ci_intervals, label, location) in enumerate(zip(values, ci_intervals, labels, locations)):
        #print(f'Plotting {label}, values: {values}, ci_intervals: {ci_intervals}')
        # lower_errors = [val_list[0] for val_list in ci_intervals]
        # upper_errors = [val_list[1] for val_list in ci_intervals]
        # lows = np.array(values) - np.array(lower_errors)
        # uppers = np.array(upper_errors) - np.array(values)
        # yerr = [lows, uppers]
        # rects = ax.bar(location, values, width, label=label, color=colors[i], yerr=yerr, error_kw=dict(ecolor='black', lw=0.5, capsize=3, capthick=0.5))
        rects = ax.bar(location, values, width, label=label, color=colors[i])


    # Add text for labels, title and custom x-axis tick labels, etc.
    ax.tick_params(axis='both', labelsize=fontsize)
    ax.set_xlabel('Lipid Class', fontsize=fontsize)
    #ax.set_title('Success Rates by Lipid Class', fontsize=fontsize)
    ax.set_xticks(x)
    ax.set_xticklabels(categories_num, rotation=-45, fontsize=fontsize, ha='left')
    if percent:
        ax.set_yticks(np.arange(0, 110, 10))
        ax.set_ylabel('Success Rate (%)', fontsize=fontsize)
        ax.set_yticklabels([f"{tick:.0f}" for tick in np.arange(0, 110, 10)], fontsize=fontsize)
    else:
        ax.set_yticks(np.arange(0, 1.1, 0.1))
        ax.set_ylabel('Frequency of predictions below 2 Å', fontsize=fontsize)
        ax.set_yticklabels([f"{tick:.1f}" for tick in np.arange(0, 1.1, 0.1)], fontsize=fontsize)
    #ax.legend(loc='center left', bbox_to_anchor=(1.02, 0.5), fontsize=fontsize, frameon=False)
    fig.tight_layout()
    plt.savefig(figname, dpi=300, bbox_inches='tight')
    plt.clf()



def plot_scatter_protein_errors(df, prefix_list, fig_path, cutoff="2023-05-20", fontsize=16):
    # Ensure df['date'] is datetime
    cutoff = pd.to_datetime(cutoff).tz_localize("UTC")
    df['PDB_release_date'] = pd.to_datetime(df['PDB_release_date'])

    colors_dict = {'AF': "#318df0", 'CHAI': "#f37019", 'RS': "#f3342d"}
     
    protein_error_metrics = ['lddt', 'tm', 'RMSD']
    for metric in protein_error_metrics:
        for prefix in prefix_list:
            x_array = df[f'{prefix}_protein_{metric}_openstc'].tolist()
            y_array = df[f'{prefix}_lipid_RMSD_spy'].tolist()
            color = colors_dict[prefix]

            with plt.rc_context({
                'font.size': fontsize,
                'axes.labelsize': fontsize,
                'axes.titlesize': fontsize + 2,
                'xtick.labelsize': fontsize,
                'ytick.labelsize': fontsize,
                'legend.fontsize': fontsize
            }):
                # Colors based on date condition (UTC-aware)
                #colors = np.where(df['PDB_release_date'] > cutoff, '#FFB347', '#AEC6CF')
                plt.figure(figsize=(8, 6))

                # Identify points with pTM > 0.8
                if f'{prefix}_ptm' in df.columns and metric == 'tm':
                    mask_highptm = df[f'{prefix}_ptm'] > 0.8
                    mask_lowptm = df[f'{prefix}_ptm'] <= 0.8

                    x_highptm = df.loc[mask_highptm, f'{prefix}_protein_{metric}_openstc']
                    y_highptm = df.loc[mask_highptm, f'{prefix}_lipid_RMSD_spy']
                    x_lowptm = df.loc[mask_lowptm, f'{prefix}_protein_{metric}_openstc']
                    y_lowptm = df.loc[mask_lowptm, f'{prefix}_lipid_RMSD_spy']

                    plt.scatter(x_highptm, y_highptm, s=10, c=color, label='pTM > 0.8')  # blue
                    plt.scatter(x_lowptm, y_lowptm, s=10, c=color, alpha=0.3, label='pTM ≤ 0.8')  # pink

                    patch_highptm = mpatches.Patch(color=color, label='pTM > 0.8')
                    patch_lowptm = mpatches.Patch(color=color, alpha=0.3, label='pTM ≤ 0.8')
                    plt.legend(handles=[patch_highptm, patch_lowptm], fontsize=fontsize, loc='upper left')
                else:
                    plt.scatter(x_array, y_array, s=10, c=color)
                    patch_all = mpatches.Patch(color=color, label='All')
                    #patch_new = mpatches.Patch(color='#FFB347', label='Untrained data')
                    #plt.legend(handles=[patch_all, patch_new], fontsize=fontsize)
                    plt.legend(handles=[patch_all], fontsize=fontsize, loc='upper left')

                #plt.title(f'Lipid RMSD vs Protein Structure Prediction Errors ({metric.upper()}, {prefix})', fontsize=fontsize + 2)
                plt.xlabel(f'Protein {metric.upper()}-score', fontsize=fontsize)
                plt.ylabel('Lipid RMSD (Å)', fontsize=fontsize)
                plt.xticks(fontsize=fontsize)
                plt.yticks(fontsize=fontsize)
                plt.axhline(y=2, color='red', linestyle='--')

                plt.tight_layout()
                plt.savefig(f'{fig_path}/Scatter_lipid_vs_protein_{metric.upper()}_{prefix}.png', dpi=300)
                plt.clf()




def plot_scoring_power(df, fig_path, fontsize=20):
    method_scoring = {"AF":"iptm", "CHAI":"iptm", "VINA":"Affinity", "RS":"pae_inter", "DD":"Confidence_Score"}

    for prefix in method_scoring.keys():
        with plt.rc_context({
            'font.size': fontsize,
            'axes.labelsize': fontsize,
            'axes.titlesize': fontsize + 2,
            'xtick.labelsize': fontsize,
            'ytick.labelsize': fontsize,
            'legend.fontsize': fontsize
        }):
            # ----------------------------------- plot scatter plot: -----------------------------------
            y_array = df[f'{prefix}_lipid_RMSD_spy'].tolist()
            x_array = df[f'{prefix}_{method_scoring[prefix]}'].tolist()
        
            plt.figure(figsize=(8, 6))
            plt.scatter(x_array, y_array, s=5)
            plt.title(f'{method_scoring[prefix]} vs Lipid RMSD (Å)', fontsize=fontsize + 2)
            plt.ylabel('Lipid RMSD (Å)', fontsize=fontsize)
            plt.xlabel(f'{method_scoring[prefix]} Score', fontsize=fontsize)
            #plt.grid(True)
            plt.ylim(bottom=0)
            if prefix == 'DD':
                plt.axvline(x=0, color='gray', linestyle='--')
                plt.axvline(x=-1.5, color='gray', linestyle='--')
            elif prefix == 'RS':
                plt.axvline(x=10, color='gray', linestyle='--')
            elif prefix == 'AF' or prefix == 'CHAI':
                plt.xlim(0, 1)
                plt.axvline(x=0.6, color='gray', linestyle='--')
                plt.axvline(x=0.8, color='gray', linestyle='--')

            plt.tight_layout()
            plt.savefig(f'{fig_path}/ScoringPower_{prefix}.png', dpi=300)
            plt.clf()

        # ----------------------------------- plot boxplot based on differnt bins of x_array -----------------------------------
        if prefix == 'DD':
            bins = [-np.inf, -6, -4.5, -3, -1.5, 0, np.inf]
            #labels = [f"({bins[i]} to {bins[i+1]}]" for i in range(len(bins)-1)]
            labels = ["< -6", "-6 ~ -4.5", "-4.5 ~ -3", "-3 ~ -1.5", "-1.5 ~ 0", "> 0"]           
        elif prefix == 'RS':
            bins = [0, 5, 10, 15, 20, 25, 30]
            labels = ["0 ~ 5", "5 ~ 10", "10 ~ 15", "15 ~ 20", "20 ~ 25", "25 ~ 30"]
            #labels = [f"({bins[i]} to {bins[i+1]}]" for i in range(len(bins)-1)]
        elif prefix == 'AF' or prefix == 'CHAI':
            bins = [0, 0.2, 0.4, 0.6, 0.8, 1.0]
            labels = ["0 ~ 0.2", "0.2 ~ 0.4", "0.4 ~ 0.6", "0.6 ~ 0.8", "0.8 ~ 1.0"]
            #labels = [f"({bins[i]} to {bins[i+1]}]" for i in range(len(bins)-1)]
        elif prefix == 'VINA':
            bins = [-20, -10, 0, 10, 20, 30, 40, 50, 60, 70, 80]
            labels = ["-20 ~ -10", "-10 ~ 0", "0 ~ 10", "10 ~ 20", "20 ~ 30", "30 ~ 40", "40 ~ 50", "50 ~ 60", "60 ~ 70", "70 ~ 80"]
            #labels = [f"({bins[i]} to {bins[i+1]}]" for i in range(len(bins)-1)]

        # Put into DataFrame
        data = pd.DataFrame({"x": x_array, "y": y_array})
        #print(f'x_array: {x_array}, y_array: {y_array}')
        data["interval"] = pd.cut(data["x"], bins=bins, labels=labels)
        data = data.dropna(subset=["interval"])
        data["interval"] = data["interval"].cat.remove_unused_categories()
        #print(f'data: {data}')


        # Create boxplot
        plt.figure(figsize=(8,5))
        #data.boxplot(column="y", by="interval", grid=False)
        sns.violinplot(
        x="interval",
        y="y",
        data=data,
        palette="pastel")

        #plt.title(f'{method_scoring[prefix]} vs Lipid RMSD (Å)')
        plt.xlabel(f'{method_scoring[prefix].replace("_", " ")}', fontsize=fontsize)
        plt.ylabel('Lipid RMSD (Å)', fontsize=fontsize)
        plt.xticks(rotation=-45, fontsize=fontsize)
        plt.yticks(fontsize=fontsize)
        plt.axhline(y=2, color='red', linestyle='--')
        plt.tight_layout() 
        plt.savefig(f'{fig_path}/ScoringPower_{prefix}_Violin.png', dpi=300)
        plt.clf()


def plot_pie_family(df, fig_name):

    #df["protein_Pfam_first"] = df["protein_Pfam"].apply(to_list).apply(lambda x: x[0] if isinstance(x, list) and len(x) > 0 else None)
    df["protein_Pfam_first"] = (
    df["protein_Pfam"]
    .apply(to_list)  # keeping your existing to_list function
    .apply(lambda x: remove_brackets(x[0]) if isinstance(x, list) and len(x) > 0 else None)
)
    #df["lipid_Lipidmaps_categories_first"] = df["lipid_Lipidmaps_categories"].apply(to_list).apply(lambda x: x[0] if isinstance(x, list) and len(x) > 0 else None)
    df["lipid_Lipidmaps_categories_first"] = (
    df["lipid_Lipidmaps_categories"]
    .apply(to_list)
    .apply(lambda x: '\n&\n'.join(x[:2]) if isinstance(x, list) and len(x) > 0 else None))

    # Count occurrences
    counts = df.groupby(["lipid_Lipidmaps_categories_first", "protein_Pfam_first"]).size().reset_index(name="count")

    # Outer ring (Protein families)
    outer_counts = counts.groupby("lipid_Lipidmaps_categories_first")["count"].sum()

    #fig, ax = plt.subplots(figsize=(8,6))
    fig = plt.figure(figsize=(9,6))
    ax = fig.add_axes([0.0, 0.1, 0.6, 0.8])  # left, bottom, width, height
    ax.axis('equal')  # keep pie circular
    ax.axis('off')    # hide axes


    # setup colors
    # num_inner = len(df["protein_Pfam_first"].unique())
    # cmap = cm.get_cmap("Pastel1")  # "Pastel2" also works
    # inner_colors = cmap(np.linspace(0, 1, num_inner))
    num_wedges = len(counts)
    #cmap = cm.get_cmap("Pastel1")
    cmaps = [cm.tab20, cm.tab20b, cm.tab20c]

    #print(f'cmap colors: {cmap}')
    #inner_colors = cmap(np.linspace(0, 1, num_wedges)) 
    #inner_colors = pastel_colors(num_wedges)
    inner_colors = [cmaps[i // 20](i % 20) for i in range(num_wedges)]

    # Outer pie
    patches, texts = ax.pie(
        outer_counts.values,
        labels=outer_counts.index,
        radius=1,
        wedgeprops=dict(width=0.3, edgecolor='w'), textprops={'fontsize': 8}, labeldistance=1.1
    )

    for i, text in enumerate(texts):
        if text.get_text() == "Fatty Acyl (FA)\n&\nGlycerolipid (GL)":
            #text.set_rotation(20)  
            pass
        elif text.get_text() == "Glycerolipid (GL)":
            #text.set_rotation(45)
            x, y = text.get_position()
            # Move lower by subtracting from y
            text.set_position((x + 0.1 , y - 0.05))  # adjust 0.05 as needed
  

    # Inner pie
    wedges, _ = ax.pie(
        counts["count"].values,
        radius=0.7,
        wedgeprops=dict(width=0.3, edgecolor='w'),
        textprops={'fontsize': 3},
        colors=inner_colors
    )
    inner_labels = counts["protein_Pfam_first"].tolist()
    ax.legend(wedges, inner_labels, title="Protein Pfam Categories", loc="center left", bbox_to_anchor=(1.2, 0.5), fontsize=8)


    plt.tight_layout()
    plt.savefig(fig_name, dpi=300)
    plt.close(fig)


def pastel_colors(n):
    np.random.seed(0)
    colors = np.random.uniform(0.6, 1.0, size=(n, 3))  # RGB values in [0.6,1.0] → soft colors
    return colors


def to_list(x):
    if isinstance(x, list):
        return x
    if pd.isna(x):
        return []
    if isinstance(x, str):
        s = x.strip()
        # Try parsing Python-list-like strings
        if s.startswith("[") and s.endswith("]"):
            try:
                v = ast.literal_eval(s)
                if isinstance(v, list):
                    return v
            except (SyntaxError, ValueError):
                pass
        # Fallback: split on commas/semicolons
        return [i.strip() for i in re.split(r"[;,]", s) if i.strip()]
    return []


def remove_brackets(s):
    if isinstance(s, str):
        return re.sub(r"\[.*?\]", "", s).strip()
    return s





def plot_scatter_protein_errors_seaborn(df, prefix_list, fig_path, cutoff="2023-05-20", fontsize=16):
    """
    Seaborn-based version of plot_scatter_protein_errors.

    Behaves similarly (same inputs and outputs). Uses seaborn.scatterplot to
    produce plots, preserves the pTM-based splitting (when available) and
    saves PNG files to the same naming scheme as the original function.
    """
    # Ensure datetime handling is identical
    cutoff = pd.to_datetime(cutoff).tz_localize("UTC")
    df['PDB_release_date'] = pd.to_datetime(df['PDB_release_date'])

    colors_dict = {'AF': "#2e8df2", 'CHAI': "#f27420", 'RS': "#f6aca8"}
    protein_error_metrics = ['lddt', 'tm', 'RMSD']

    for metric in protein_error_metrics:
        for prefix in prefix_list:
            x_col = f'{prefix}_protein_{metric}_openstc'
            y_col = f'{prefix}_lipid_RMSD_spy'
            color = colors_dict.get(prefix, '#999999')

            with plt.rc_context({
                'font.size': fontsize,
                'axes.labelsize': fontsize,
                'axes.titlesize': fontsize + 2,
                'xtick.labelsize': fontsize,
                'ytick.labelsize': fontsize,
                'legend.fontsize': fontsize
            }):
                plt.figure(figsize=(8, 6))

                # If pTM exists and metric is 'tm', split by pTM threshold
                if f'{prefix}_ptm' in df.columns and metric == 'tm':
                    mask_high = df[f'{prefix}_ptm'] > 0.8
                    df_high = df.loc[mask_high, [x_col, y_col]].dropna()
                    df_low = df.loc[~mask_high, [x_col, y_col]].dropna()

                    if not df_high.empty:
                        sns.set_style("darkgrid")
                        sns.scatterplot(data=df_high, x=x_col, y=y_col, color=color, s=40, label='pTM > 0.8')
                        
                    if not df_low.empty:
                        sns.set_style("darkgrid")
                        sns.scatterplot(data=df_low, x=x_col, y=y_col, color=color, alpha=0.3, s=40, label='pTM ≤ 0.8')
                        


                    # Create legend patches similar to original
                    patch_high = mpatches.Patch(color=color, label='pTM > 0.8')
                    patch_low = mpatches.Patch(color=color, alpha=0.3, label='pTM ≤ 0.8')
                    plt.legend(handles=[patch_high, patch_low], fontsize=fontsize, loc='upper left')
                else:
                    # Use seaborn for the whole set
                    tmp = df[[x_col, y_col]].dropna()
                    if not tmp.empty:
                        sns.set_style("darkgrid")
                        sns.scatterplot(data=tmp, x=x_col, y=y_col, color=color, s=40)
                    patch_all = mpatches.Patch(color=color, label='All')
                    plt.legend(handles=[patch_all], fontsize=fontsize, loc='upper left')
                    


                plt.xlabel(f'Protein {metric.upper()}-score', fontsize=fontsize)
                plt.ylabel('Lipid RMSD (Å)', fontsize=fontsize)
                plt.xticks(fontsize=fontsize)
                plt.yticks(fontsize=fontsize)
                plt.gca().xaxis.set_major_formatter(FuncFormatter(lambda x, pos: f"{x:.1f}"))
                plt.axhline(y=2, color='red', linestyle='--')

                plt.tight_layout()
                plt.savefig(f'{fig_path}/Scatter_lipid_vs_protein_{metric.upper()}_{prefix}_seaborn.png', dpi=300)
                plt.clf()
    

def plot_bar_protein_errors_seaborn(df, prefix_list, fig_path, fontsize=16):
    """ 
    Plot the TM score distributions for different methods using seaborn.
    x-axis: methods (mapped from prefix)
    y-axis: TM score distributions
    """
    
    # Mapping from prefix to full method names
    prefix_to_name = {'AF': 'AlphaFold 3', 'CHAI': 'Chai-1', 'RS': 'RoseTTAFoldAA'}
    colors_dict = {'AF': '#a1c9f4', 'CHAI': '#ffb482', 'RS': '#ff9f9b'}

    # Collect TM score data from each prefix column
    data_list = []
    for prefix in prefix_list:
        col = f"{prefix}_protein_tm_openstc"
        if col in df.columns:
            temp_df = pd.DataFrame({
                'Method': prefix_to_name.get(prefix, prefix),
                'TM_Score': df[col]
            })
            data_list.append(temp_df)
    
    # Combine all data into one long-form dataframe for seaborn
    combined_df = pd.concat(data_list, ignore_index=True)
    
    # Create color palette that matches each method
    palette = {prefix_to_name[p]: colors_dict[p] for p in prefix_list if p in colors_dict}


    with plt.rc_context({
        'font.size': fontsize,
        'axes.labelsize': fontsize,
        'axes.titlesize': fontsize + 2,
        'xtick.labelsize': fontsize,
        'ytick.labelsize': fontsize,
        'legend.fontsize': fontsize
    }):
        plt.figure(figsize=(8, 6))
        sns.set_style("darkgrid")

        # Create boxplot to show distribution
        sns.boxplot(
            data=combined_df,
            y='Method',
            x='TM_Score',
            palette=palette
        )
        
        plt.tick_params(axis='y', which='major', pad=16)
        plt.ylabel('Method', fontsize=fontsize)
        plt.xlabel('Protein TM Score', fontsize=fontsize)

        plt.xlim(0, 1)
        plt.tight_layout()
        plt.savefig(f'{fig_path}/Boxplot_protein_TM_seaborn.png', dpi=300)
        plt.clf()
