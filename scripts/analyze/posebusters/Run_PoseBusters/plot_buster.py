import matplotlib.pyplot as plt 
import pandas as pd
import os
import seaborn as sns




def Calculate_stats(df):
    stats_dict = {}

    # Columns to analyze
    cols = ["File_Loads", "sanitization","inchi_convertible","all_atoms_connected","molecular_formula","molecular_bonds","double_bond_stereochemistry","tetrahedral_chirality","bond_lengths","bond_angles","internal_steric_clash","aromatic_ring_flatness","non-aromatic_ring_non-flatness","double_bond_flatness","internal_energy","protein-ligand_maximum_distance","minimum_distance_to_protein","minimum_distance_to_organic_cofactors","minimum_distance_to_inorganic_cofactors","minimum_distance_to_waters","volume_overlap_with_protein","volume_overlap_with_organic_cofactors","volume_overlap_with_inorganic_cofactors","volume_overlap_with_waters"]
    # Add column to check if the columns are all true
    df["File_Loads"] = (df["mol_pred_loaded"] & df["mol_true_loaded"] & df["mol_cond_loaded"])
    for col in cols:
        percent_true = df[col].fillna(0).mean() * 100
        stats_dict[col] = percent_true

    stats_categories_dict = {}
    cols_chem_valid = ["File_Loads", "sanitization","inchi_convertible","all_atoms_connected","molecular_formula","molecular_bonds"]
    cols_intra_valid = ["double_bond_stereochemistry","tetrahedral_chirality","bond_lengths","bond_angles","internal_steric_clash","aromatic_ring_flatness","non-aromatic_ring_non-flatness","double_bond_flatness","internal_energy"]
    cols_inter_valid = ["protein-ligand_maximum_distance","minimum_distance_to_protein","minimum_distance_to_organic_cofactors","minimum_distance_to_inorganic_cofactors","minimum_distance_to_waters","volume_overlap_with_protein","volume_overlap_with_organic_cofactors","volume_overlap_with_inorganic_cofactors","volume_overlap_with_waters"]
    for category in ["Chemical_Validity", "Intramolecular_Validity", "Intermolecular_Validity"]:
        if category == "Chemical_Validity":
            # calculate the percentage of rows where all values in the columns are True
            pct = (df[cols_chem_valid].fillna(0).all(axis=1).mean()) * 100
        elif category == "Intramolecular_Validity":
            pct = (df[cols_intra_valid].fillna(0).all(axis=1).mean()) * 100
        else:
            category == "Intermolecular_Validity"
            pct = (df[cols_inter_valid].fillna(0).all(axis=1).mean()) * 100
        stats_categories_dict[category] = pct           
            
    return (stats_dict, stats_categories_dict)




def PlotBuster(list_of_stats_dicts, prefix_list, colors=None):
    # items = list(stats_dict.keys())
    # percentages = list(stats_dict.values())

    assert len(list_of_stats_dicts) == len(prefix_list), "Length of stats_dicts and prefix_list must match"

    # Default color palette if not provided
    if colors is None:
        colors = sns.color_palette("Set2", len(list_of_stats_dicts))

    # Make grid style
    sns.set_style("whitegrid", {
        'grid.color': '0.85',
        'grid.linestyle': '-',
        'axes.edgecolor': '0.3',
        'axes.linewidth': 1,
        'xtick.bottom': True,
        'ytick.left': True})

    plt.figure(figsize=(12, 7)) 
    plt.rcParams.update({'font.size': 16})

    # Convert each dict into a DataFrame with a 'Group' column
    all_df = []
    for stats_dict, group_name in zip(list_of_stats_dicts, prefix_list):
        df = pd.DataFrame(list(stats_dict.items()), columns=["Item", "Percentage"])
        df["Group"] = group_name
        all_df.append(df)

    # Concatenate all DataFrames
    plot_df = pd.concat(all_df, ignore_index=True)
    rename_map = {
    "AF": "AlphaFold3",
    "CHAI": "Chai-1",
    "autodock": "AutoDock Vina",
    "rfaa": "RoseTTAFoldAA",
    "DD": "DiffDock-L" 
    }

    plot_df["Group"] = plot_df["Group"].replace(rename_map)

    # Items to keep for plotting
    selected_items = ["File_Loads", "sanitization", "molecular_formula",  "molecular_bonds","tetrahedral_chirality", "double_bond_stereochemistry", "bond_lengths","bond_angles","aromatic_ring_flatness", "double_bond_flatness", "internal_steric_clash","internal_energy", "minimum_distance_to_protein","volume_overlap_with_protein"]

    # Filter
    plot_df = plot_df[plot_df["Item"].isin(selected_items)]

    # Rename each item individually
    item_rename_map = {
        "File_Loads": "File loads",
        "sanitization": "Sanitization",
        "molecular_formula": "Molecular formula",
        "molecular_bonds": "Bonds",
        "tetrahedral_chirality": "Tetrahedral chirality",
        "double_bond_stereochemistry": "Double bond stereochemistry",
        "bond_lengths": "Bond lengths",
        "bond_angles": "Bond angles",
        "aromatic_ring_flatness": "Planar aromatic rings",
        "double_bond_flatness": "Planer double bonds",
        "internal_steric_clash": "Internal steric clash",
        "internal_energy": "Energy ratio",
        "minimum_distance_to_protein": "Minimum protein-ligand distance",
        "volume_overlap_with_protein": "Volume overlap with protein"
    }

    plot_df["Item"] = plot_df["Item"].replace(item_rename_map)

    # Plot grouped barplot with hue=Group
    ax = sns.barplot(x="Item", y="Percentage", hue="Group", data=plot_df, palette=colors)
    sns.despine(top=True, right=True)
    ax.set_axisbelow(True)
    ax.tick_params(axis='both', which='major', length=6, width=1.2, color='black')

    plt.xlabel("")
    plt.ylabel("Passing Percentage (%)", fontsize=18)
    plt.ylim(0, 100)
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()

    # Move legend outside the plot
    plt.legend(title="Method", bbox_to_anchor=(1.05, 1), loc='upper left')

    # Ensure directory exists
    os.makedirs("./PoseBuster_Plots/", exist_ok=True)
    plt.savefig(f'./PoseBuster_Plots/PoseBuster_Bar_Combined.png', bbox_inches="tight")
    plt.close()



def Plot_Buster_Categories(dict_category_dict):
    # Extract categories and metrics
    categories = list(dict_category_dict.keys())
    metrics = list(next(iter(dict_category_dict.values())).keys())

    rename_map = {
        "AF": "AlphaFold 3",
        "CHAI": "Chai-1",
        "autodock": "AutoDock Vina",
        "rfaa": "RoseTTAFoldAA",
        "DD": "DiffDock-L"
    }
    x_labels = [rename_map.get(cat, cat) for cat in categories]

    # Pastel style: colors + markers for each metric
    styles = {
        "Chemical_Validity": ("o", "#DAA7B6"),
        "Intramolecular_Validity": ("s", "#DBCDF0"),
        "Intermolecular_Validity": ("^", "#F7D9C4")
    }

    # Convert data to a DataFrame for Seaborn
    data = []
    for cat in categories:
        for metric in metrics:
            data.append({
                "Category": rename_map.get(cat, cat),
                "Metric": metric.replace("_", " "),
                "Value": dict_category_dict[cat][metric],
                "Marker": styles[metric][0],
                "Color": styles[metric][1]
            })
    df = pd.DataFrame(data)

    # Seaborn style setup
    sns.set_style("whitegrid", {
        'grid.color': '0.85',
        'grid.linestyle': '-',
        'axes.edgecolor': '0.3',
        'axes.linewidth': 1,
        'xtick.bottom': True,
        'ytick.left': True})

    with plt.rc_context({
        'font.size': 16,
        'axes.titlesize': 20,
        'axes.labelsize': 18,
        'xtick.labelsize': 16,
        'ytick.labelsize': 16,
        'legend.fontsize': 16,
        'legend.title_fontsize': 16
    }):
        width = max(12, len(categories) * 2)
        plt.figure(figsize=(width, 6), dpi=300)

        #  Plot with seaborn
        ax = sns.scatterplot(
            data=df,
            x="Category",
            y="Value",
            hue="Metric",
            style="Metric",
            markers={m.replace("_", " "): styles[m][0] for m in metrics},
            palette={m.replace("_", " "): styles[m][1] for m in metrics},
            s=120,
            edgecolor='none'
        )

        # Clean up the plot
        sns.despine(top=True, right=True)
        ax.set_axisbelow(True)
        ax.tick_params(axis='both', which='major', length=6, width=1.2, color='black')
        plt.ylim(0, 110)
        plt.xlabel("", fontsize=18)
        plt.ylabel("Passing Percentage (%)", fontsize=18)
        plt.xticks(rotation=30, ha="right")
        plt.grid(True, linestyle="--", alpha=0.6)

        # Legend outside
        plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0., fontsize=16, title=None)
        plt.tight_layout()

        os.makedirs("./PoseBuster_Plots/", exist_ok=True)
        plt.savefig('./PoseBuster_Plots/PoseBuster_Category.png', bbox_inches="tight", dpi=500)
        plt.close()

if __name__ == "__main__":

    list_of_stats_dicts, dict_category_dict, prefix_list = [], {}, ["AF", "CHAI", "autodock", "rfaa", "DD"]
    colors = ['#a1c9f4', '#ffb482', '#8de5a1', '#ff9f9b', '#c9a0dc']

    for prefix in prefix_list:
        df = pd.read_csv(f'./{prefix}_results_short.csv')
        (stats_dict, stats_categories_dict) = Calculate_stats(df)
        list_of_stats_dicts.append(stats_dict)
        dict_category_dict[prefix] = stats_categories_dict

    PlotBuster(list_of_stats_dicts, prefix_list, colors)
    Plot_Buster_Categories(dict_category_dict)
    print('Finished plotting')
        