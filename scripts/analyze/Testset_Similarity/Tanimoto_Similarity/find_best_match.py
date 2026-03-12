
import re
import pandas as pd

def Read_Cluster_Results():
    file_path = "../ProteinCartography/demo/cluster_benchmark/output/final_results/cluster_benchmark_aggregated_features.tsv"
    cluster_df = pd.read_csv(file_path, sep="\t")
    cluster_df = cluster_df[["protid", "LeidenCluster", "Dataset"]]
    return cluster_df


def Read_Lipid_Results():
    file_path = "CCD_similarity_atompair.csv"
    lipid_df = pd.read_csv(file_path)
    return lipid_df
    

def GetMatrix(cluster_df, lipid_df):
    # Get a list of LeidenCluster from the "LeidenCluster" column of cluster_df, remove duplicates
    list_clusters = sorted(cluster_df["LeidenCluster"].dropna().unique().tolist())
    # Loop over each LeidenCluster in the list
    for lei_cluster in list_clusters:
        df_processed = ProcessCluster(lei_cluster, cluster_df, lipid_df)
        df_processed.to_csv(f'{lei_cluster}_lipsimscore.csv', index=False)


def _extract_ccd(protid):
    token = str(protid).split("-")[-1].strip().upper()
    return re.sub(r"\d+$", "", token)


def ProcessCluster(lei_cluster, cluster_df, lipid_df):
    cluster_subset = cluster_df[cluster_df["LeidenCluster"] == lei_cluster]
    cluster_subset_test = cluster_subset[cluster_subset["Dataset"] == "Test"].copy()
    cluster_subset_general = cluster_subset[cluster_subset["Dataset"] == "General"].copy()

    if cluster_subset_test.empty:
        cluster_subset_test["Closest_tanimoto"] = pd.NA
        cluster_subset_test["ClosestData"] = pd.NA
        return cluster_subset_test

    sim_df = lipid_df.copy()
    if "CCD_test" in sim_df.columns:
        sim_df = sim_df.set_index("CCD_test")

    sim_df.index = sim_df.index.astype(str).str.strip().str.upper()
    sim_df.columns = sim_df.columns.astype(str).str.strip().str.upper()

    # Keep general rows (protid + CCD) in original order so ties resolve deterministically.
    general_pairs = (
        cluster_subset_general[["protid"]]
        .assign(CCD=lambda d: d["protid"].map(_extract_ccd).astype(str).str.strip().str.upper())
    )
    general_pairs = general_pairs[general_pairs["CCD"].isin(sim_df.columns)].copy()

    if general_pairs.empty:
        cluster_subset_test["Closest_tanimoto"] = pd.NA
        cluster_subset_test["ClosestData"] = pd.NA
        return cluster_subset_test

    def _get_score_and_data(protid):
        test_ccd = _extract_ccd(protid)
        if test_ccd not in sim_df.index:
            return pd.NA, pd.NA

        general_scores = pd.to_numeric(sim_df.loc[test_ccd, general_pairs["CCD"]], errors="coerce")
        if hasattr(general_scores, "to_numpy"):
            general_scores = general_scores.to_numpy()

        best_score = None
        best_protid = pd.NA
        for i, score in enumerate(general_scores):
            if pd.isna(score):
                continue
            if best_score is None or score > best_score:
                best_score = float(score)
                best_protid = general_pairs.iloc[i]["protid"]

        if best_score is None:
            return pd.NA, pd.NA
        return best_score, best_protid

    result = cluster_subset_test["protid"].map(_get_score_and_data)
    cluster_subset_test["Closest_tanimoto"] = result.map(lambda x: x[0])
    cluster_subset_test["ClosestData"] = result.map(lambda x: x[1])
    return cluster_subset_test




if __name__ == "__main__":
    # Create a new csv file with all test data and their accordingly most simliar lipid similarity score in the cluster
    cluster_df = Read_Cluster_Results()
    lipid_df = Read_Lipid_Results()
    GetMatrix(cluster_df, lipid_df)
