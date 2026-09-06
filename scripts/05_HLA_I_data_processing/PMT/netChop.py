import os
import pandas as pd


# ============================================================
# 1. Merge NetChop web outputs
# ============================================================

os.chdir("./MS_unspecific/circRNA/result/PMT/netchop/result")

netChop_1 = pd.read_excel("./web_1.xlsx")
netChop_2 = pd.read_excel("./web_2.xlsx")

all_netChop = pd.concat(
    [netChop_1, netChop_2],
    ignore_index=True
)

all_netChop.to_csv(
    "all_netChop.txt",
    sep="\t",
    index=False
)


# ============================================================
# 2. Parse NetChop output into a clean table
# ============================================================

def parse_protein_data(file_path):
    data = []

    with open(file_path, "r") as f:
        for line in f:
            if line.strip() and line.strip()[0].isdigit():
                parts = line.strip().split()

                if len(parts) >= 5:
                    try:
                        pos = int(parts[0])
                        aa = parts[1]
                        c = parts[2]
                        score = float(parts[3])
                        ident = parts[4]
                        data.append([pos, aa, c, score, ident])
                    except Exception as e:
                        print(
                            f"Skipping line: {line.strip()} because of error: {e}"
                        )

    return pd.DataFrame(
        data,
        columns=["Position", "AA", "C", "Score", "ORF_ID"]
    )


df_1_clean = parse_protein_data("./all_netChop.txt")

df_1_clean.to_csv(
    "./MS_unspecific/circRNA/result/PMT/netchop/result/filter_result/deal_HLA_1.txt",
    sep="\t",
    index=False,
    header=True
)


# ============================================================
# 3. Restore full ORF IDs
# ============================================================

os.chdir(
    "./MS_unspecific/circRNA/result/PMT/netchop/result/filter_result"
)

deal_HLA_1 = pd.read_csv(
    "deal_HLA_1.txt",
    sep="\t",
    header=0
)

mapping_ID = pd.read_csv(
    "./MS_unspecific/circRNA/result/PMT/netchop/NetChop_ID_mapping.csv"
)

new_ID_deal_HLA_1 = deal_HLA_1.copy().merge(
    mapping_ID,
    left_on="ORF_ID",
    right_on="short_id",
    how="left"
)

# Rename the original ORF ID column for compatibility with NetMHCpan results
new_ID_deal_HLA_1 = new_ID_deal_HLA_1.rename(
    columns={"original_id": "new_id"}
)

new_ID_deal_HLA_1.to_csv(
    "./MS_unspecific/circRNA/result/PMT/netchop/result/filter_result/deal_HLA_1_with_full_ID.txt",
    sep="\t",
    index=False,
    header=True
)


# ============================================================
# 4. Prepare HLA-I peptide coordinates
# ============================================================

peptide_path = (
    "./MS_unspecific/circRNA/result/PMT/netchop/"
    "peptide_Cloud_HLA_1.txt"
)

peptide_seq = pd.read_csv(
    peptide_path,
    header=0,
    sep="\t"
)

peptide_df = peptide_seq.copy()

peptide_df["start"] = peptide_df["Pos"] + 1
peptide_df["Length"] = peptide_df["Peptide"].str.len()
peptide_df["end"] = peptide_df["Pos"] + peptide_df["Length"]

peptide_df.to_csv(
    "./MS_unspecific/circRNA/result/PMT/netchop/result/filter_result/HLA_1_peptide.txt",
    index=False,
    sep="\t",
    header=True
)


# ============================================================
# 5. Evaluate precise and flanking NetChop cleavage support
# ============================================================

def analyze_peptide_precise_and_fuzzy_cleavage(
    peptide_df,
    netchop_df,
    flank=5,
    threshold=0.5
):
    result = []

    nc = netchop_df.copy()
    nc["Position"] = pd.to_numeric(
        nc["Position"],
        errors="coerce"
    ).astype("Int64")
    nc["Score"] = pd.to_numeric(
        nc["Score"],
        errors="coerce"
    )

    nc = nc.dropna(
        subset=["Position", "Score"]
    ).copy()
    nc["Position"] = nc["Position"].astype(int)
    nc["Score"] = nc["Score"].astype(float)

    groups = {
        key: group
        for key, group in nc.groupby("new_id", sort=False)
    }

    for _, row in peptide_df.iterrows():
        orf_id = row["ID"]
        peptide = row["Peptide"]
        start = int(row["start"])
        end = int(row["end"])
        hla_type = row["type"]
        hla_rank = row["X_EL_Rank"]

        sub = groups.get(orf_id)

        if sub is None or sub.empty:
            result.append({
                "ORF_ID": orf_id,
                "Peptide": peptide,
                "HLA_type": hla_type,
                "HLA_Rank": hla_rank,
                "Start": start,
                "End": end,
                "Fuzzy_Left": False,
                "Fuzzy_Right": False,
                "Precise_Cleavage": False,
                "Best_Left_Score": None,
                "Best_Right_Score": None,
                "Match_Type": "No_NetChop"
            })
            continue

        orf_max_pos = int(sub["Position"].max())

        high_cuts = set(
            sub.loc[
                sub["Score"] > threshold,
                "Position"
            ].tolist()
        )

        left_window = list(
            range(
                max(1, start - flank),
                start
            )
        )

        right_window = list(
            range(
                end,
                min(end + flank + 1, orf_max_pos + 1)
            )
        )

        left_match = any(
            pos in high_cuts
            for pos in left_window
        )

        right_match = any(
            pos in high_cuts
            for pos in right_window
        )

        left_scores = sub.loc[
            sub["Position"].isin(left_window),
            "Score"
        ]

        right_scores = sub.loc[
            sub["Position"].isin(right_window),
            "Score"
        ]

        # Precise cleavage requires cleavage at both peptide boundaries:
        # N-terminal boundary at start - 1 and C-terminal boundary at end.
        precise = (
            (start - 1) in high_cuts
            and end in high_cuts
        )

        if precise:
            match_type = "Precise"
        elif left_match or right_match:
            match_type = "Fuzzy"
        else:
            match_type = "None"

        result.append({
            "ORF_ID": orf_id,
            "Peptide": peptide,
            "HLA_type": hla_type,
            "HLA_Rank": hla_rank,
            "Start": start,
            "End": end,
            "Fuzzy_Left": left_match,
            "Fuzzy_Right": right_match,
            "Precise_Cleavage": precise,
            "Best_Left_Score": (
                float(left_scores.max())
                if not left_scores.empty
                else None
            ),
            "Best_Right_Score": (
                float(right_scores.max())
                if not right_scores.empty
                else None
            ),
            "Match_Type": match_type
        })

    return pd.DataFrame(result)


HLA_1_result = analyze_peptide_precise_and_fuzzy_cleavage(
    peptide_df,
    new_ID_deal_HLA_1,
    flank=5,
    threshold=0.5
)

u_HLA_1_result = HLA_1_result.drop_duplicates()

u_HLA_1_result.to_csv(
    "./MS_unspecific/circRNA/result/PMT/netchop/result/filter_result/netChop_result_1.txt",
    sep="\t",
    index=False,
    header=True
)
