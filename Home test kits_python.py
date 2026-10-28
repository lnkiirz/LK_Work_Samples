import pandas as pd

# ============================================================
# INPUT FILES
# ============================================================

DIRECTORY = "//a70adom.bcbssc.com/R1/fs/bcbs2/data2/MEDICARE PART D/MEDBLUE_RX/STARS/Data Analyst Stuff/Tableau/Dashboards/Home test kits/"

input_everly = pd.read_excel(
    DIRECTORY + "input_everly.xlsx",
    sheet_name="Cumulative Summary",
    usecols="B:H",
    skiprows=22,  # Start at row 23
    nrows=8       # Through row 30
)

input_cologuard_2026 = pd.read_excel(
    DIRECTORY + "input_cologuard_2026.xlsx",
    sheet_name="Program Summary",
    usecols="B:F",
    skiprows=7,
    nrows=16
)

input_cologuard_2025 = pd.read_excel(
    DIRECTORY + "input_cologuard_2025.xlsx",
    sheet_name="Program Summary",
    usecols="B:F",
    skiprows=7,
    nrows=16
)

# ============================================================
# COLOGUARD 2026
# ============================================================

cologuard_2026_clean = (
    input_cologuard_2026[["Unnamed: 1", "TOTAL"]]
    .rename(columns={"Unnamed: 1": "Metric"})
)

cologuard_2026_clean = cologuard_2026_clean[
    cologuard_2026_clean["Metric"].notna()
].copy()

cologuard_2026_clean["Metric"] = (
    cologuard_2026_clean["Metric"].replace({
        "SHIPMENTS": "2026 Kits Shipped",
        "KITS RECEIVED": "2026 Returned Kits",
    })
)

# ============================================================
# COLOGUARD 2025
# ============================================================

cologuard_2025_clean = (
    input_cologuard_2025[["Unnamed: 1", "TOTAL"]]
    .rename(columns={"Unnamed: 1": "Metric"})
)

cologuard_2025_clean = cologuard_2025_clean[
    cologuard_2025_clean["Metric"].notna()
].copy()

cologuard_2025_clean["Metric"] = (
    cologuard_2025_clean["Metric"].replace({
        "SHIPMENTS": "2025 Kits Shipped",
        "KITS RECEIVED": "2025 Returned Kits",
    })
)

# ============================================================
# COMBINED COLOGUARD
# ============================================================

cologuard_clean = pd.concat(
    [cologuard_2025_clean, cologuard_2026_clean],
    ignore_index=True
)

cologuard_clean = cologuard_clean[
    cologuard_clean["Metric"].isin([
        "2026 Kits Shipped",
        "2026 Returned Kits",
        "2025 Kits Shipped",
        "2025 Returned Kits",
    ])
]

cologuard_clean = (
    cologuard_clean
    .set_index("Metric")["TOTAL"]
    .to_frame()
    .T
    .reset_index(drop=True)
)

numeric_cols = [
    "2026 Kits Shipped",
    "2026 Returned Kits",
    "2025 Kits Shipped",
    "2025 Returned Kits",
]

for column in numeric_cols:
    cologuard_clean[column] = pd.to_numeric(
        cologuard_clean[column],
        errors="coerce"
    )

cologuard_clean["Kit Type"] = "Cologuard"
cologuard_clean["Category"] = "COL"

cologuard_clean["2026 Return Rate"] = (
    cologuard_clean["2026 Returned Kits"] / cologuard_clean["2026 Kits Shipped"] * 100
)

cologuard_clean["2025 Return Rate"] = (
    cologuard_clean["2025 Returned Kits"] / cologuard_clean["2025 Kits Shipped"] * 100
)

cologuard_clean = cologuard_clean[
    [
        "Kit Type",
        "Category",
        "2026 Kits Shipped",
        "2026 Returned Kits",
        "2026 Return Rate",
        "2025 Kits Shipped",
        "2025 Returned Kits",
        "2025 Return Rate",
    ]
]

# ============================================================
# COL
# ============================================================

col = (
    input_everly.loc[input_everly["Kit Type"] == "Bag FIT"]
    .copy()
)

col["Category"] = "COL"

col["2026 Return Rate"] = col["2026 Return Rate"] * 100
col["2025 Return Rate"] = col["2025 Return Rate"] * 100

cols = col.columns.tolist()
cols.insert(1, cols.pop(cols.index("Category")))
col = col[cols]

col = pd.concat(
    [col, cologuard_clean],
    ignore_index=True
)

# ============================================================
# TOTAL COL
# ============================================================

kit_cols = [c for c in col.columns if "Kits" in c]

total_col = (
    col[["Kit Type", "Category", *kit_cols]]
    .select_dtypes(include="number")
    .sum()
    .to_frame()
    .T
)

total_col["Kit Type"] = "COL"
total_col["Category"] = "Total"

total_col["2025 Return Rate"] = (
    total_col["2025 Returned Kits"] / total_col["2025 Kits Shipped"] * 100
)

total_col["2026 Return Rate"] = (
    total_col["2026 Returned Kits"] / total_col["2026 Kits Shipped"] * 100
)

total_col = total_col[
    [
        "Kit Type",
        "Category",
        *[
            c
            for c in total_col.columns
            if c not in ["Kit Type", "Category"]
        ]
    ]
]

# ============================================================
# KED
# ============================================================

ked = input_everly.loc[
    (
        input_everly["Kit Type"].str.contains(
            r"UACR|eGFR|eGFR \+ UACR",
            regex=True,
            na=False
        )
    )
    &
    (
        ~input_everly["Kit Type"].str.contains(
            "A1c",
            na=False
        )
    )
].copy()

numeric_cols = ked.select_dtypes(include="number").columns

ked_final = (
    ked.groupby("Kit Type", as_index=False)[numeric_cols]
    .sum(min_count=1)
)

ked_final["Category"] = "KED"

ked_final["2025 Return Rate"] = ked_final["2025 Return Rate"] * 100
ked_final["2026 Return Rate"] = ked_final["2026 Return Rate"] * 100

cols = ked_final.columns.tolist()
cols.insert(1, cols.pop(cols.index("Category")))
ked_final = ked_final[cols]

# ============================================================
# TOTAL KED
# ============================================================

kit_cols = [c for c in ked_final.columns if "Kits" in c]

total_ked = (
    ked_final[["Kit Type", "Category", *kit_cols]]
    .select_dtypes(include="number")
    .sum()
    .to_frame()
    .T
)

total_ked["Kit Type"] = "KED"
total_ked["Category"] = "Total"

total_ked["2025 Return Rate"] = (
    total_ked["2025 Returned Kits"] / total_ked["2025 Kits Shipped"] * 100
)

total_ked["2026 Return Rate"] = (
    total_ked["2026 Returned Kits"] / total_ked["2026 Kits Shipped"] * 100
)

total_ked = total_ked[
    [
        "Kit Type",
        "Category",
        *[
            c
            for c in total_ked.columns
            if c not in ["Kit Type", "Category"]
        ]
    ]
]

# ============================================================
# A1C
# ============================================================

a1c = input_everly.loc[
    input_everly["Kit Type"].str.contains(
        "A1c",
        na=False
    )
].copy()

numeric_cols = a1c.select_dtypes(include="number").columns

a1c_final = (
    a1c.groupby("Kit Type", as_index=False)[numeric_cols]
    .sum(min_count=1)
)

a1c_final["Category"] = "A1C"

a1c_final["2025 Return Rate"] = a1c_final["2025 Return Rate"] * 100
a1c_final["2026 Return Rate"] = a1c_final["2026 Return Rate"] * 100

cols = a1c_final.columns.tolist()
cols.insert(1, cols.pop(cols.index("Category")))
a1c_final = a1c_final[cols]

# ============================================================
# TOTAL A1C
# ============================================================

kit_cols = [c for c in a1c_final.columns if "Kits" in c]

total_a1c = (
    a1c_final[["Kit Type", "Category", *kit_cols]]
    .select_dtypes(include="number")
    .sum()
    .to_frame()
    .T
)

total_a1c["Kit Type"] = "A1C"
total_a1c["Category"] = "Total"

total_a1c["2025 Return Rate"] = (
    total_a1c["2025 Returned Kits"] / total_a1c["2025 Kits Shipped"] * 100
)

total_a1c["2026 Return Rate"] = (
    total_a1c["2026 Returned Kits"] / total_a1c["2026 Kits Shipped"] * 100
)

total_a1c = total_a1c[
    [
        "Kit Type",
        "Category",
        *[
            c
            for c in total_a1c.columns
            if c not in ["Kit Type", "Category"]
        ]
    ]
]

# ============================================================
# TOTAL ALL
# ============================================================

combined = pd.concat(
    [
        col,
        ked_final,
        a1c_final,
    ],
    ignore_index=True
)

kit_cols = [c for c in combined.columns if "Kits" in c]

total_all = (
    combined[["Kit Type", "Category", *kit_cols]]
    .select_dtypes(include="number")
    .sum()
    .to_frame()
    .T
)

total_all["Kit Type"] = "All"
total_all["Category"] = "Total"

total_all["2025 Return Rate"] = (
    total_all["2025 Returned Kits"] / total_all["2025 Kits Shipped"] * 100
)

total_all["2026 Return Rate"] = (
    total_all["2026 Returned Kits"] / total_all["2026 Kits Shipped"] * 100
)

total_all = total_all[
    [
        "Kit Type",
        "Category",
        *[
            c
            for c in total_all.columns
            if c not in ["Kit Type", "Category"]
        ]
    ]
]

# ============================================================
# FINAL OUTPUT
# ============================================================

output = pd.concat(
    [
        col,
        ked_final,
        a1c_final,
        total_ked,
        total_col,
        total_a1c,
        total_all,
    ],
    ignore_index=True
)

value_cols = output.loc[
    :,
    "2026 Kits Shipped":"2025 Return Rate"
].columns

output = output.melt(
    id_vars=[
        c
        for c in output.columns
        if c not in value_cols
    ],
    value_vars=value_cols,
    var_name="name",
    value_name="value"
)

output["value"] = output["value"].round(2)

# ============================================================
# WRITE OUTPUT
# ============================================================

output.to_excel(
    DIRECTORY + "Test kits for Tableau.xlsx",
    index=False
)