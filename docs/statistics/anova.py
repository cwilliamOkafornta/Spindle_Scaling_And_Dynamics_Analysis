import pandas as pd
import statsmodels.api as sm

df = pd.read_csv("./Supplementary_table_WT.csv")
df = df.rename(columns={"Final pole_pole length (µm)": "Length"})
df = df[["Stage", "Cell", "Length"]].dropna()
df["Stage"] = df["Stage"].astype(str)
df["Cell"]  = df["Cell"].astype(str)
df["All"]   = 1  # single dummy group; we only want variance components

# Random intercepts for Stage and for Cell (no interaction term)
m = sm.MixedLM.from_formula(
    "Length ~ 1",
    groups="All",           # dummy grouping
    re_formula="0",
    vc_formula={"Stage": "0 + C(Stage)",
                "Cell":  "0 + C(Cell)"},
    data=df
)
res = m.fit(reml=True, method="lbfgs", maxiter=1000)
print(res.summary())

# Extract variance components (Stage, Cell) and residual (measurement noise)
print("Residual variance (measurement noise):", float(res.scale))
# In recent statsmodels, the variance components live here:
print("Variance components:", getattr(res, "vcomp", "upgrade statsmodels to access res.vcomp reliably"))
