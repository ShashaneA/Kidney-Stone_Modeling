# /// script
# requires-python = ">=3.13"
# dependencies = [
#     "marimo>=0.23.3",
#     "matplotlib>=3.10.9",
#     "numpy>=2.4.4",
#     "pandas>=3.0.3",
# ]
# ///

import marimo

__generated_with = "0.23.5"
app = marimo.App(width="medium")


@app.cell
def _():
    # KIDNEY STONE MODELING 
    # What this notebook does:
    # 1. Read 24-hour urine chemistry data
    # 2. Clean and convert units
    # 3. Send the chemistry data to EQUIL2 through R
    # 4. Use EQUIL2 supersaturation outputs to create a risk score
    # 5. Split samples into Low / Moderate / High risk groups
    # 6. Make clear figures for presentation
    import marimo as mo
    from pathlib import Path
    import shutil
    import subprocess

    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt

    return Path, np, pd, plt, shutil, subprocess


@app.cell
def _(Path):
    # SETTINGS

    # Put this notebook, the Excel data file, and run_equil2_batch.R in the same folder.

    # Try both names depending on how S1 Data is daved on your folder.

    if Path("S1_Data.xls").exists():
        DATA_FILE = "S1_Data.xls"
        DATA_SHEET = "Raw data"
    elif Path("S1_Data.xlsx").exists():
        DATA_FILE = "S1_Data.xlsx"
        DATA_SHEET = None
    else:
        raise FileNotFoundError(
            "could not find S1 Data.xls or S1_Data.xlsx in this folder."
        )

    R_SCRIPT_FILE = "run_equil2_batch.R"
    INPUT_FILE = "inputs.csv"
    OUTPUT_FILE = "outputs.csv"

    print("Using data file:", DATA_FILE)
    print("Using R script:", R_SCRIPT_FILE)
    return DATA_FILE, DATA_SHEET, INPUT_FILE, OUTPUT_FILE, R_SCRIPT_FILE


@app.cell
def _(DATA_FILE, DATA_SHEET, pd):
    # LOAD AND RENAME DATA

    # The original dataset uses longer column names. We rename only the columns needed for this project.

    if DATA_SHEET is None:
        raw_data = pd.read_excel(DATA_FILE)
    else:
        raw_data = pd.read_excel(DATA_FILE, sheet_name=DATA_SHEET)

    column_cleaned_data = raw_data.copy()
    column_cleaned_data.columns = (
        column_cleaned_data.columns
        .str.lower()
        .str.strip()
    )

    rename_map = {
        "volume (ml)": "volume_ml",
        "ca(mmol/24h)": "calcium_day",
        "ox(mmol/24h)": "oxalate_day",
        "cit(mmol/24h)": "citrate_day",
        "mg(mmol/24h)": "magnesium_day",
        "na(mmol/24h)": "sodium_day",
        "k(mmol/24h)": "potassium_day",
        "cl(mmol/24h)": "chloride_day",
        "p(mmol/24h)": "phosphate_day",
        "ua(mmol/24h)": "urate_day",
        "urine ph value": "ph",
    # These are original saturation columns from the data file. We rename them so they do not conflict with EQUIL2 outputs.
        "aps(caox)": "ss_caox_raw",
        "aps(cap)": "ss_cap_raw",
    }

    renamed_data = column_cleaned_data.rename(columns=rename_map)

    print("Available columns after cleaning:")
    print(renamed_data.columns.tolist())

    renamed_data.head()
    return (renamed_data,)


@app.cell
def _():
    # REQUIRED COLUMNS

    # These are the values needed to build the EQUIL2 input file.

    required_columns = [
        "volume_ml",
        "calcium_day",
        "oxalate_day",
        "citrate_day",
        "magnesium_day",
        "sodium_day",
        "potassium_day",
        "chloride_day",
        "phosphate_day",
        "urate_day",
        "ph",
    ]
    return (required_columns,)


@app.cell
def _(np, pd, renamed_data, required_columns):
    # CLEAN ROWS

    # This keeps only rows with the required chemistry values while also filters out impossible urine volume and unrealistic pH.

    def make_clean_urine_data(dataframe, columns_needed):
        cleaned = dataframe.copy()

        missing = [
            required_name
            for required_name in columns_needed
            if required_name not in cleaned.columns
        ]

        if len(missing) > 0:
            raise ValueError(f"Missing required columns: {missing}")

        cleaned = cleaned.replace([np.inf, -np.inf], np.nan)

        for required_name in columns_needed:
            cleaned[required_name] = pd.to_numeric(
                cleaned[required_name],
                errors="coerce"
            )

        cleaned = cleaned.dropna(subset=columns_needed)

        # Convert urine volume from mL to L.
        cleaned["volume_L"] = cleaned["volume_ml"] / 1000

        # Remove impossible volume and unrealistic pH values.
        cleaned = cleaned[cleaned["volume_L"] > 0]
        cleaned = cleaned[(cleaned["ph"] > 3) & (cleaned["ph"] < 10)]

        return cleaned

    clean_urine_data = make_clean_urine_data(renamed_data, required_columns)

    print(f"Valid rows after cleaning: {len(clean_urine_data)}")
    clean_urine_data.head()
    return (clean_urine_data,)


@app.cell
def _(clean_urine_data):
    # UNIT CONVERSION: mmol/day -> mmol/L

    # The dataset gives many ions as mmol per 24 hours. EQUIL2 needs concentration (either mmol/L or mmol/dL), so we divide by urine volume in L.

    def add_mmol_per_liter_columns(dataframe):
        converted_data = dataframe.copy()
        urine_volume_l = converted_data["volume_L"]

        converted_data["sodium_mmol_L"] = converted_data["sodium_day"] / urine_volume_l
        converted_data["potassium_mmol_L"] = converted_data["potassium_day"] / urine_volume_l
        converted_data["chloride_mmol_L"] = converted_data["chloride_day"] / urine_volume_l
        converted_data["calcium_mmol_L"] = converted_data["calcium_day"] / urine_volume_l
        converted_data["magnesium_mmol_L"] = converted_data["magnesium_day"] / urine_volume_l
        converted_data["phosphate_mmol_L"] = converted_data["phosphate_day"] / urine_volume_l
        converted_data["oxalate_mmol_L"] = converted_data["oxalate_day"] / urine_volume_l
        converted_data["citrate_mmol_L"] = converted_data["citrate_day"] / urine_volume_l
        converted_data["urate_mmol_L"] = converted_data["urate_day"] / urine_volume_l

    # These were not measured in the dataset, so the current model uses fixed values. This is an assumption and should be listed as a limitation.
        converted_data["ammonia_mmol_L"] = 30.0
        converted_data["sulfate_mmol_L"] = 10.0

        return converted_data

    model_input_data = add_mmol_per_liter_columns(clean_urine_data)

    model_input_data.head()
    return (model_input_data,)


@app.cell
def _(INPUT_FILE, model_input_data, pd):
    # CREATE EQUIL2 INPUT FILE

    # This CSV is the bridge between python and R where python prepares the data, and R runs EQUIL2.

    equil2_input_data = pd.DataFrame({
        "sodium_mmol_L": model_input_data["sodium_mmol_L"],
        "potassium_mmol_L": model_input_data["potassium_mmol_L"],
        "ammonia_mmol_L": model_input_data["ammonia_mmol_L"],
        "chloride_mmol_L": model_input_data["chloride_mmol_L"],
        "calcium_mmol_L": model_input_data["calcium_mmol_L"],
        "magnesium_mmol_L": model_input_data["magnesium_mmol_L"],
        "sulfate_mmol_L": model_input_data["sulfate_mmol_L"],
        "phosphate_mmol_L": model_input_data["phosphate_mmol_L"],
        "oxalate_mmol_L": model_input_data["oxalate_mmol_L"],
        "citrate_mmol_L": model_input_data["citrate_mmol_L"],
        "urate_mmol_L": model_input_data["urate_mmol_L"],
        "pH": model_input_data["ph"],
        "volume_L": model_input_data["volume_L"],
    })

    equil2_input_data.to_csv(INPUT_FILE, index=False)

    print(f"Saved {INPUT_FILE} with {len(equil2_input_data)} rows.")
    equil2_input_data.head()
    return


@app.cell
def _(OUTPUT_FILE, Path, R_SCRIPT_FILE, shutil, subprocess):

    # RUN EQUIL2 THROUGH R

    # This calls run_equil2_batch.R. The R file reads inputs.csv and writes outputs.csv.

    def find_rscript_executable():
        path_from_system = shutil.which("Rscript")
        if path_from_system is not None:
            return path_from_system

        r_folder = Path(r"C:\Program Files\R")
        matches = sorted(r_folder.glob(r"R-*\bin\Rscript.exe"))

        if len(matches) == 0:
            raise FileNotFoundError(
                "Could not find Rscript.exe. Make sure R is installed."
            )

        # Use the newest-looking R installation if more than one exists (since we made multiple for testing)
        return str(matches[-1])

    rscript_path = find_rscript_executable()
    script_path = Path(R_SCRIPT_FILE)

    if not script_path.exists():
        raise FileNotFoundError(
            f"Could not find {R_SCRIPT_FILE}. Put it in the same folder as this notebook."
        )

    print("Using Rscript:", rscript_path)
    print("Running R file:", script_path)

    r_result = subprocess.run(
        [rscript_path, "--vanilla", str(script_path)],
        check=True,
        capture_output=True,
        text=True,
    )

    print("R output:")
    print(r_result.stdout)

    if r_result.stderr.strip() != "":
        print("R messages/warnings:")
        print(r_result.stderr)

    if not Path(OUTPUT_FILE).exists():
        raise FileNotFoundError(f"R ran, but {OUTPUT_FILE} was not created.")
    return


@app.cell
def _(OUTPUT_FILE, model_input_data, pd):

    # LOAD EQUIL2 OUTPUTS AND COMBINE WITH ORIGINAL DATA

    equil2_output_data = pd.read_csv(OUTPUT_FILE)

    # Older R scripts may use this name so keeping this here to make the notebook robust.
    if "calcium_oxalate_ss" in equil2_output_data.columns:
        equil2_output_data = equil2_output_data.rename(
            columns={"calcium_oxalate_ss": "ss_caox"}
        )

    # Remove old or duplicate supersaturation columns from the original data.
    old_ss_columns = [
        current_name
        for current_name in model_input_data.columns
        if current_name.startswith("ss_") or current_name in ["ss_caox_raw", "ss_cap_raw"]
    ]

    data_without_old_ss = model_input_data.drop(
        columns=old_ss_columns,
        errors="ignore",
    )

    final_model_data = pd.concat(
        [
            data_without_old_ss.reset_index(drop=True),
            equil2_output_data.reset_index(drop=True),
        ],
        axis=1,
    )

    final_model_data = final_model_data.loc[:, ~final_model_data.columns.duplicated()]

    # Make sure these core output columns exist.
    if "ss_caox" not in final_model_data.columns:
        raise ValueError("EQUIL2 output is missing ss_caox.")

    if "ss_cap" not in final_model_data.columns:
        # If calcium phosphate was not returned directly, use brushite if available.
        if "ss_brushite" in final_model_data.columns:
            final_model_data["ss_cap"] = final_model_data["ss_brushite"]
        else:
            final_model_data["ss_cap"] = 0.0

    print("Final columns:")
    print(final_model_data.columns.tolist())

    print("\nMain supersaturation summary:")
    print(final_model_data[["ss_caox", "ss_cap"]].describe())

    final_model_data.head()
    return (final_model_data,)


@app.cell
def _(final_model_data, np, pd):

    # CREATE RISK SCORE AND RISK LABELS

    # This is a rule-based classification workflow, not a trained machine learning classifier yet.
    # The score uses EQUIL2 supersaturation plus urine chemistry values.

    def normalize_for_stone_model(values):
        numeric_values = pd.to_numeric(values, errors="coerce").fillna(0)
        value_range = numeric_values.max() - numeric_values.min()

        if value_range == 0:
            return np.zeros(len(numeric_values))

        return (numeric_values - numeric_values.min()) / value_range

    def build_risk_table(dataframe):
        risk_data = dataframe.copy()

    # Some correction terms: citrate can bind calcium, and magnesium can bind oxalate.
        risk_data["ca_free"] = (
            risk_data["calcium_mmol_L"] /
            (1 + 0.5 * risk_data["citrate_mmol_L"])
        )

        risk_data["ox_effective"] = (
            risk_data["oxalate_mmol_L"] /
            (1 + 0.1 * risk_data["magnesium_mmol_L"])
        )

    # Log-transform supersaturation so very large values do not dominate.
        risk_data["log_ss_caox"] = np.log1p(
            pd.to_numeric(risk_data["ss_caox"], errors="coerce")
            .fillna(0)
            .clip(lower=0)
        )

        risk_data["log_ss_cap"] = np.log1p(
            pd.to_numeric(risk_data["ss_cap"], errors="coerce")
            .fillna(0)
            .clip(lower=0)
        )

    # Risk terms: higher values increase the score.
        risk_data["calcium_risk"] = normalize_for_stone_model(risk_data["calcium_mmol_L"])
        risk_data["oxalate_risk"] = normalize_for_stone_model(risk_data["oxalate_mmol_L"])
        risk_data["caox_ss_risk"] = normalize_for_stone_model(risk_data["log_ss_caox"])
        risk_data["cap_ss_risk"] = normalize_for_stone_model(risk_data["log_ss_cap"])

    # Protection terms: higher values lower the score.
        risk_data["citrate_protection"] = normalize_for_stone_model(risk_data["citrate_mmol_L"])
        risk_data["magnesium_protection"] = normalize_for_stone_model(risk_data["magnesium_mmol_L"])

    # Weighted score
        risk_data["stone_risk_score"] = np.clip(
            (
                0.30 * risk_data["caox_ss_risk"] +
                0.20 * risk_data["calcium_risk"] +
                0.20 * risk_data["oxalate_risk"] +
                0.15 * risk_data["cap_ss_risk"] -
                0.10 * risk_data["citrate_protection"] -
                0.05 * risk_data["magnesium_protection"]
            ),
            0,
            1,
        )

        lower_cutoff = risk_data["stone_risk_score"].quantile(0.50)
        upper_cutoff = risk_data["stone_risk_score"].quantile(0.75)

        def label_one_score(score_value):
            if score_value >= upper_cutoff:
                return "High Risk"
            if score_value >= lower_cutoff:
                return "Moderate Risk"
            return "Low Risk"

        risk_data["stone_risk_label"] = risk_data["stone_risk_score"].apply(label_one_score)
        risk_data["patient_id"] = np.arange(1, len(risk_data) + 1)

        selected_columns = [
            "patient_id",
            "ph",
            "calcium_mmol_L",
            "oxalate_mmol_L",
            "citrate_mmol_L",
            "magnesium_mmol_L",
            "ss_caox",
            "ss_cap",
            "stone_risk_score",
            "stone_risk_label",
        ]

        sorted_results = risk_data[selected_columns].sort_values(
            "stone_risk_score",
            ascending=False,
        )

        return risk_data, sorted_results, lower_cutoff, upper_cutoff

    risk_model_data, risk_results, p50_score, p75_score = build_risk_table(final_model_data)

    risk_results.to_csv("patient_kidney_stone_risk_results.csv", index=False)

    print("Risk score distribution:")
    print(risk_model_data["stone_risk_score"].describe().round(4))
    print(f"\nThresholds: p50 = {p50_score:.4f}, p75 = {p75_score:.4f}")
    print("\nRisk counts:")
    print(risk_model_data["stone_risk_label"].value_counts())

    risk_results.head(20)
    return p50_score, p75_score, risk_model_data


@app.cell
def _(p50_score, p75_score, plt, risk_model_data):

    # FIGURE 1: RISK SCORE DISTRIBUTION

    fig, ax = plt.subplots(figsize=(10, 6))

    scores = risk_model_data["stone_risk_score"]

    ax.axvspan(0, p50_score, alpha=0.12, label=f"Low Risk (< {p50_score:.2f})")
    ax.axvspan(p50_score, p75_score, alpha=0.12, label=f"Moderate ({p50_score:.2f}–{p75_score:.2f})")
    ax.axvspan(p75_score, 1.0, alpha=0.12, label=f"High Risk (≥ {p75_score:.2f})")

    ax.axvline(p50_score, linestyle="--", linewidth=2)
    ax.axvline(p75_score, linestyle="--", linewidth=2)

    ax.hist(
        scores,
        bins=40,
        edgecolor="white",
        linewidth=0.6,
        alpha=0.85,
    )

    ax.set_title("Risk Scores Define the Risk Groups", fontsize=15)
    ax.set_xlabel("Stone Risk Score", fontsize=12)
    ax.set_ylabel("Number of Samples", fontsize=12)
    ax.legend()
    ax.grid(alpha=0.20)

    plt.tight_layout()
    plt.savefig("risk_score_distribution.png", dpi=150)
    plt.show()

    print("Saved: risk_score_distribution.png")
    return


@app.cell
def _(plt, risk_model_data):
    def _():

        # FIGURE 2: PATIENTS BY RISK CATEGORY

        category_order = ["High Risk", "Moderate Risk", "Low Risk"]
        risk_counts = risk_model_data["stone_risk_label"].value_counts()

        plot_labels = [label for label in category_order if label in risk_counts.index]
        plot_counts = [risk_counts[label] for label in plot_labels]

        fig, ax = plt.subplots(figsize=(9, 6))

        bars = ax.bar(plot_labels, plot_counts, width=0.55)

        for bar in bars:
            height = bar.get_height()
            percent = 100 * height / len(risk_model_data)
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                height + max(plot_counts) * 0.02,
                f"{int(height)}\n({percent:.1f}%)",
                ha="center",
                va="bottom",
                fontsize=11,
                fontweight="bold",
            )

        ax.set_title("Most Samples Were Classified as Low Risk", fontsize=15)
        ax.set_xlabel("Risk Category", fontsize=12)
        ax.set_ylabel("Number of Samples", fontsize=12)
        ax.set_ylim(0, max(plot_counts) * 1.20)
        ax.grid(axis="y", alpha=0.20)

        plt.tight_layout()
        plt.savefig("patients_by_risk_category.png", dpi=150)
        plt.show()
        return print("Saved: patients_by_risk_category.png")


    _()
    return


@app.cell
def _(pd, plt, risk_model_data):
    def _():
        # FIGURE 3: CaOx SUPERSATURATION BY pH AND RISK GROUP
        # The dashed SS = 1 line is a chemistry reference, not the class cutoff.

        plot_data = risk_model_data.copy()
        plot_data["ph"] = pd.to_numeric(plot_data["ph"], errors="coerce")
        plot_data["ss_caox"] = pd.to_numeric(plot_data["ss_caox"], errors="coerce")
        plot_data = plot_data.dropna(subset=["ph", "ss_caox", "stone_risk_label"])

        fig, ax = plt.subplots(figsize=(10, 6))

        risk_order = ["Low Risk", "Moderate Risk", "High Risk"]

        for risk_name in risk_order:
            group_data = plot_data[plot_data["stone_risk_label"] == risk_name]
            ax.scatter(
                group_data["ph"],
                group_data["ss_caox"],
                label=risk_name,
                alpha=0.65,
                s=45,
            )

        ax.axhline(
            y=1,
            linestyle="--",
            linewidth=2,
            label="SS = 1 reference line",
        )

        ax.set_title("Higher-Risk Samples Cluster at Higher CaOx Supersaturation", fontsize=15)
        ax.set_xlabel("Urine pH", fontsize=12)
        ax.set_ylabel("EQUIL2 Calcium Oxalate Supersaturation", fontsize=12)
        ax.legend()
        ax.grid(alpha=0.25)

        plt.tight_layout()
        plt.savefig("caox_supersaturation_by_risk.png", dpi=150)
        plt.show()
        return print("Saved: caox_supersaturation_by_risk.png")


    _()
    return


if __name__ == "__main__":
    app.run()
