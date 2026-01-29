import marimo

__generated_with = "0.19.6"
app = marimo.App(width="columns")


@app.cell(column=0)
def _():
    import marimo as mo
    import polars as pl
    from pathlib import Path
    return Path, mo


@app.cell
def _(mo):
    mo.md("""
    # Heart Transplant Recipients - CLIF Data Analysis

    This notebook loads and explores CLIF (Common Longitudinal ICU data Format) tables
    for characterizing heart transplant recipient hospitalizations.
    """)
    return


@app.cell
def _(Path):
    # Define data directory
    data_dir = Path("data")
    filetype = "parquet"
    timezone = "US/Central"
    return data_dir, filetype, timezone


@app.cell
def _(mo):
    mo.md("""
    ## Load CLIF Tables using clifpy
    """)
    return


@app.cell
def _():
    # Import all CLIF table classes from clifpy
    from clifpy.tables import (
        Patient,
        Hospitalization,
        Adt,
        Vitals,
        Labs,
        RespiratorySupport,
        Position,
        MedicationAdminContinuous,
        MedicationAdminIntermittent,
        PatientAssessments,
        HospitalDiagnosis,
        CodeStatus,
        CrrtTherapy,
        EcmoMcs,
        MicrobiologyCulture,
        MicrobiologyNonculture,
        PatientProcedures,
    )
    return (
        Adt,
        CodeStatus,
        CrrtTherapy,
        EcmoMcs,
        HospitalDiagnosis,
        Hospitalization,
        Labs,
        MedicationAdminContinuous,
        MedicationAdminIntermittent,
        MicrobiologyCulture,
        MicrobiologyNonculture,
        Patient,
        PatientAssessments,
        PatientProcedures,
        Position,
        RespiratorySupport,
        Vitals,
    )


@app.cell
def _(Patient, data_dir, filetype, timezone):
    # Load Patient table
    patient = Patient.from_file(
        data_directory=str(data_dir),
        filetype=filetype,
        timezone=timezone
    )
    patient.df
    return (patient,)


@app.cell
def _(Hospitalization, data_dir, filetype, timezone):
    # Load Hospitalization table
    hospitalization = Hospitalization.from_file(
        data_directory=str(data_dir),
        filetype=filetype,
        timezone=timezone
    )
    hospitalization.df
    return (hospitalization,)


@app.cell
def _(Adt, data_dir, filetype, timezone):
    # Load ADT (Admission-Discharge-Transfer) table
    adt = Adt.from_file(
        data_directory=str(data_dir),
        filetype=filetype,
        timezone=timezone
    )
    adt.df
    return (adt,)


@app.cell
def _(Vitals, data_dir, filetype, timezone):
    # Load Vitals table
    vitals = Vitals.from_file(
        data_directory=str(data_dir),
        filetype=filetype,
        timezone=timezone
    )
    vitals.df
    return (vitals,)


@app.cell
def _(Labs, data_dir, filetype, timezone):
    # Load Labs table
    labs = Labs.from_file(
        data_directory=str(data_dir),
        filetype=filetype,
        timezone=timezone
    )
    labs.df
    return (labs,)


@app.cell
def _(RespiratorySupport, data_dir, filetype, timezone):
    # Load Respiratory Support table
    respiratory_support = RespiratorySupport.from_file(
        data_directory=str(data_dir),
        filetype=filetype,
        timezone=timezone
    )
    respiratory_support.df
    return (respiratory_support,)


@app.cell
def _(Position, data_dir, filetype, timezone):
    # Load Position table
    position = Position.from_file(
        data_directory=str(data_dir),
        filetype=filetype,
        timezone=timezone
    )
    position.df
    return (position,)


@app.cell
def _(MedicationAdminContinuous, data_dir, filetype, timezone):
    # Load Continuous Medication Administration table
    medication_continuous = MedicationAdminContinuous.from_file(
        data_directory=str(data_dir),
        filetype=filetype,
        timezone=timezone
    )
    medication_continuous.df
    return (medication_continuous,)


@app.cell
def _(MedicationAdminIntermittent, data_dir, filetype, timezone):
    # Load Intermittent Medication Administration table
    medication_intermittent = MedicationAdminIntermittent.from_file(
        data_directory=str(data_dir),
        filetype=filetype,
        timezone=timezone
    )
    medication_intermittent.df
    return (medication_intermittent,)


@app.cell
def _(PatientAssessments, data_dir, filetype, timezone):
    # Load Patient Assessments table
    patient_assessments = PatientAssessments.from_file(
        data_directory=str(data_dir),
        filetype=filetype,
        timezone=timezone
    )
    patient_assessments.df
    return (patient_assessments,)


@app.cell
def _(HospitalDiagnosis, data_dir, filetype, timezone):
    # Load Hospital Diagnosis table
    hospital_diagnosis = HospitalDiagnosis.from_file(
        data_directory=str(data_dir),
        filetype=filetype,
        timezone=timezone
    )
    hospital_diagnosis.df
    return (hospital_diagnosis,)


@app.cell
def _(CodeStatus, data_dir, filetype, timezone):
    # Load Code Status table
    code_status = CodeStatus.from_file(
        data_directory=str(data_dir),
        filetype=filetype,
        timezone=timezone
    )
    code_status.df
    return (code_status,)


@app.cell
def _(CrrtTherapy, data_dir, filetype, timezone):
    # Load CRRT Therapy table
    crrt_therapy = CrrtTherapy.from_file(
        data_directory=str(data_dir),
        filetype=filetype,
        timezone=timezone
    )
    crrt_therapy.df
    return (crrt_therapy,)


@app.cell
def _(EcmoMcs, data_dir, filetype, timezone):
    # Load ECMO/MCS table
    ecmo_mcs = EcmoMcs.from_file(
        data_directory=str(data_dir),
        filetype=filetype,
        timezone=timezone
    )
    ecmo_mcs.df
    return (ecmo_mcs,)


@app.cell
def _(MicrobiologyCulture, data_dir, filetype, timezone):
    # Load Microbiology Culture table
    microbiology_culture = MicrobiologyCulture.from_file(
        data_directory=str(data_dir),
        filetype=filetype,
        timezone=timezone
    )
    microbiology_culture.df
    return (microbiology_culture,)


@app.cell
def _(MicrobiologyNonculture, data_dir, filetype, timezone):
    # Load Microbiology Nonculture table
    microbiology_nonculture = MicrobiologyNonculture.from_file(
        data_directory=str(data_dir),
        filetype=filetype,
        timezone=timezone
    )
    microbiology_nonculture.df
    return (microbiology_nonculture,)


@app.cell
def _(PatientProcedures, data_dir, filetype, timezone):
    # Load Patient Procedures table
    patient_procedures = PatientProcedures.from_file(
        data_directory=str(data_dir),
        filetype=filetype,
        timezone=timezone
    )
    patient_procedures.df
    return (patient_procedures,)


@app.cell
def _(mo):
    mo.md("""
    ## Summary of Loaded Tables
    """)
    return


@app.cell
def _(
    adt,
    code_status,
    crrt_therapy,
    ecmo_mcs,
    hospital_diagnosis,
    hospitalization,
    labs,
    medication_continuous,
    medication_intermittent,
    microbiology_culture,
    microbiology_nonculture,
    mo,
    patient,
    patient_assessments,
    patient_procedures,
    position,
    respiratory_support,
    vitals,
):
    # Create summary of all loaded tables
    table_summary = {
        "Table": [
            "Patient",
            "Hospitalization",
            "ADT",
            "Vitals",
            "Labs",
            "Respiratory Support",
            "Position",
            "Medication (Continuous)",
            "Medication (Intermittent)",
            "Patient Assessments",
            "Hospital Diagnosis",
            "Code Status",
            "CRRT Therapy",
            "ECMO/MCS",
            "Microbiology Culture",
            "Microbiology Nonculture",
            "Patient Procedures",
        ],
        "Rows": [
            len(patient.df),
            len(hospitalization.df),
            len(adt.df),
            len(vitals.df),
            len(labs.df),
            len(respiratory_support.df),
            len(position.df),
            len(medication_continuous.df),
            len(medication_intermittent.df),
            len(patient_assessments.df),
            len(hospital_diagnosis.df),
            len(code_status.df),
            len(crrt_therapy.df),
            len(ecmo_mcs.df),
            len(microbiology_culture.df),
            len(microbiology_nonculture.df),
            len(patient_procedures.df),
        ],
        "Columns": [
            len(patient.df.columns),
            len(hospitalization.df.columns),
            len(adt.df.columns),
            len(vitals.df.columns),
            len(labs.df.columns),
            len(respiratory_support.df.columns),
            len(position.df.columns),
            len(medication_continuous.df.columns),
            len(medication_intermittent.df.columns),
            len(patient_assessments.df.columns),
            len(hospital_diagnosis.df.columns),
            len(code_status.df.columns),
            len(crrt_therapy.df.columns),
            len(ecmo_mcs.df.columns),
            len(microbiology_culture.df.columns),
            len(microbiology_nonculture.df.columns),
            len(patient_procedures.df.columns),
        ],
    }
    mo.ui.table(table_summary)
    return


@app.cell(column=1)
def _():
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Identify the time of heart transplant

    Heart transplant recipients receive methylprednisolone 125mg q8h on their first post-operative day.
    We can use this steroid pattern to identify the transplant timing.
    """)
    return


@app.cell
def _(medication_intermittent):
    # Check what medication categories are available (clifpy uses pandas)
    med_categories = sorted(medication_intermittent.df["med_category"].unique())
    med_categories
    return


@app.cell
def _(medication_intermittent):
    # Filter to steroids only (using pandas syntax)
    steroids = medication_intermittent.df[
        medication_intermittent.df["med_category"] == "methylprednisolone"
    ]
    steroids
    return (steroids,)


@app.cell
def _(steroids):
    # Look at the dose distribution for methylprednisolone
    steroids.groupby("med_dose").size().sort_values(ascending=False)
    return


@app.cell
def _(steroids):
    # Filter to the 1000mg (1g) OR dose - this marks the transplant time
    methylpred_1g = steroids[steroids["med_dose"] == 1000].copy()
    methylpred_1g
    return (methylpred_1g,)


@app.cell
def _(methylpred_1g):
    # Get the first 1000mg dose per hospitalization as transplant cross clamp time
    transplant_times = (
        methylpred_1g
        .sort_values("admin_dttm")
        .groupby("hospitalization_id")
        .first()
        .reset_index()[["hospitalization_id", "admin_dttm"]]
        .rename(columns={"admin_dttm": "transplant_cross_clamp"})
    )
    transplant_times
    return (transplant_times,)


@app.cell
def _(hospitalization, mo, patient, transplant_times):
    n_patients = len(patient.df)
    n_hospitalizations = len(hospitalization.df)
    n_transplants = len(transplant_times)

    mo.md(f"""
    ## Transplant Times Identified

    | Count | Value |
    |-------|-------|
    | Patients | **{n_patients}** |
    | Hospitalizations | **{n_hospitalizations}** |
    | Transplants identified (1g methylpred) | **{n_transplants}** |
    | Coverage | **{n_transplants / n_patients * 100:.1f}%** of patients |
    """)
    return


@app.cell
def _(steroids, transplant_times):
    # Filter steroids to only hospitalizations with 1g push
    hosp_with_transplant = transplant_times["hospitalization_id"].tolist()
    steroids_transplant = steroids[
        steroids["hospitalization_id"].isin(hosp_with_transplant)
    ].copy()

    # Merge to get transplant time and calculate hours relative to transplant
    steroids_transplant = steroids_transplant.merge(transplant_times, on="hospitalization_id")
    steroids_transplant["hours_from_transplant"] = (
        (steroids_transplant["admin_dttm"] - steroids_transplant["transplant_cross_clamp"])
        .dt.total_seconds() / 3600
    )
    steroids_transplant
    return (steroids_transplant,)


@app.cell
def _(mo, steroids_transplant):
    # Dropdown to select hospitalization
    hosp_ids = sorted(steroids_transplant["hospitalization_id"].unique())
    hosp_selector = mo.ui.dropdown(
        options=hosp_ids,
        value=hosp_ids[0],
        label="Select Hospitalization"
    )
    hosp_selector
    return (hosp_selector,)


@app.cell
def _(hosp_selector, steroids_transplant):
    import altair as alt

    # Filter to selected hospitalization
    plot_data = steroids_transplant[
        steroids_transplant["hospitalization_id"] == hosp_selector.value
    ].copy()

    # Create scatter plot of methylpred dosing
    chart = alt.Chart(plot_data).mark_circle(size=100).encode(
        x=alt.X("hours_from_transplant:Q", title="Hours from Transplant"),
        y=alt.Y("med_dose:Q", title="Methylprednisolone Dose (mg)"),
        tooltip=["admin_dttm:T", "med_dose:Q", "med_name:N"]
    ).properties(
        title=f"Methylprednisolone Dosing - {hosp_selector.value}",
        width=600,
        height=300
    )

    # Add vertical line at transplant time (hour 0)
    vline = alt.Chart().mark_rule(color="red", strokeDash=[5,5]).encode(
        x=alt.datum(0)
    )

    chart + vline
    return (alt,)


@app.cell
def _(mo):
    mo.md("""
    ## Post-Transplant ICU Admission

    Find the first ICU ADT event after transplant_cross_clamp to identify when patient leaves OR.
    """)
    return


@app.cell
def _(adt):
    # Check ADT location categories
    adt_locations = sorted(adt.df["location_category"].unique())
    adt_locations
    return


@app.cell
def _(adt, transplant_times):
    # Filter ADT to ICU locations and merge with transplant times
    icu_adt = adt.df[adt.df["location_category"] == "icu"].copy()
    icu_adt = icu_adt.merge(transplant_times, on="hospitalization_id")

    # Filter to ICU admissions AFTER transplant cross clamp
    icu_adt = icu_adt[icu_adt["in_dttm"] > icu_adt["transplant_cross_clamp"]]

    # Get first ICU admission after transplant for each hospitalization
    post_transplant_icu = (
        icu_adt
        .sort_values("in_dttm")
        .groupby("hospitalization_id")
        .first()
        .reset_index()[["hospitalization_id", "in_dttm", "transplant_cross_clamp"]]
        .rename(columns={"in_dttm": "post_transplant_ICU_in_dttm"})
    )
    post_transplant_icu
    return (post_transplant_icu,)


@app.cell
def _(post_transplant_icu):
    # Calculate time from cross clamp to ICU (should be a few hours for OR time)
    post_transplant_icu["hours_in_OR"] = (
        (post_transplant_icu["post_transplant_ICU_in_dttm"] - post_transplant_icu["transplant_cross_clamp"])
        .dt.total_seconds() / 3600
    )
    post_transplant_icu[["hospitalization_id", "transplant_cross_clamp", "post_transplant_ICU_in_dttm", "hours_in_OR"]]
    return


@app.cell
def _(mo, post_transplant_icu):
    # Verify: OR time should typically be 2-8 hours
    median_or_time = post_transplant_icu["hours_in_OR"].median()
    min_or_time = post_transplant_icu["hours_in_OR"].min()
    max_or_time = post_transplant_icu["hours_in_OR"].max()

    mo.md(f"""
    ### Verification: Time in OR (Cross Clamp to ICU Arrival)

    | Metric | Hours |
    |--------|-------|
    | Median | **{median_or_time:.1f}** |
    | Min | **{min_or_time:.1f}** |
    | Max | **{max_or_time:.1f}** |

    *Expected: 2-8 hours for typical heart transplant surgery*
    """)
    return


@app.cell
def _(mo):
    mo.md("""
    ## Inotrope Analysis: Dopamine & Dobutamine

    Plot median hourly dose for first 72 hours after ICU admission.
    """)
    return


@app.cell
def _(medication_continuous):
    # Check continuous medication categories
    cont_med_categories = sorted(medication_continuous.df["med_category"].unique())
    cont_med_categories
    return


@app.cell
def _(medication_continuous, post_transplant_icu):
    # Filter to dopamine and dobutamine
    inotropes = medication_continuous.df[
        medication_continuous.df["med_category"].isin(["dopamine", "dobutamine"])
    ].copy()

    # Merge with ICU admission times
    inotropes = inotropes.merge(
        post_transplant_icu[["hospitalization_id", "post_transplant_ICU_in_dttm"]],
        on="hospitalization_id"
    )

    # Calculate hours from ICU admission
    inotropes["hours_from_icu"] = (
        (inotropes["admin_dttm"] - inotropes["post_transplant_ICU_in_dttm"])
        .dt.total_seconds() / 3600
    )

    # Filter to first 72 hours
    inotropes_72h = inotropes[
        (inotropes["hours_from_icu"] >= 0) & (inotropes["hours_from_icu"] <= 72)
    ].copy()

    # Create hour bins
    inotropes_72h["hour_bin"] = inotropes_72h["hours_from_icu"].astype(int)
    inotropes_72h
    return (inotropes_72h,)


@app.cell
def _(inotropes_72h):
    # Calculate mean hourly dose by medication
    hourly_mean = (
        inotropes_72h
        .groupby(["hour_bin", "med_category"])["med_dose"]
        .mean()
        .reset_index()
        .rename(columns={"med_dose": "mean_dose"})
    )
    hourly_mean
    return (hourly_mean,)


@app.cell
def _(alt, hourly_mean):
    # Plot mean hourly dopamine and dobutamine
    inotrope_chart_mean = alt.Chart(hourly_mean).mark_line(point=True).encode(
        x=alt.X("hour_bin:Q", title="Hours from ICU Admission"),
        y=alt.Y("mean_dose:Q", title="Mean Dose (mcg/kg/min)"),
        color=alt.Color("med_category:N", title="Medication"),
        tooltip=["hour_bin:Q", "med_category:N", "mean_dose:Q"]
    ).properties(
        title="Mean Hourly Dopamine & Dobutamine - First 72h Post-ICU",
        width=700,
        height=400
    )

    inotrope_chart_mean
    return


if __name__ == "__main__":
    app.run()
