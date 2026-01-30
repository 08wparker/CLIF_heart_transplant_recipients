import marimo

__generated_with = "0.19.6"
app = marimo.App(width="columns")


@app.cell(column=0)
def _():
    import marimo as mo
    import polars as pl
    import pandas as pd
    import numpy as np
    from pathlib import Path
    return Path, mo, np, pd


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
    # Filter to methylprednisolone
    steroids = medication_intermittent.df[
        medication_intermittent.df["med_category"] == "methylprednisolone"
    ]
    steroids
    return (steroids,)


@app.cell
def _(medication_intermittent):
    # Filter to prednisone (patients transition from methylpred to prednisone)
    prednisone = medication_intermittent.df[
        medication_intermittent.df["med_category"] == "prednisone"
    ]
    prednisone
    return (prednisone,)


@app.cell
def _(steroids):
    # Look at the dose distribution for methylprednisolone
    steroids.groupby("med_dose").size().sort_values(ascending=False)
    return


@app.cell
def _(hospitalization, steroids):
    # Filter to methylpred > 300mg (high-dose steroids indicating transplant)
    methylpred_high = steroids[steroids["med_dose"] > 300].copy()

    # Join with hospitalization to get patient_id
    methylpred_high = methylpred_high.merge(
        hospitalization.df[["hospitalization_id", "patient_id"]],
        on="hospitalization_id"
    )
    methylpred_high
    return (methylpred_high,)


@app.cell
def _(methylpred_high):
    # For each patient, find the first hospitalization with high-dose methylpred
    # Step 1: Get earliest high-dose methylpred time per hospitalization
    first_per_hosp = (
        methylpred_high
        .sort_values("admin_dttm")
        .groupby("hospitalization_id")
        .first()
        .reset_index()[["hospitalization_id", "patient_id", "admin_dttm"]]
    )

    # Step 2: For each patient, select the hospitalization with the earliest first dose
    transplant_times = (
        first_per_hosp
        .sort_values("admin_dttm")
        .groupby("patient_id")
        .first()
        .reset_index()[["patient_id", "hospitalization_id", "admin_dttm"]]
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
    | Transplants identified (methylpred >300mg) | **{n_transplants}** |
    | Match | **{"✓" if n_transplants == n_patients else "✗"}** (should equal patients) |
    """)
    return (n_patients,)


@app.cell
def _(steroids, transplant_times):
    # Filter steroids to only transplant hospitalizations (one per patient)
    transplant_hosp_ids = transplant_times["hospitalization_id"].tolist()
    steroids_transplant = steroids[
        steroids["hospitalization_id"].isin(transplant_hosp_ids)
    ].copy()

    # Merge to get transplant time and calculate hours relative to transplant
    steroids_transplant = steroids_transplant.merge(
        transplant_times[["hospitalization_id", "transplant_cross_clamp"]],
        on="hospitalization_id"
    )
    steroids_transplant["hours_from_transplant"] = (
        (steroids_transplant["admin_dttm"] - steroids_transplant["transplant_cross_clamp"])
        .dt.total_seconds() / 3600
    )
    steroids_transplant
    return steroids_transplant, transplant_hosp_ids


@app.cell
def _(mo, steroids_transplant):
    # Dropdown to select hospitalization
    _hosp_list = sorted(steroids_transplant["hospitalization_id"].unique())
    hosp_selector = mo.ui.dropdown(
        options=_hosp_list,
        value=_hosp_list[0],
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
    ## ADT Timeline with Methylprednisolone Dosing

    Visualize patient location changes overlaid with methylprednisolone administration.
    """)
    return


@app.cell
def _(mo, transplant_times):
    _hosp_with_transplant = transplant_times["hospitalization_id"].tolist()
    adt_hosp_selector = mo.ui.dropdown(
        options=_hosp_with_transplant,
        value=_hosp_with_transplant[0],
        label="Select Hospitalization ID"
    )
    adt_hosp_selector
    return (adt_hosp_selector,)


@app.cell
def _(mo):
    adt_time_range = mo.ui.range_slider(
        start=-200,
        stop=800,
        value=[-100, 200],
        step=10,
        label="Time Range (hours from transplant)"
    )
    adt_time_range
    return (adt_time_range,)


@app.cell
def _(
    adt,
    adt_hosp_selector,
    adt_time_range,
    alt,
    prednisone,
    steroids,
    transplant_times,
):
    _selected_hosp = adt_hosp_selector.value
    _time_min, _time_max = adt_time_range.value

    _adt_data = adt.df[adt.df["hospitalization_id"] == _selected_hosp].copy()
    _steroid_data = steroids[steroids["hospitalization_id"] == _selected_hosp].copy()
    _prednisone_data = prednisone[prednisone["hospitalization_id"] == _selected_hosp].copy()
    _transplant_time = transplant_times[
        transplant_times["hospitalization_id"] == _selected_hosp
    ]["transplant_cross_clamp"].iloc[0]

    _adt_data["hours_from_transplant"] = (
        (_adt_data["in_dttm"] - _transplant_time).dt.total_seconds() / 3600
    )
    _adt_data["out_hours_from_transplant"] = (
        (_adt_data["out_dttm"] - _transplant_time).dt.total_seconds() / 3600
    )
    _steroid_data["hours_from_transplant"] = (
        (_steroid_data["admin_dttm"] - _transplant_time).dt.total_seconds() / 3600
    )
    _prednisone_data["hours_from_transplant"] = (
        (_prednisone_data["admin_dttm"] - _transplant_time).dt.total_seconds() / 3600
    )

    # Filter to time range (keep ADT if it overlaps with the window)
    _adt_data = _adt_data[
        (_adt_data["out_hours_from_transplant"] >= _time_min) &
        (_adt_data["hours_from_transplant"] <= _time_max)
    ]
    _steroid_data = _steroid_data[
        (_steroid_data["hours_from_transplant"] >= _time_min) &
        (_steroid_data["hours_from_transplant"] <= _time_max)
    ]
    _prednisone_data = _prednisone_data[
        (_prednisone_data["hours_from_transplant"] >= _time_min) &
        (_prednisone_data["hours_from_transplant"] <= _time_max)
    ]

    _location_order = ["or", "icu", "ward", "stepdown", "pacu", "ed", "procedural", "other"]
    _existing_locations = [loc for loc in _location_order if loc in _adt_data["location_category"].unique()]
    _other_locations = [loc for loc in _adt_data["location_category"].unique() if loc not in _location_order]
    _all_locations = _existing_locations + _other_locations
    _all_locations.extend(["methylpred", "prednisone"])

    _x_scale = alt.Scale(domain=[_time_min, _time_max])

    _adt_bars = alt.Chart(_adt_data).mark_bar(height=25).encode(
        x=alt.X("hours_from_transplant:Q", title="Hours from Transplant", scale=_x_scale),
        x2="out_hours_from_transplant:Q",
        y=alt.Y("location_category:N", title="Location", sort=_all_locations),
        color=alt.Color("location_category:N", legend=None),
        tooltip=[
            alt.Tooltip("location_category:N", title="Location"),
            alt.Tooltip("in_dttm:T", title="In Time"),
            alt.Tooltip("out_dttm:T", title="Out Time"),
            alt.Tooltip("hours_from_transplant:Q", title="In (hours)", format=".1f"),
            alt.Tooltip("out_hours_from_transplant:Q", title="Out (hours)", format=".1f")
        ]
    )

    _steroid_plot_data = _steroid_data.copy()
    _steroid_plot_data["location_category"] = "methylpred"

    _steroid_points = alt.Chart(_steroid_plot_data).mark_point(
        size=100,
        shape="diamond",
        filled=True
    ).encode(
        x=alt.X("hours_from_transplant:Q", scale=_x_scale),
        y=alt.Y("location_category:N", sort=_all_locations),
        color=alt.Color("med_dose:Q", scale=alt.Scale(scheme="reds"), title="Methylpred (mg)"),
        tooltip=[
            alt.Tooltip("admin_dttm:T", title="Admin Time"),
            alt.Tooltip("med_dose:Q", title="Dose (mg)"),
            alt.Tooltip("hours_from_transplant:Q", title="Hours", format=".1f")
        ]
    )

    _prednisone_plot_data = _prednisone_data.copy()
    _prednisone_plot_data["location_category"] = "prednisone"

    _prednisone_points = alt.Chart(_prednisone_plot_data).mark_point(
        size=100,
        shape="circle",
        filled=True
    ).encode(
        x=alt.X("hours_from_transplant:Q", scale=_x_scale),
        y=alt.Y("location_category:N", sort=_all_locations),
        color=alt.Color("med_dose:Q", scale=alt.Scale(scheme="blues"), title="Prednisone (mg)"),
        tooltip=[
            alt.Tooltip("admin_dttm:T", title="Admin Time"),
            alt.Tooltip("med_dose:Q", title="Dose (mg)"),
            alt.Tooltip("hours_from_transplant:Q", title="Hours", format=".1f")
        ]
    )

    _transplant_rule = alt.Chart().mark_rule(
        color="darkgreen",
        strokeWidth=2,
        strokeDash=[5, 5]
    ).encode(x=alt.datum(0))

    _combined_chart = alt.layer(
        _adt_bars,
        _steroid_points,
        _prednisone_points,
        _transplant_rule
    ).properties(
        title=f"ADT Timeline with Steroids - {_selected_hosp}",
        width=800,
        height=400
    ).resolve_scale(color='independent')

    _combined_chart
    return


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

    Plot median hourly dose for first week (168 hours) after ICU admission.
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

    # Filter to first week (168 hours)
    inotropes_168h = inotropes[
        (inotropes["hours_from_icu"] >= 0) & (inotropes["hours_from_icu"] <= 168)
    ].copy()

    # Create hour bins
    inotropes_168h["hour_bin"] = inotropes_168h["hours_from_icu"].astype(int)
    inotropes_168h
    return


@app.cell
def _(patient_hourly):
    # Calculate mean hourly dose by medication from imputed data
    hourly_mean_imputed = (
        patient_hourly
        .groupby(["hour_bin", "med_category"])["avg_dose_locf"]
        .mean()
        .reset_index()
        .rename(columns={"avg_dose_locf": "mean_dose"})
    )
    hourly_mean_imputed
    return (hourly_mean_imputed,)


@app.cell
def _(alt, hourly_mean_imputed):
    # Plot mean hourly dopamine and dobutamine from imputed data
    inotrope_chart_mean = alt.Chart(hourly_mean_imputed).mark_line(point=True).encode(
        x=alt.X("hour_bin:Q", title="Hours from ICU Admission"),
        y=alt.Y("mean_dose:Q", title="Mean Dose (mcg/kg/min)"),
        color=alt.Color("med_category:N", title="Medication"),
        tooltip=["hour_bin:Q", "med_category:N", alt.Tooltip("mean_dose:Q", format=".2f")]
    ).properties(
        title="Mean Hourly Vasoactive Medication Dose - First Week Post-ICU",
        width=800,
        height=400
    )

    inotrope_chart_mean
    return


@app.cell
def _(patient_hourly):
    # Calculate median hourly dose by medication from imputed data
    hourly_median_imputed = (
        patient_hourly
        .groupby(["hour_bin", "med_category"])["avg_dose_locf"]
        .median()
        .reset_index()
        .rename(columns={"avg_dose_locf": "median_dose"})
    )
    hourly_median_imputed
    return (hourly_median_imputed,)


@app.cell
def _(alt, hourly_median_imputed):
    # Plot median hourly vasoactive medication doses
    inotrope_chart_median = alt.Chart(hourly_median_imputed).mark_line(point=True).encode(
        x=alt.X("hour_bin:Q", title="Hours from ICU Admission"),
        y=alt.Y("median_dose:Q", title="Median Dose"),
        color=alt.Color("med_category:N", title="Medication"),
        tooltip=["hour_bin:Q", "med_category:N", alt.Tooltip("median_dose:Q", format=".3f")]
    ).properties(
        title="Median Hourly Vasoactive Medication Dose - First Week Post-ICU",
        width=800,
        height=400
    )

    inotrope_chart_median
    return


@app.cell
def _(patient_hourly):
    # Calculate percentage of patients on each medication per hour
    # A patient is "on" a medication if dose > 0
    patient_hourly_on = patient_hourly.copy()
    patient_hourly_on["on_med"] = patient_hourly_on["avg_dose_locf"] > 0

    # Count patients on each medication per hour
    hourly_counts = (
        patient_hourly_on
        .groupby(["hour_bin", "med_category"])
        .agg(
            n_on=("on_med", "sum"),
            n_total=("on_med", "count")
        )
        .reset_index()
    )
    hourly_counts["pct_on"] = 100 * hourly_counts["n_on"] / hourly_counts["n_total"]
    hourly_counts
    return (hourly_counts,)


@app.cell
def _(alt, hourly_counts):
    # Plot percentage of patients on each vasoactive medication
    vasoactive_pct_chart = alt.Chart(hourly_counts).mark_line(point=True).encode(
        x=alt.X("hour_bin:Q", title="Hours from ICU Admission"),
        y=alt.Y("pct_on:Q", title="% Patients on Medication"),
        color=alt.Color("med_category:N", title="Medication"),
        tooltip=[
            "hour_bin:Q",
            "med_category:N",
            alt.Tooltip("pct_on:Q", title="% On", format=".1f"),
            alt.Tooltip("n_on:Q", title="N On"),
            alt.Tooltip("n_total:Q", title="N Total")
        ]
    ).properties(
        title="Percentage of Patients on Vasoactive Medications - First Week Post-ICU",
        width=800,
        height=400
    )

    vasoactive_pct_chart
    return


@app.cell
def _(mo):
    mo.md("""
    ## Patient-Level Imputed Inotrope Dataset

    Create a dataset with 168 hourly observations (1 week) per patient using LOCF imputation.
    - If MAR action is "not given", dose = 0
    - Use LOCF to fill missing hours
    """)
    return


@app.cell
def _(medication_continuous):
    # Check columns in continuous medication table
    medication_continuous.df.columns.tolist()
    return


@app.cell
def _(medication_continuous):
    # Check MAR action categories
    if "mar_action_category" in medication_continuous.df.columns:
        mar_actions = sorted(medication_continuous.df["mar_action_category"].dropna().unique())
    else:
        mar_actions = "Column not found"
    mar_actions
    return


@app.cell
def _():
    # Define vasoactive medications to track
    vasoactive_meds = [
        "dopamine",
        "dobutamine",
        "epinephrine",
        "norepinephrine",
        "vasopressin",
        "milrinone",
        "nitric_oxide"
    ]
    vasoactive_meds
    return (vasoactive_meds,)


@app.cell
def _(medication_continuous, np, post_transplant_icu, vasoactive_meds):
    # Filter vasoactive medications and merge with ICU times
    vasoactive_raw = medication_continuous.df[
        medication_continuous.df["med_category"].isin(vasoactive_meds)
    ].copy()

    vasoactive_raw = vasoactive_raw.merge(
        post_transplant_icu[["hospitalization_id", "post_transplant_ICU_in_dttm"]],
        on="hospitalization_id"
    )

    # Calculate hours from ICU and create hour bin
    vasoactive_raw["hours_from_icu"] = (
        (vasoactive_raw["admin_dttm"] - vasoactive_raw["post_transplant_ICU_in_dttm"])
        .dt.total_seconds() / 3600
    )
    vasoactive_raw["hour_bin"] = vasoactive_raw["hours_from_icu"].apply(lambda x: int(x) if x >= 0 else -1)

    # Filter to first week (0-167 hours)
    vasoactive_raw = vasoactive_raw[(vasoactive_raw["hour_bin"] >= 0) & (vasoactive_raw["hour_bin"] < 168)]

    # Handle MAR action - if "not given" or similar, set dose to 0
    if "mar_action_category" in vasoactive_raw.columns:
        vasoactive_raw["effective_dose"] = np.where(
            vasoactive_raw["mar_action_category"].str.lower().str.contains("not given|held|stopped", na=False),
            0,
            vasoactive_raw["med_dose"]
        )
    else:
        vasoactive_raw["effective_dose"] = vasoactive_raw["med_dose"]

    vasoactive_raw[["hospitalization_id", "hour_bin", "med_category", "med_dose", "effective_dose"]].head(20)
    return (vasoactive_raw,)


@app.cell
def _(pd, transplant_hosp_ids, vasoactive_meds, vasoactive_raw):
    # Create skeleton: all patients x 168 hours (1 week) x all vasoactive medications
    hours = list(range(168))

    skeleton = pd.DataFrame([
        {"hospitalization_id": h, "hour_bin": hr, "med_category": med}
        for h in transplant_hosp_ids
        for hr in hours
        for med in vasoactive_meds
    ])

    # Calculate mean effective dose per patient-hour-medication from raw data
    hourly_doses = (
        vasoactive_raw
        .groupby(["hospitalization_id", "hour_bin", "med_category"])["effective_dose"]
        .mean()
        .reset_index()
        .rename(columns={"effective_dose": "avg_dose"})
    )

    # Merge with skeleton
    patient_hourly = skeleton.merge(hourly_doses, on=["hospitalization_id", "hour_bin", "med_category"], how="left")

    # Imputation: backfill first, then forward fill (LOCF)
    # This ensures early hours before first observation get backfilled from first dose
    patient_hourly = patient_hourly.sort_values(["hospitalization_id", "med_category", "hour_bin"])

    # Step 1: Backfill (fills early NAs from first observation)
    patient_hourly["avg_dose_imputed"] = (
        patient_hourly
        .groupby(["hospitalization_id", "med_category"])["avg_dose"]
        .bfill()
    )

    # Step 2: Forward fill (LOCF for any remaining gaps)
    patient_hourly["avg_dose_locf"] = (
        patient_hourly
        .groupby(["hospitalization_id", "med_category"])["avg_dose_imputed"]
        .ffill()
    )

    # Fill remaining NaN with 0 (patient never received this inotrope)
    patient_hourly["avg_dose_locf"] = patient_hourly["avg_dose_locf"].fillna(0)
    patient_hourly = patient_hourly.drop(columns=["avg_dose_imputed"])

    patient_hourly
    return (patient_hourly,)


@app.cell
def _(patient_hourly):
    # Pivot to wide format: one row per patient-hour, columns for each medication
    vasoactive_wide = patient_hourly.pivot_table(
        index=["hospitalization_id", "hour_bin"],
        columns="med_category",
        values="avg_dose_locf"
    ).reset_index()

    vasoactive_wide.columns.name = None
    # Add _dose suffix to all medication columns
    vasoactive_wide.columns = [
        f"{col}_dose" if col not in ["hospitalization_id", "hour_bin"] else col
        for col in vasoactive_wide.columns
    ]

    vasoactive_wide
    return (vasoactive_wide,)


@app.cell
def _(mo, n_patients, vasoactive_wide):
    n_patients_vaso = vasoactive_wide["hospitalization_id"].nunique()
    n_rows_vaso = len(vasoactive_wide)
    expected_rows_vaso = n_patients * 168

    mo.md(f"""
    ### Imputed Vasoactive Dataset Summary

    | Metric | Value |
    |--------|-------|
    | Patients | **{n_patients_vaso}** |
    | Rows | **{n_rows_vaso}** |
    | Expected (patients × 168h) | **{expected_rows_vaso}** |
    | Match | **{"✓" if n_rows_vaso == expected_rows_vaso else "✗"}** |
    """)
    return


if __name__ == "__main__":
    app.run()
