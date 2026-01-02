# ✈️ Kirchhoff-Love Plate Solver

A Python/Streamlit tool for analyzing aerospace structures (rectangular plates) based on the "Aircraft Structures" curriculum from the University of Liège (2013-2014).

[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](REPLACE_WITH_YOUR_STREAMLIT_APP_URL_HERE)

## 📖 Overview
This tool reproduces exact analytic solutions for:
1. **Pure Bending**: Navier solution for simply supported plates under constant pressure.
2. **Tension-Bending Coupling**: The effect of axial tension on bending stiffness (Geometric Stiffness).
3. **Initial Curvature**: Deformation of initially curved plates under tension.

## 🛠️ Installation (Local)
If you prefer to run this locally:

```bash
pip install -r requirements.txt
streamlit run plate_analysis_app.py