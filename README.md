# Antenna Array Analysis (Element × Array Factor)

Interactive Streamlit app to explore antenna patterns from an **element pattern** multiplied by a **linear array factor** (ideal model, no coupling/parasitics).

---

## Quick Start

```bash
# 1) Clone & enter
git clone https://github.com/your-username/antenna-array-analysis.git
cd antenna-array-analysis

# 2) Install (Python 3.9+)
pip install -r requirements.txt
# if you don't have a file yet, minimal set:
# pip install streamlit numpy matplotlib plotly pandas

# 3) Run
streamlit run app.py
```

Open the browser link (usually `http://localhost:8501`).

---

## How to Use (UI Cheatsheet)

**Element Pattern**
- *Half-wave Dipole* (axis = z) or *Isotropic*.

**Array Preset**
- Pick a preset (Single Dipole, Two-Element, Broadside, End-Fire, Steer 30°, Binomial) **or**
- **Manual**: set `N` (elements), `d/λ` (spacing), `β` (deg/element), and `taper` (Uniform/Binomial).

**Display**
- 2D polar scale: **dB** or **linear power** (with adjustable dB floor).
- 3D radius: Power, |E|, or dB.
- Output shows **HPBW** and **Directivity** estimates.

> Geometry in this build: **array along z**, dipole axis = **z**, 3D pattern made by revolving the elevation cut.

---

## Typical Presets (What to Observe)

- **Single Dipole**: classic figure-8 in the horizontal plane; toroidal 3D “donut”.
- **Two-Element In-Phase (broadside)**: main lobe broadside; narrower beam than a single element.
- **Two-Element Phased**: introduces end-fire behavior with β ≈ −90° at d ≈ 0.25λ.
- **End-Fire (+z / −z)**: beam steered along the array axis; sign of β selects direction.
- **Binomial Broadside**: sidelobes suppressed; wider main beam.

---

## Math (short & simple)

- **Element (half-wave dipole)**  
  P(θ) ∝ |cos((π/2)cosθ) / sinθ|²  
  *Why:* dipole has nulls on axis, max broadside.

- **Array Factor (linear, along z)**  
  AF(θ) = |Σ wₙ e^{j(kdn cosθ + βn)}|² , k=2π, d in λ  
  *Why:* accounts for path difference + steering phase.*

- **Total Pattern**  
  P(θ) = P_elem(θ) × AF(θ) (normalized)  
  *Why:* element × array decomposition.*

- **HPBW**  
  Width at −3 dB around the main lobe.

- **Directivity (estimate)**  
  D ≈ 4π / ∫∫ Pn(θ) sinθ dθ dφ  
  *Why:* peak/average power over sphere.*

- **Steering Rule**  
  To steer toward θ₀: β = −kd cosθ₀.

---

## Screenshots (placeholders)

- **2D Elevation Cut (θ-polar)**  
  _Paste image here_

- **3D Pattern**  
  _Paste image here_

---

## Folder Hints

```
app.py                # Streamlit app (this code)
requirements.txt      # numpy, matplotlib, plotly, streamlit, pandas
README.md             # this file
screenshots/          # add your PNGs here
```

---

## Notes & Limits

- **Ideal** array-factor model: no ground, polarization plots, or mutual coupling.
- 3D surface = revolve the elevation pattern about z (fast + good for teaching).
- HPBW is from the **elevation cut**; if the lobe isn’t on that cut, HPBW may not reflect the true 3D width.
