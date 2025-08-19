# Antenna Array Analysis (Element × Array Factor)

Interactive Streamlit app to explore antenna patterns from an **element pattern** multiplied by a **linear array factor**  
(ideal model, no coupling/parasitics).

---

## Quick Start

### Clone repository
```bash
git clone https://github.com/jlc2161/Antenna-Array-Pattern-Tool

cd antenna-array-pattern-tool

# Create virtual environment
python -m venv .venv

source .venv/Scripts/activate   # Windows
# or
source .venv/bin/activate       # macOS/Linux

# Install dependencies
pip install -r requirements.txt

# Run Application
streamlit run array_pattern_tool.py

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

## Screenshots

- **2D Elevation Cut (θ-polar)**
![2D Dipole Pattern](screenshots/2D%20Dipole.png)


- **3D Pattern**
![3D Dipole Pattern](screenshots/3D%20Dipole.png)

---

## Notes & Limits

- **Ideal** array-factor model: no ground, polarization plots, or mutual coupling.
- 3D surface = revolve the elevation pattern about z (fast + good for teaching).
- HPBW is from the **elevation cut**; if the lobe isn’t on that cut, HPBW may not reflect the true 3D width.
