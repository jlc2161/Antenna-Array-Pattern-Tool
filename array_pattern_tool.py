import numpy as np
import streamlit as st
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from math import comb

st.set_page_config(page_title="Antenna Array Analysis", layout="wide")
st.title("Antenna Array Analysis (Element × Array Factor)")
st.markdown(
    "Select an **element pattern** and configure a **linear array** (N, spacing, phase). "
    "This multiplies the element pattern by the array factor.\n\n"
    "Note: This is an array-factor model (no mutual coupling/parasitics)."
)

# ---------- helpers ----------
def undb(x_db): 
    return 10**(x_db/10.0)

def dipole_halfwave(theta):
    """Half-wave thin dipole (axis = z)."""
    eps = 1e-6
    sth = np.sin(theta)
    sth = np.where(np.abs(sth) < eps, eps, sth)
    E = np.cos(0.5*np.pi*np.cos(theta)) / sth
    P = np.abs(E)**2
    return P/np.max(P)

def isotropic(theta):
    return np.ones_like(theta)

def array_factor(theta, N, d_lam, beta_rad, weights):
    """Linear array along z; AF = Σ w_n e^{j[k d n cosθ + β n]} (normalized)."""
    k = 2*np.pi
    n = np.arange(N) - (N-1)/2
    af = np.zeros_like(theta, dtype=complex)
    for i, ni in enumerate(n):
        af += weights[i] * np.exp(1j*(k*d_lam*ni*np.cos(theta) + beta_rad*ni))
    AF = np.abs(af)**2
    return AF/np.max(AF)

def hpbw(theta, P):
    Pn = P/np.max(P); idx = np.argmax(Pn); half = 0.5
    L = idx
    while L>0 and Pn[L]>=half: L -= 1
    R = idx
    while R<len(Pn)-1 and Pn[R]>=half: R += 1
    def edge(a,b):
        x1,y1 = theta[a], Pn[a]; x2,y2 = theta[b], Pn[b]
        return x1 + (half-y1)*(x2-x1)/((y2-y1)+1e-12)
    thL = theta[0] if L==0 else edge(L, L+1)
    thR = theta[-1] if R==len(Pn)-1 else edge(R-1, R)
    return np.degrees(thR-thL)

def directivity_est(theta, P):
    Pn = P/np.max(P)
    dth = theta[1]-theta[0]
    omega = 2*np.pi * np.sum(Pn*np.sin(theta))*dth
    return 4*np.pi/omega

def map_radius(P, mode, floor_db=-30.0):
    if mode == "Power (P)":   return P
    if mode == "Field (|E|)": return np.sqrt(P)
    Pdb = 10*np.log10(np.maximum(P, 1e-12))
    return np.clip((Pdb - floor_db)/(0 - floor_db), 0, 1)  # 0..1 from floor..0 dB

# ---------- sidebar ----------
with st.sidebar:
    st.header("Element Pattern")
    elem_choice = st.selectbox("Element", ["Half-wave Dipole", "Isotropic"])

    st.header("Array Preset")
    preset = st.selectbox(
        "Preset",
        [
            "Manual",
            "Single Dipole ",
            "Two-Element Phased ",
            "Two-Element In-Phase",
            "Four-Element Broadside",
            "8-Element Broadside",
            "End-Fire (+z, ordinary)",
            "End-Fire (−z, ordinary)",
            "End-Fire (+z, Hansen–Woodyard)",
            "Steer to 30° off broadside",
            "Binomial Broadside (no sidelobes)",
        ],
    )

    st.header("Display")
    floor_db   = st.slider("dB floor (for dB modes)", -60, -10, -30, 1)
    polar_mode = st.selectbox("2D polar scale", ["dB", "linear power"])
    radius_mode= st.selectbox("3D radius uses", ["Power (P)", "Field (|E|)", "dB (floor -30 dB)"])

    # Parameters (unrestricted N, like before)
    taper = "Uniform"
    if preset == "Single Dipole ":
        N, d_lam, beta_deg = 1, 0.5, 0.0
    elif preset == "Two-Element Phased ":
        N, d_lam, beta_deg = 2, 0.25, -90.0
    elif preset == "Two-Element In-Phase":
        N, d_lam, beta_deg = 2, 0.5, 0.0
    elif preset == "Four-Element Broadside":
        N, d_lam, beta_deg = 4, 0.5, 0.0
    elif preset == "8-Element Broadside":
        N, d_lam, beta_deg = 8, 0.5, 0.0
    elif preset == "End-Fire (+z, ordinary)":
        N, d_lam = 8, 0.25
        beta_deg = -360.0 * d_lam   # kd with θ0=0° ⇒ β = −kd
    elif preset == "End-Fire (−z, ordinary)":
        N, d_lam = 8, 0.25
        beta_deg = +360.0 * d_lam
    elif preset == "End-Fire (+z, Hansen–Woodyard)":
        N, d_lam = 8, 0.25
        beta_deg = -(360.0*d_lam + 180.0/N)
    elif preset == "Steer to 30° off broadside":
        N, d_lam = 8, 0.5
        theta0 = np.deg2rad(60.0)  # 30° off broadside (broadside=90°)
        beta_deg = -360.0 * d_lam * np.cos(theta0)
    elif preset == "Binomial Broadside (no sidelobes)":
        N, d_lam, beta_deg = 8, 0.5, 0.0
        taper = "Binomial"
    else:
        N       = st.number_input("Elements (N)", 1, 32, 4)
        d_lam   = st.slider("Spacing (λ)", 0.05, 1.50, 0.50, 0.01)
        beta_deg= st.slider("Progressive phase β (deg/element)", -180.0, 180.0, 0.0, 1.0)
        taper   = st.selectbox("Amplitude taper", ["Uniform", "Binomial"])

# ---------- compute (elevation cut) ----------
theta = np.linspace(1e-3, np.pi-1e-3, 2001)  # 0..π

Pelem = dipole_halfwave(theta) if elem_choice == "Half-wave Dipole" else isotropic(theta)

if taper == "Uniform":
    w = np.ones(N)
else:
    w = np.array([comb(N-1,k) for k in range(N)], float)
w = w/np.max(w)

AF = array_factor(theta, N, d_lam, np.radians(beta_deg), w)
P  = (Pelem * AF)
P  = P/np.max(P)  # always normalized

# ---------- plots ----------
c1, c2 = st.columns([1,1])

with c1:
    st.subheader("2D Elevation Cut ")
    theta_plot = np.concatenate([theta, 2*np.pi - theta[::-1]])
    P_plot     = np.concatenate([P,     P[::-1]])

    fig, ax = plt.subplots(figsize=(5,5), subplot_kw={'projection':'polar'})
    if polar_mode == "linear power":
        ax.plot(theta_plot, P_plot)
        ax.set_rticks([0.25, 0.5, 0.75, 1.0])
        ax.set_title("Linear power", pad=12)
    else:
        Pdb = 10*np.log10(np.maximum(P_plot, 1e-12))
        r = np.clip(Pdb - floor_db, 0, None)
        ax.plot(theta_plot, r)
        ticks = [0,10,20,30]
        ax.set_rticks(ticks)
        ax.set_yticklabels([f"{floor_db + t:.0f}" for t in ticks])
        ax.set_title("dB (normalized)", pad=12)
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_rlabel_position(135)
    st.pyplot(fig); plt.close(fig)

    st.markdown(f"**HPBW (approx):** {hpbw(theta,P):.1f}°  &nbsp;&nbsp; "
                f"**Directivity (est.):** {directivity_est(theta,P):.1f}")

with c2:
    st.subheader("3D Pattern")
    if radius_mode == "Power (P)":
        Rline = P
    elif radius_mode == "Field (|E|)":
        Rline = np.sqrt(P)
    else:
        Pdb = 10*np.log10(np.maximum(P, 1e-12))
        floor = -30.0
        Rline = np.clip((Pdb - floor)/(0 - floor), 0, 1)

    phi = np.linspace(0, 2*np.pi, 220)
    TH, PH = np.meshgrid(theta, phi)
    Rm = np.tile(Rline, (phi.size, 1))
    X = Rm * np.sin(TH) * np.cos(PH)
    Y = Rm * np.sin(TH) * np.sin(PH)
    Z = Rm * np.cos(TH)

    fig3d = go.Figure(data=[go.Surface(x=X, y=Y, z=Z, colorscale="Viridis", showscale=True)])
    fig3d.update_scenes(aspectmode="data")
    fig3d.update_layout(
        margin=dict(l=0,r=0,t=30,b=0),
        title="Normalized 3D Pattern (rotate/zoom)",
        scene_camera=dict(eye=dict(x=2.0, y=0.0, z=0.25), up=dict(x=0, y=0, z=1))
    )
    st.plotly_chart(fig3d, use_container_width=True)

st.caption(
    "Presets include NEETS examples (Single Dipole, Two-Element Phased) plus broadside, end-fire, "
    "steering, and binomial options. Patterns are normalized by default."
)
