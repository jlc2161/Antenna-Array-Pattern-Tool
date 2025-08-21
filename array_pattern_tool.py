import numpy as np
import streamlit as st
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from math import comb

# ──────────────────────────────────────────────────────────────────────────────
# App setup
# ──────────────────────────────────────────────────────────────────────────────
st.set_page_config(page_title="Antenna Array Analysis", layout="wide")
st.title("Antenna Array Analysis (Element × Array Factor)")
st.markdown(
    "Pick an **element pattern** and configure a **linear array** (N, spacing, phase). "
    "The total pattern = element pattern × array factor. "
    "*This is an array‑factor model (no mutual coupling/parasitics).*"
)

# ──────────────────────────────────────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────────────────────────────────────
def dipole_halfwave(theta):
    """Half‑wave thin dipole element (axis = z). Power pattern, normalized."""
    eps = 1e-6
    sth = np.sin(theta)
    sth = np.where(np.abs(sth) < eps, eps, sth)
    E = np.cos(0.5*np.pi*np.cos(theta)) / sth
    P = np.abs(E)**2
    return P / np.max(P)

def isotropic(theta):
    """Isotropic element: constant power."""
    return np.ones_like(theta)

def array_factor(theta, N, d_lam, beta_rad, weights):
    """
    Linear array along z:
      AF(θ) = | Σ w_n · exp{ j[k d n cosθ + β n] } |^2  (normalized)
    """
    k = 2*np.pi  # since d is in λ
    n = np.arange(N) - (N-1)/2
    af = np.zeros_like(theta, dtype=complex)
    for i, ni in enumerate(n):
        af += weights[i] * np.exp(1j*(k*d_lam*ni*np.cos(theta) + beta_rad*ni))
    AF = np.abs(af)**2
    return AF / np.max(AF)

def hpbw(theta, P):
    """Half‑power beamwidth from a 2D elevation cut."""
    Pn = P/np.max(P)
    idx = np.argmax(Pn)
    half = 0.5
    L = idx
    while L > 0 and Pn[L] >= half:
        L -= 1
    R = idx
    while R < len(Pn)-1 and Pn[R] >= half:
        R += 1
    def edge(a,b):
        x1,y1 = theta[a], Pn[a]; x2,y2 = theta[b], Pn[b]
        return x1 + (half-y1)*(x2-x1)/((y2-y1)+1e-12)
    thL = theta[0] if L==0 else edge(L, L+1)
    thR = theta[-1] if R==len(Pn)-1 else edge(R-1, R)
    return np.degrees(thR-thL)

def directivity_est(theta, P):
    """Crude directivity estimate from a single‑cut pattern."""
    Pn = P/np.max(P)
    dth = theta[1]-theta[0]
    omega = 2*np.pi * np.sum(Pn*np.sin(theta)) * dth
    return 4*np.pi/omega

# ──────────────────────────────────────────────────────────────────────────────
# Sidebar
# ──────────────────────────────────────────────────────────────────────────────
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

    # NEW: 3D smoothing controls
    st.subheader("3D quality")
    phi_samples = st.slider("Azimuth samples (φ)", 120, 540, 300, 20)
    theta_upscale = st.slider("Elevation upsample ×", 1, 4, 2, 1)
    colorscale_choice = st.selectbox(
        "Colors (3D)",
        ["Viridis (classic)", "Directivity: Red→Green→Blue"]
    )

    # Parameters (presets)
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

    st.divider()
    st.subheader("Physical overlay")
    show_phys_2d = st.checkbox("Show physical elements (2D)", value=True)
    dip_len_lam  = st.slider("Element length (λ) for dipole", 0.3, 1.2, 0.5, 0.01)
    scale_lam_to_r = st.slider("Overlay scale: 1.0 r ≈ ? λ", 0.5, 2.0, 1.0, 0.01)
    stick_radius = st.slider("Element radius (λ, visual)", 0.005, 0.05, 0.02, 0.005)

# ──────────────────────────────────────────────────────────────────────────────
# Compute (elevation cut, with optional upsampling for smoother 3D)
# ──────────────────────────────────────────────────────────────────────────────
theta_base = np.linspace(1e-3, np.pi-1e-3, 2001)            # base resolution for 2D
theta_fine = np.linspace(theta_base[0], theta_base[-1],
                         theta_base.size * int(theta_upscale))  # upsample for 3D smoothness

# Element power pattern
if elem_choice == "Half-wave Dipole":
    Pelem_base = dipole_halfwave(theta_base)
    Pelem_fine = dipole_halfwave(theta_fine)
else:
    Pelem_base = isotropic(theta_base)
    Pelem_fine = isotropic(theta_fine)

# Weights
if taper == "Uniform":
    w = np.ones(N)
else:  # binomial
    w = np.array([comb(N-1,k) for k in range(N)], float)
w = w / np.max(w)

# Array factor
AF_base = array_factor(theta_base, N, d_lam, np.radians(beta_deg), w)
AF_fine = array_factor(theta_fine, N, d_lam, np.radians(beta_deg), w)

# Total (normalized)
P2D = (Pelem_base * AF_base); P2D /= np.max(P2D)
P3Dline = (Pelem_fine * AF_fine); P3Dline /= np.max(P3Dline)

# ──────────────────────────────────────────────────────────────────────────────
# Layout columns
# ──────────────────────────────────────────────────────────────────────────────
c1, c2 = st.columns([1,1])

# ──────────────────────────────────────────────────────────────────────────────
# 2D Polar (with physical element overlay for BOTH linear and dB)
# ──────────────────────────────────────────────────────────────────────────────
with c1:
    st.subheader("2D Elevation Cut")
    theta_plot = np.concatenate([theta_base, 2*np.pi - theta_base[::-1]])
    P_plot     = np.concatenate([P2D,        P2D[::-1]])

    fig, ax = plt.subplots(figsize=(6.2,6.2), subplot_kw={'projection':'polar'})

    if polar_mode == "linear power":
        ax.plot(theta_plot, P_plot, lw=2)
        ax.set_rticks([0.25, 0.5, 0.75, 1.0])
        ax.set_title("Linear power", pad=12)
    else:
        Pdb = 10*np.log10(np.maximum(P_plot, 1e-12))
        r = np.clip(Pdb - floor_db, 0, None)
        ax.plot(theta_plot, r, lw=2)
        ticks = [0,10,20,30]
        ax.set_rticks(ticks)
        ax.set_yticklabels([f"{floor_db + t:.0f}" for t in ticks])
        ax.set_title("dB (normalized)", pad=12)

    # styling
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_rlabel_position(135)

    # Physical overlay (dipole sticks along z) — shown for both linear & dB
    if show_phys_2d:
        # draw array axis at θ=0 and θ=π
        ax.plot([0,0], [0,1.05], color="k", lw=2, alpha=0.65)
        ax.plot([np.pi,np.pi], [0,1.05], color="k", lw=2, alpha=0.65)

        n_idx = np.arange(N) - (N-1)/2
        z_lam = n_idx * d_lam
        r_centers = np.abs(z_lam/scale_lam_to_r)
        r_max = 1.05
        half_len_r = (dip_len_lam/2) / scale_lam_to_r
        for rc in r_centers:
            if rc < r_max - 1e-3:
                for th in (0.0, np.pi):
                    r1 = max(rc - half_len_r, 0.0)
                    r2 = min(rc + half_len_r, r_max)
                    ax.plot([th, th], [r1, r2], color="k", lw=6, alpha=0.85, solid_capstyle="round")
                    ax.plot([th], [rc], marker="o", ms=8, color="k", alpha=0.9)

    st.pyplot(fig); plt.close(fig)

    st.markdown(
        f"**HPBW (approx):** {hpbw(theta_base,P2D):.1f}°  &nbsp;&nbsp; "
        f"**Directivity (est.):** {directivity_est(theta_base,P2D):.1f}"
    )

# ──────────────────────────────────────────────────────────────────────────────
# 3D Surface (smoother, with lighting & colorscale options)
# ──────────────────────────────────────────────────────────────────────────────
with c2:
    st.subheader("3D Pattern")

    # Choose how to map radius
    if radius_mode == "Power (P)":
        Rline = P3Dline
    elif radius_mode == "Field (|E|)":
        Rline = np.sqrt(P3Dline)
    else:
        Pdb = 10*np.log10(np.maximum(P3Dline, 1e-12))
        floor = -30.0
        Rline = np.clip((Pdb - floor)/(0 - floor), 0, 1)

    # Build a smoother surface: φ dense, θ from fine grid (upsampled)
    phi = np.linspace(0, 2*np.pi, int(phi_samples))
    TH, PH = np.meshgrid(theta_fine, phi)
    Rm = np.tile(Rline, (phi.size, 1))

    X = Rm * np.sin(TH) * np.cos(PH)
    Y = Rm * np.sin(TH) * np.sin(PH)
    Z = Rm * np.cos(TH)

    # Colorscale choice
    if colorscale_choice.startswith("Viridis"):
        colorscale = "Viridis"
    else:
        # Red (max) -> Green -> Blue (min)  for teaching “directivity”
        colorscale = [
            [0.00, "rgb(0, 0, 255)"],   # blue (low)
            [0.50, "rgb(0, 255, 0)"],  # green (mid)
            [1.00, "rgb(255, 0, 0)"],  # red (high)
        ]

    surf = go.Surface(
        x=X, y=Y, z=Z,
        colorscale=colorscale,
        showscale=True,
        # Nice lighting for smooth appearance
        lighting=dict(ambient=0.65, diffuse=0.75, specular=0.25, roughness=0.8, fresnel=0.1),
        lightposition=dict(x=1.0, y=0.0, z=2.0)
    )

    fig3d = go.Figure(data=[surf])

    # Add physical dipoles as sticks (same geometry as 2D overlay)
    n_idx = np.arange(N) - (N-1)/2
    z_lam = n_idx * d_lam
    for zi in z_lam:
        z1 = (zi - 0.5*dip_len_lam) / scale_lam_to_r
        z2 = (zi + 0.5*dip_len_lam) / scale_lam_to_r
        fig3d.add_trace(go.Scatter3d(
            x=[0,0], y=[0,0], z=[z1,z2],
            mode="lines",
            line=dict(width=8, color="black"),
            showlegend=False
        ))

    fig3d.update_scenes(
        aspectmode="data",
        xaxis_title="x(λ)", yaxis_title="y(λ)", zaxis_title="z(λ)",
        xaxis=dict(showgrid=False), yaxis=dict(showgrid=False), zaxis=dict(showgrid=False)
    )
    fig3d.update_layout(
        margin=dict(l=0,r=0,t=30,b=0),
        title="Normalized 3D Pattern (rotate/zoom)",
        scene_camera=dict(eye=dict(x=1.7, y=1.7, z=1.45), up=dict(x=0, y=0, z=1)),
    )

    st.plotly_chart(fig3d, use_container_width=True)

st.caption(
    "2D overlay shows a **linear dipole array** placed along the z‑axis with your chosen spacing (λ) and length. "
    "3D shows both the pattern surface and the element sticks. Use **3D quality** controls in the sidebar to make the surface smoother. "
    "In the **Directivity** palette, **red = stronger/higher**, **blue = weaker/lower**."
)
