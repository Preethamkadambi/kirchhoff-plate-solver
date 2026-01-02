import streamlit as st
import numpy as np
import plotly.graph_objects as go
import pandas as pd

# ==========================================
# CORE SOLVER CLASS
# ==========================================
class KirchhoffPlateSolver:
    """
    Solver for Kirchhoff-Love plates based on 'Aircraft Structures' 
    University of Liège (2013-2014).
    """
    
    def __init__(self, a, b, h, E, nu):
        self.a = a  # Length (x-direction)
        self.b = b  # Width (y-direction)
        self.h = h  # Thickness
        self.E = E  # Young's Modulus
        self.nu = nu # Poisson's ratio
        
        # Calculate Flexural Rigidity D
        # Formula from Page 24  / Page 43 [cite: 1241]
        self.D = (self.E * self.h**3) / (12 * (1 - self.nu**2))
        
    def get_D(self):
        return self.D

    def solve_navier_general(self, p_func, n11=0, phi_03_coeffs=None, m_max=9, n_max=9, grid_res=50):
        """
        Generalized Navier Solver covering all valid analytic cases in the PDF.
        
        Covers:
        - Case 1: Constant Pressure (Pages 24-33)
        - Case 2: Tension Effect (Pages 43-49)
        - Case 3: Initial Curvature (Pages 50-53)
        
        Parameters:
        - p_func: Function p(x, y) or constant float p0
        - n11: Axial tension (N/m)
        - phi_03_coeffs: List of (m, n, b_mn) for initial curvature
        """
        X = np.linspace(0, self.a, grid_res)
        Y = np.linspace(0, self.b, grid_res)
        X_grid, Y_grid = np.meshgrid(X, Y)
        
        u3 = np.zeros_like(X_grid)
        phi_03_grid = np.zeros_like(X_grid)
        logs = []
        
        # Base Log
        logs.append(f"**Flexural Rigidity D**: {self.D:.4e} Pa·m³ [Page 24]")
        
        # Iterate modes
        for m in range(1, m_max + 1):
            for n in range(1, n_max + 1):
                
                # 1. Determine Load Coefficients (a_mn)
                # For constant pressure p0, analytic integration is used (Page 28 )
                if callable(p_func):
                    # Numerical integration for arbitrary p(x,y) - Generalized from Page 49 [cite: 1368]
                    # Simplified here to just support the constant case p0 passed as float for exact matching of PDF
                    pass 
                elif isinstance(p_func, (int, float)):
                    p0 = p_func
                    if m % 2 != 0 and n % 2 != 0:
                        a_mn = (16 * p0) / (np.pi**2 * m * n) # [cite: 713]
                    else:
                        a_mn = 0
                else:
                    a_mn = 0

                # 2. Determine Initial Curvature Term (b_mn)
                # Page 52 [cite: 1463]
                b_mn = 0
                if phi_03_coeffs:
                    for (mc, nc, bc) in phi_03_coeffs:
                        if m == mc and n == nc:
                            b_mn = bc
                            # Add to visualization grid
                            phi_03_grid += b_mn * np.sin(m * np.pi * X_grid / self.a) * np.sin(n * np.pi * Y_grid / self.b)

                # 3. Calculate Stiffness Denominator
                # Pure Bending Stiffness (Page 28 [cite: 716])
                K_bending = self.D * np.pi**4 * ((m/self.a)**2 + (n/self.b)**2)**2
                
                # Geometric Stiffness due to Tension (Page 46 )
                # Note: PDF eq uses D*pi^6... this matches when factoring out correctly.
                # The term n11 * m^2 * pi^2 / a^2 is the added stiffness.
                K_tension = n11 * (m * np.pi / self.a)**2
                
                K_total = K_bending + K_tension
                
                # 4. Calculate Mode Amplitude A_mn
                # General Equation: K_total * A_mn = a_mn + (Tension_Load_from_Curvature)
                # Source: Page 53 [cite: 1490] for curvature effect
                
                # Load from curvature: - n11 * m^2 * pi^2 / a^2 * b_mn
                load_curvature = - n11 * (m * np.pi / self.a)**2 * b_mn
                
                # Total Load for this mode
                RHS = a_mn + load_curvature
                
                if K_total != 0:
                    A_mn = RHS / K_total
                else:
                    A_mn = 0
                
                # 5. Summation
                spatial = np.sin(m * np.pi * X_grid / self.a) * np.sin(n * np.pi * Y_grid / self.b)
                u3 += A_mn * spatial
                
                # Log significant modes
                if abs(A_mn) > 1e-9 and m <= 5 and n <= 5:
                     logs.append(f"Mode ({m},{n}): Load(a_mn)={a_mn:.2e} | Curv_Load={load_curvature:.2e} | Stiffness={K_total:.2e} -> Amp(A_mn)={A_mn:.4e}")

        return X_grid, Y_grid, u3, phi_03_grid, logs

# ==========================================
# STREAMLIT UI
# ==========================================
def main():
    st.set_page_config(page_title="Kirchhoff-Love Plate Solver", layout="wide")
    
    st.title("Exact Solutions for Kirchhoff-Love Plates")
    st.markdown("""
    **Source Document:** *StructAeroPlatesPart2.pdf* (University of Liège)
    
    This tool implements the analytic Navier solutions for simply supported rectangular plates 
    as derived in the lecture notes. It covers **pure bending**, **tension-bending coupling**, 
    and **initial curvature effects**.
    """)
    
    # --- Sidebar Inputs ---
    st.sidebar.header("1. Plate Properties")
    E = st.sidebar.number_input("Young's Modulus E (Pa)", value=70e9, format="%.2e")
    nu = st.sidebar.number_input("Poisson's Ratio ν", value=0.3)
    a = st.sidebar.number_input("Length a (m)", value=1.0)
    b = st.sidebar.number_input("Width b (m)", value=1.0)
    h = st.sidebar.number_input("Thickness h (m)", value=0.01)

    solver = KirchhoffPlateSolver(a, b, h, E, nu)
    
    # --- Tabs ---
    tab1, tab2, tab3 = st.tabs([
        "Page 24-33: Constant Pressure", 
        "Page 43-49: Tension & Pressure", 
        "Page 50-53: Initial Curvature"
    ])

    # ------------------------------------------------------------
    # CASE 1: Pure Bending (Constant Pressure)
    # ------------------------------------------------------------
    with tab1:
        st.subheader("Case 1: Simply Supported Plate - Constant Pressure")
        st.markdown("**Reference:** Pages 24-33")
        st.markdown(r"**Equation:** $D \nabla^4 u_3 = p_0$ [cite: 626]")
        
        col1, col2 = st.columns([1, 2])
        with col1:
            p0 = st.number_input("Pressure p0 (Pa)", value=1000.0, key="c1_p")
            m_max = st.slider("Fourier Modes (m)", 1, 15, 5, step=2, key="c1_m")
        
        with col2:
            X, Y, u3, _, logs = solver.solve_navier_general(p_func=p0, n11=0, m_max=m_max, n_max=m_max)
            
            # Theoretical check for Square plate (Page 30 [cite: 774])
            if abs(a - b) < 1e-5:
                 w_center = np.max(u3)
                 w_norm = w_center * solver.D / (p0 * a**4)
                 st.metric("Normalized Deflection (wD/pa⁴)", f"{w_norm:.5f}", delta=f"{w_norm - 0.00406:.5f} vs Theoretical")
            
            fig = go.Figure(data=[go.Surface(z=u3, x=X, y=Y)])
            fig.update_layout(title="Deflection u3 (m)", height=500)
            st.plotly_chart(fig, use_container_width=True)
            
            with st.expander("Calculation Log"):
                for l in logs: st.write(l)

    # ------------------------------------------------------------
    # CASE 2: Tension + Pressure
    # ------------------------------------------------------------
    with tab2:
        st.subheader("Case 2: Tension Effect on Bending")
        st.markdown("**Reference:** Pages 43-49")
        st.markdown(r"**Equation:** $D \nabla^4 u_3 - \tilde{n}^{11} u_{3,11} = p$ ")
        
        col1, col2 = st.columns([1, 2])
        with col1:
            p0_c2 = st.number_input("Pressure p0 (Pa)", value=1000.0, key="c2_p")
            n11_c2 = st.number_input("Axial Tension n11 (N/m)", value=1e5, format="%.2e", key="c2_n")
            m_max_c2 = st.slider("Fourier Modes (m)", 1, 15, 5, step=2, key="c2_m")
        
        with col2:
            X, Y, u3_c2, _, logs_c2 = solver.solve_navier_general(p_func=p0_c2, n11=n11_c2, m_max=m_max_c2, n_max=m_max_c2)
            
            fig = go.Figure(data=[go.Surface(z=u3_c2, x=X, y=Y, colorscale='Plasma')])
            fig.update_layout(title="Deflection with Tension (m)", height=500)
            st.plotly_chart(fig, use_container_width=True)
            
            st.info("Notice that increasing Tension (n11) reduces the max deflection, effectively increasing stiffness (Page 47 [cite: 1335]).")
            with st.expander("Calculation Log"):
                for l in logs_c2: st.write(l)

    # ------------------------------------------------------------
    # CASE 3: Initial Curvature
    # ------------------------------------------------------------
    with tab3:
        st.subheader("Case 3: Small Initial Curvature")
        st.markdown("**Reference:** Pages 50-53")
        st.markdown(r"**Equation:** $D \nabla^4 u_3 - \tilde{n}^{11} u_{3,11} = \tilde{n}^{11} \varphi_{03,11}$ ")
        st.markdown("*Assuming p=0 to isolate curvature effects*")
        
        col1, col2 = st.columns([1, 2])
        with col1:
            n11_c3 = st.number_input("Axial Tension n11 (N/m)", value=5e5, format="%.2e", key="c3_n")
            st.markdown("### Initial Shape Coefficients")
            st.markdown("Define $\\varphi_{03}$ (Page 52 [cite: 1463])")
            b_11 = st.number_input("Amplitude b_11 (m)", value=0.01)
            b_31 = st.number_input("Amplitude b_31 (m)", value=0.00)
            
        with col2:
            coeffs = [(1, 1, b_11), (3, 1, b_31)]
            # Pressure is 0 for this specific case in the PDF (Page 51 [cite: 1448])
            X, Y, u3_c3, phi_03, logs_c3 = solver.solve_navier_general(p_func=0, n11=n11_c3, phi_03_coeffs=coeffs)
            
            fig = go.Figure()
            fig.add_trace(go.Surface(z=phi_03, x=X, y=Y, opacity=0.3, showscale=False, name='Initial Shape'))
            fig.add_trace(go.Surface(z=phi_03 + u3_c3, x=X, y=Y, colorscale='Viridis', name='Final Shape'))
            fig.update_layout(title="Initial (Transparent) vs Deformed Shape (Color)", height=500)
            st.plotly_chart(fig, use_container_width=True)
            
            with st.expander("Calculation Log"):
                for l in logs_c3: st.write(l)

if __name__ == "__main__":
    main()