from manim import *

# ================================================================
# Linear Algebra Animations for Jupyter Book (Manim CE)
# Author: Prepared for Stanford Numerical Linear Algebra course
# ================================================================
# Notes:
# - Compatible with Manim Community Edition.
# - Each scene is self-contained; render with:
#     manim -pqh linear_transforms_manim.py GridDeformationScene
#   or replace the scene name accordingly.
# - Matrices and vectors are ASCII-only.
# - Feel free to modify matrices below to explore different effects.
# ================================================================

# -----------------------------
# Helpers
# -----------------------------
def matrix_tex(M, label=None):
    """Return a MathTex showing a 2x2 matrix M with an optional label."""
    a,b = M[0]
    c,d = M[1]
    if label is None:
        tex = MathTex(
            "\\begin{bmatrix}"
            f"{a:.1f} & {b:.1f} \\\\ {c:.1f} & {d:.1f}"
            "\\end{bmatrix}"
        )
    else:
        tex = MathTex(
            label, "=",
            "\\begin{bmatrix}"
            f"{a:.1f} & {b:.1f} \\\\ {c:.1f} & {d:.1f}"
            "\\end{bmatrix}"
        )
    return tex

def vector_tex(v, name=None):
    """Return a MathTex column vector with optional name."""
    x, y = v
    if name is None:
        return MathTex(
            "\\begin{bmatrix}"
            f"{x:.1f} \\\\ {y:.1f}"
            "\\end{bmatrix}"
        )
    return MathTex(
        name, "=",
        "\\begin{bmatrix}"
        f"{x:.1f} \\\\ {y:.1f}"
        "\\end{bmatrix}"
    )

def arrow_vec(start, vec, **kwargs):
    """Draw an arrow from start to start+vec."""
    return Arrow(start, start + vec, buff=0, **kwargs)

# -----------------------------
# 1) Grid Deformation: T(x) = A x
# -----------------------------
class GridDeformationScene(Scene):
    """
    Visualize T(x) = A x by deforming a grid (NumberPlane) and basis vectors.
    Modify A to explore different linear transformations.
    """
    def construct(self):
        # Matrix A (edit these entries to explore other maps)
        A = [[2, 0],
             [0, -1]]  # scale x by 2, flip y

        title = Text("Grid deformation: T(x) = A x", font_size=36)
        title.to_edge(UP)
        self.play(FadeIn(title))

        plane = NumberPlane(
            x_range=[-4, 4, 1],
            y_range=[-3, 3, 1],
            background_line_style={"stroke_opacity": 0.3}
        )
        self.play(Create(plane))

        # Basis vectors before transform
        e1 = arrow_vec(ORIGIN, RIGHT, color=YELLOW)
        e2 = arrow_vec(ORIGIN, UP, color=GREEN)
        e1_label = MathTex("e_1").next_to(e1.get_end(), DOWN)
        e2_label = MathTex("e_2").next_to(e2.get_end(), LEFT)
        self.play(GrowArrow(e1), GrowArrow(e2), Write(e1_label), Write(e2_label))

        # Show matrix A
        A_tex = matrix_tex(A, label="A")
        A_tex.to_corner(UR).shift(LEFT*0.5 + DOWN*0.2)
        self.play(Write(A_tex))

        # Apply the linear map to the plane and the basis vectors
        self.wait(0.5)
        self.play(
            plane.animate.apply_matrix(A),
            e1.animate.apply_matrix(A),
            e2.animate.apply_matrix(A),
            run_time=2
        )

        # Update basis labels to match new positions
        self.play(
            e1_label.animate.next_to(e1.get_end(), DOWN),
            e2_label.animate.next_to(e2.get_end(), LEFT)
        )
        self.wait(1.5)

# -----------------------------
# 2) Column Picture: Ax as linear combo of columns
# -----------------------------
class ColumnPictureScene(Scene):
    """
    Show Ax = x1 a1 + x2 a2 in R^2 by constructing Ax from the columns of A.
    """
    def construct(self):
        A = [[1, 2],
             [0, 1]]  # shear + shift of second column

        a1 = np.array([A[0][0], A[1][0], 0.0])
        a2 = np.array([A[0][1], A[1][1], 0.0])

        # Choose an x
        x = np.array([1.5, -1.0, 0.0])
        Ax = a1 * x[0] + a2 * x[1]

        title = Text("Column picture: Ax = x1 a1 + x2 a2", font_size=36)
        title.to_edge(UP)
        self.play(FadeIn(title))

        axes = Axes(x_range=[-4, 4, 1], y_range=[-3, 3, 1])
        self.play(Create(axes))

        # Draw columns a1, a2
        a1_arrow = Arrow(ORIGIN, a1, buff=0, color=YELLOW)
        a2_arrow = Arrow(ORIGIN, a2, buff=0, color=GREEN)
        a1_label = MathTex("a_1").next_to(a1_arrow.get_end(), UP)
        a2_label = MathTex("a_2").next_to(a2_arrow.get_end(), RIGHT)
        self.play(GrowArrow(a1_arrow), Write(a1_label))
        self.play(GrowArrow(a2_arrow), Write(a2_label))

        # Scale columns by x1, x2 and show partial sums
        x1_tex = MathTex("x_1 = " + f"{x[0]:.1f}")
        x2_tex = MathTex("x_2 = " + f"{x[1]:.1f}")
        x1_tex.to_corner(UL).shift(DOWN*0.5 + RIGHT*0.5)
        x2_tex.next_to(x1_tex, DOWN)
        self.play(Write(x1_tex), Write(x2_tex))

        x1a1 = Arrow(ORIGIN, a1 * x[0], buff=0, color=YELLOW)
        self.play(GrowArrow(x1a1))
        partial = a1 * x[0]

        x2a2_from_partial = Arrow(partial, partial + a2 * x[1], buff=0, color=GREEN)
        self.play(GrowArrow(x2a2_from_partial))

        # Resulting Ax
        Ax_arrow = Arrow(ORIGIN, Ax, buff=0, color=RED)
        Ax_label = MathTex("A x").next_to(Ax, RIGHT)
        self.play(GrowArrow(Ax_arrow), Write(Ax_label))
        self.wait(1.5)

# -----------------------------
# 3) Composition: C = A B as successive transforms
# -----------------------------
class CompositionScene(Scene):
    """
    Visualize composition T2(T1(x)) as AB acting on a shape and the grid.
    """
    def construct(self):
        # First transform B: rotation by 45 deg
        theta = PI/4
        B = [[np.cos(theta), -np.sin(theta)],
             [np.sin(theta),  np.cos(theta)]]

        # Second transform A: anisotropic scaling
        A = [[1.5, 0.0],
             [0.0, 0.75]]

        title = Text("Composition: C = A B (apply B, then A)", font_size=36)
        title.to_edge(UP)
        self.play(FadeIn(title))

        plane = NumberPlane(
            x_range=[-4, 4, 1],
            y_range=[-3, 3, 1],
            background_line_style={"stroke_opacity": 0.3}
        )
        self.play(Create(plane))

        # A polygon to watch deformation
        poly = Polygon(LEFT, RIGHT, RIGHT+UP, LEFT+UP, color=YELLOW)
        self.play(Create(poly))

        # Show matrices
        B_tex = matrix_tex(B, label="B").to_corner(UR).shift(DOWN*0.5)
        A_tex = matrix_tex(A, label="A").next_to(B_tex, DOWN, aligned_edge=RIGHT)
        self.play(Write(B_tex), Write(A_tex))

        # Apply B then A
        self.play(
            plane.animate.apply_matrix(B),
            poly.animate.apply_matrix(B),
            run_time=1.8
        )
        self.play(
            plane.animate.apply_matrix(A),
            poly.animate.apply_matrix(A),
            run_time=1.8
        )
        self.wait(1.0)

# -----------------------------
# 4) Change of Basis demo
# -----------------------------
class ChangeOfBasisScene(Scene):
    """
    Show a new basis V, coordinates alpha = V^{-1} x, and similarity A_std = V A_new V^{-1}.
    """
    def construct(self):
        # New basis vectors (columns of V)
        v1 = np.array([1.0, 1.0, 0.0])
        v2 = np.array([-1.0, 1.0, 0.0])
        V = np.array([[1.0, -1.0],
                      [1.0,  1.0]])

        # A_new acts in the v-basis (example: simple scaling)
        A_new = np.array([[2.0, 0.0],
                          [0.0, 0.5]])

        # Compute standard basis representation: A_std = V A_new V^{-1}
        Vinv = np.linalg.inv(V)
        A_std = V @ A_new @ Vinv

        title = Text("Change of basis: A_std = V A_new V^{-1}", font_size=36)
        title.to_edge(UP)
        self.play(FadeIn(title))

        axes = Axes(x_range=[-4, 4, 1], y_range=[-3, 3, 1])
        self.play(Create(axes))

        # Draw e1, e2
        e1 = Arrow(ORIGIN, RIGHT, buff=0, color=BLUE)
        e2 = Arrow(ORIGIN, UP, buff=0, color=BLUE)
        self.play(GrowArrow(e1), GrowArrow(e2))

        # Draw v1, v2
        v1_arrow = Arrow(ORIGIN, v1, buff=0, color=GREEN)
        v2_arrow = Arrow(ORIGIN, v2, buff=0, color=GREEN)
        v1_label = MathTex("v_1").next_to(v1, RIGHT)
        v2_label = MathTex("v_2").next_to(v2, UP)
        self.play(GrowArrow(v1_arrow), Write(v1_label))
        self.play(GrowArrow(v2_arrow), Write(v2_label))

        # Show matrices
        V_tex = MathTex("V = ",
                        "\\begin{bmatrix}"
                        f"{V[0,0]} & {V[0,1]} \\\\ {V[1,0]} & {V[1,1]}"
                        "\\end{bmatrix}")
        Anew_tex = MathTex("A_{new} = ",
                           "\\begin{bmatrix}"
                           f"{A_new[0,0]} & {A_new[0,1]} \\\\ {A_new[1,0]} & {A_new[1,1]}"
                           "\\end{bmatrix}")
        V_tex.to_corner(UR).shift(DOWN*0.25)
        Anew_tex.next_to(V_tex, DOWN, aligned_edge=RIGHT)
        self.play(Write(V_tex), Write(Anew_tex))

        # Show an x and its coordinates alpha = V^{-1} x
        x = np.array([2.0, 1.0, 0.0])
        x_arrow = Arrow(ORIGIN, x, buff=0, color=YELLOW)
        x_label = MathTex("x").next_to(x, RIGHT)
        self.play(GrowArrow(x_arrow), Write(x_label))

        alpha = Vinv @ x[:2]
        alpha_tex = MathTex("\\alpha = V^{-1} x = ",
                            "\\begin{bmatrix}"
                            f"{alpha[0]:.1f} \\ {alpha[1]:.1f}"
                            "\\end{bmatrix}")
        alpha_tex.to_corner(UL).shift(DOWN*0.25)
        self.play(Write(alpha_tex))

        # Apply A_new in v-basis: x_new -> y_new, then convert back
        y_new = A_new @ alpha
        y_std = V @ y_new
        y_arrow = Arrow(ORIGIN, np.array([y_std[0], y_std[1], 0.0]), buff=0, color=RED)
        y_label = MathTex("y = A_{std} x").next_to(y_arrow.get_end(), RIGHT)

        # Animate the effect in two conceptual steps
        self.play(TransformFromCopy(x_arrow, y_arrow), run_time=1.5)
        self.play(Write(y_label))
        self.wait(1.0)

# -----------------------------
# 5) Orthogonal projection example
# -----------------------------
class ProjectionScene(Scene):
    """
    Orthogonal projection onto the x-axis in R^2, represented by matrix P.
    """
    def construct(self):
        P = [[1, 0],
             [0, 0]]  # projection onto x-axis

        title = Text("Projection: P onto x-axis", font_size=36)
        title.to_edge(UP)
        self.play(FadeIn(title))

        plane = NumberPlane(
            x_range=[-4, 4, 1],
            y_range=[-3, 3, 1],
            background_line_style={"stroke_opacity": 0.3}
        )
        self.play(Create(plane))

        # Random vector x
        x = np.array([2.5, 1.5, 0.0])
        x_arrow = Arrow(ORIGIN, x, buff=0, color=YELLOW)
        x_label = MathTex("x").next_to(x, RIGHT)
        self.play(GrowArrow(x_arrow), Write(x_label))

        # Its projection Px
        Px = np.array([x[0], 0.0, 0.0])
        Px_arrow = Arrow(ORIGIN, Px, buff=0, color=RED)
        Px_label = MathTex("P x").next_to(Px, DOWN)
        self.play(GrowArrow(Px_arrow), Write(Px_label))

        # Draw dashed line from x to its projection
        drop = DashedLine(x, Px, color=WHITE)
        self.play(Create(drop))
        self.wait(1.0)
