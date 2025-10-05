from manim import *
import numpy as np

class EigenvaluesSingularValues(Scene):
    def construct(self):
        # --- Parameters ---
        epsilon = 0.1 # A small epsilon to show ill-conditioning
        
        # --- Matrices ---
        X = np.array([[1, 1], [epsilon, -epsilon]])
        D = np.array([[1, 0], [0, -1]])
        X_inv = np.linalg.inv(X)
        A = X @ D @ X_inv
        
        # Manim uses np.array directly, but for display might need to format
        A_display = Tex(
            r"$A = X D X^{-1} = \begin{pmatrix} 0 & \frac{1}{\epsilon} \\ \epsilon & 0 \end{pmatrix}$"
        ).scale(0.8).to_edge(UP + LEFT)
        
        A_val_display = Tex(
            r"$A = \begin{pmatrix} 0 & " + f"{1/epsilon:.2f}" + r" \\ " + f"{epsilon:.2f}" + r" & 0 \end{pmatrix}$"
        ).scale(0.8).next_to(A_display, DOWN)

        A_sq_display = Tex(r"$A^2 = I = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}$").scale(0.8).to_edge(UP + RIGHT)

        # --- Eigen Stuff ---
        eigen_values, eigen_vectors = np.linalg.eig(A)
        # Manim's vector class is 3D, so convert 2D vectors
        eigen_vector1 = eigen_vectors[:, 0]
        eigen_vector2 = eigen_vectors[:, 1]
        
        # Scale for display
        eigen_vector1_scaled = np.append(eigen_vector1 * 2 / np.linalg.norm(eigen_vector1), 0)
        eigen_vector2_scaled = np.append(eigen_vector2 * 2 / np.linalg.norm(eigen_vector2), 0)

        # --- Singular Value Decomposition ---
        U, s, Vh = np.linalg.svd(A)
        singular_values = s
        # Vh is V transpose, so the columns of V are rows of Vh
        singular_vector_v1 = Vh[0, :] 
        singular_vector_v2 = Vh[1, :]
        
        # Scale for display
        singular_vector_v1_scaled = np.append(singular_vector_v1 * 2, 0)
        singular_vector_v2_scaled = np.append(singular_vector_v2 * 2, 0)

        # --- Scene Setup ---
        plane = NumberPlane(
            x_range=[-5, 5, 1],
            y_range=[-5, 5, 1],
            x_length=10,
            y_length=10,
            background_line_style={
                "stroke_color": GREEN_A,
                "stroke_width": 2,
                "stroke_opacity": 0.6
            }
        ).add_coordinates().shift(LEFT * 2)
        
        title = Tex("Eigenvalues vs. Singular Values for $A = \\begin{pmatrix} 0 & 1/\\epsilon \\ \\epsilon & 0 \\end{pmatrix}$", font_size=40).to_edge(UP)

        self.add(plane)
        self.play(Write(title))
        self.wait(0.5)
        self.play(Write(A_display), Write(A_val_display), Write(A_sq_display))
        self.wait(1)

        # --- Part 1: Eigenvectors and Eigenvalues ---
        self.next_section("Eigenvector Demonstration")
        eigen_text = Text("Eigenvectors: Directions unchanged (scaled)", font_size=30).to_edge(DOWN)
        self.play(Write(eigen_text))
        
        # Display eigenvectors
        vec_e1 = Arrow(ORIGIN, eigen_vector1_scaled, buff=0, color=RED)
        label_e1 = MathTex(r"v_1").next_to(vec_e1, RIGHT, buff=0.1).set_color(RED)
        vec_e2 = Arrow(ORIGIN, eigen_vector2_scaled, buff=0, color=BLUE)
        label_e2 = MathTex(r"v_2").next_to(vec_e2, RIGHT, buff=0.1).set_color(BLUE)

        self.play(GrowArrow(vec_e1), Write(label_e1), GrowArrow(vec_e2), Write(label_e2))
        self.wait(1)

        # Apply A once to eigenvectors
        vec_e1_transformed_pos = np.append((A @ eigen_vector1_scaled[:2]), 0)
        vec_e2_transformed_pos = np.append((A @ eigen_vector2_scaled[:2]), 0)
        
        vec_e1_transformed = Arrow(ORIGIN, vec_e1_transformed_pos, buff=0, color=RED)
        vec_e2_transformed = Arrow(ORIGIN, vec_e2_transformed_pos, buff=0, color=BLUE)

        self.play(
            Transform(vec_e1, vec_e1_transformed),
            Transform(vec_e2, vec_e2_transformed)
        )
        self.wait(1)
        self.play(FadeOut(vec_e1, label_e1, vec_e2, label_e2, eigen_text))

        # --- Part 2: Arbitrary Vector and its transformation ---
        self.next_section("Arbitrary Vector Transformation")
        arbitrary_text = Text("Applying A to a general vector", font_size=30).to_edge(DOWN)
        self.play(Write(arbitrary_text))

        initial_vector = np.array([1, 1, 0]) # A simple vector for demonstration
        initial_vector_scaled = initial_vector * 1.5 / np.linalg.norm(initial_vector)
        
        vec_orig = Arrow(ORIGIN, initial_vector_scaled, buff=0, color=WHITE)
        label_orig = MathTex(r"x").next_to(vec_orig, UP).set_color(WHITE)
        self.play(GrowArrow(vec_orig), Write(label_orig))
        self.wait(1)
        
        # Apply A once
        vec_A_once_pos = np.append((A @ initial_vector_scaled[:2]), 0)
        vec_A_once = Arrow(ORIGIN, vec_A_once_pos, buff=0, color=YELLOW)
        label_A_once = MathTex(r"Ax").next_to(vec_A_once, UP).set_color(YELLOW)
        
        self.play(
            Transform(vec_orig, vec_A_once),
            Transform(label_orig, label_A_once)
        )
        self.add_foreground_mobject(label_A_once) # Keep label on top
        self.wait(1)

        # Apply A twice (A^2 = I)
        vec_A_twice_pos = np.append((A @ A @ initial_vector_scaled[:2]), 0) # Should be close to initial_vector_scaled
        vec_A_twice = Arrow(ORIGIN, vec_A_twice_pos, buff=0, color=GREEN)
        label_A_twice = MathTex(r"A^2x = Ix").next_to(vec_A_twice, UP).set_color(GREEN)

        self.play(
            Transform(vec_orig, vec_A_twice),
            Transform(label_orig, label_A_twice)
        )
        self.wait(2)
        self.play(FadeOut(vec_orig, label_orig, arbitrary_text))

        # --- Part 3: Singular Vectors and Singular Values ---
        self.next_section("Singular Vector Demonstration")
        singular_text = Text("Singular Vectors: Directions of Max/Min Stretch", font_size=30).to_edge(DOWN)
        self.play(Write(singular_text))

        # Display singular vectors (input basis)
        vec_sv1_in = Arrow(ORIGIN, singular_vector_v1_scaled, buff=0, color=ORANGE)
        label_sv1_in = MathTex(r"v_1").next_to(vec_sv1_in, RIGHT, buff=0.1).set_color(ORANGE)
        vec_sv2_in = Arrow(ORIGIN, singular_vector_v2_scaled, buff=0, color=PURPLE)
        label_sv2_in = MathTex(r"v_2").next_to(vec_sv2_in, RIGHT, buff=0.1).set_color(PURPLE)

        self.play(GrowArrow(vec_sv1_in), Write(label_sv1_in), GrowArrow(vec_sv2_in), Write(label_sv2_in))
        self.wait(1)
        
        # Transform by A
        vec_sv1_out_pos = np.append((A @ singular_vector_v1_scaled[:2]), 0)
        vec_sv2_out_pos = np.append((A @ singular_vector_v2_scaled[:2]), 0)

        # Calculate expected length for display
        expected_len_sv1 = singular_values[0] * np.linalg.norm(singular_vector_v1_scaled)
        expected_len_sv2 = singular_values[1] * np.linalg.norm(singular_vector_v2_scaled)

        # Ensure the arrow scaling reflects the singular value
        vec_sv1_out = Arrow(ORIGIN, vec_sv1_out_pos, buff=0, color=ORANGE, max_stroke_width_to_length_ratio=0.03)
        vec_sv2_out = Arrow(ORIGIN, vec_sv2_out_pos, buff=0, color=PURPLE, max_stroke_width_to_length_ratio=0.03)

        # Display singular values
        sv_text_large = MathTex(r"\|Av_1\| \approx \sigma_1 \|v_1\| = {:.2f}".format(singular_values[0]), substrings_to_isolate=["\\sigma_1"]).scale(0.7).to_corner(DR).set_color(ORANGE)
        sv_text_small = MathTex(r"\|Av_2\| \approx \sigma_2 \|v_2\| = {:.2f}".format(singular_values[1]), substrings_to_isolate=["\\sigma_2"]).scale(0.7).next_to(sv_text_large, UP).set_color(PURPLE)

        self.play(
            Transform(vec_sv1_in, vec_sv1_out),
            Transform(vec_sv2_in, vec_sv2_out),
            FadeOut(label_sv1_in, label_sv2_in)
        )
        self.wait(1)
        self.play(Write(sv_text_large), Write(sv_text_small))
        self.wait(3)

        self.play(FadeOut(Group(*self.mobjects)))
        self.wait(1)
