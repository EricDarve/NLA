from manim import *

class GramSchmidtProcess2D(Scene):
    """
    A Manim scene to visualize the Gram-Schmidt process for two 2D vectors.
    It demonstrates how to construct an orthogonal basis {u1, u2} from
    a set of linearly independent vectors {v1, v2}.
    """
    def construct(self):
        # --- 1. Setup the Scene ---
        x_range_len = 8  # 7 - (-1)
        y_range_len = 6  # 5 - (-1)
        
        axes = Axes(
            x_range=[-1, 7, 1],
            y_range=[-1, 5, 1],
            x_length=x_range_len,
            y_length=y_range_len,
            axis_config={"color": GRAY},
        ).to_edge(LEFT, buff=0.5)

        axes_labels = axes.get_axis_labels(x_label="x_1", y_label="x_2")
        self.add(axes, axes_labels)

        title = Tex("Gram-Schmidt Process").to_edge(UP)
        self.play(Write(title))

        # --- 2. Define and Show Initial Vectors ---
        v1_coords = np.array([4, 1, 0])
        v2_coords = np.array([2, 3, 0])

        v1 = Arrow(axes.c2p(0,0), axes.c2p(*v1_coords[:2]), buff=0, color=BLUE)
        v1_label = MathTex("v_1", color=BLUE).next_to(v1.get_end(), RIGHT)

        v2 = Arrow(axes.c2p(0,0), axes.c2p(*v2_coords[:2]), buff=0, color=RED)
        v2_label = MathTex("v_2", color=RED).next_to(v2.get_end(), UP)

        self.play(
            LaggedStart(
                GrowArrow(v1), Write(v1_label),
                GrowArrow(v2), Write(v2_label),
                lag_ratio=0.5
            )
        )
        self.wait(1)

        # --- 3. Step 1: Normalize v1 to get q1 ---
        step1_text = MathTex(r"1.\; q_1 = \frac{v_1}{\|v_1\|}").to_edge(RIGHT, buff=1).shift(2 * LEFT + UP*2)
        
        # Normalize v1
        q1_coords = v1_coords / np.linalg.norm(v1_coords)
        q1 = Arrow(axes.c2p(0,0), axes.c2p(*q1_coords[:2]), buff=0, color=GREEN)
        q1_label = MathTex("q_1", color=GREEN).next_to(q1.get_end(), 1.6 * DOWN)
        
        self.play(Write(step1_text))
        # Animate v1 shrinking to become the unit vector q1
        self.play(Transform(v1, q1), Transform(v1_label, q1_label))
        self.wait(1)

        # --- 4. Step 2: Get q2 from v2 ---
        step2_text = MathTex(r"2.\; u_2 = v_2 - (v_2 \cdot q_1)q_1").next_to(step1_text, DOWN, buff=0.5, aligned_edge=LEFT)
        step3_text = MathTex(r"3.\; q_2 = \frac{u_2}{\|u_2\|}").next_to(step2_text, DOWN, buff=0.5, aligned_edge=LEFT)

        # Calculate and visualize the projection onto the unit vector q1
        proj_scalar = np.dot(v2_coords, q1_coords)
        proj_coords = proj_scalar * q1_coords
        
        proj_vec = Arrow(axes.c2p(0,0), axes.c2p(*proj_coords[:2]), buff=0, color=YELLOW)
        proj_label = MathTex(r"(v_2 \cdot q_1)q_1", color=YELLOW).next_to(proj_vec.get_end(), RIGHT, buff=0.2)
        
        # Dashed line showing the perpendicular component
        dashed_line = DashedLine(start=axes.c2p(*v2_coords[:2]), end=axes.c2p(*proj_coords[:2]), color=WHITE)
        
        self.play(Write(step2_text))
        self.play(GrowArrow(proj_vec), Write(proj_label))
        self.play(Create(dashed_line))
        self.wait(1)

        # Visualize the vector subtraction to get u2
        u2_coords = v2_coords - proj_coords
        u2 = Arrow(axes.c2p(0,0), axes.c2p(*u2_coords[:2]), buff=0, color=PURPLE)
        u2_label = MathTex("u_2", color=PURPLE).next_to(u2.get_end(), LEFT)
        
        # Animate the subtraction vector moving to the origin
        subtraction_arrow = Arrow(
            start=axes.c2p(*proj_coords[:2]),
            end=axes.c2p(*v2_coords[:2]),
            buff=0,
            color=PURPLE
        )
        self.play(Transform(subtraction_arrow, u2))
        self.play(Write(u2_label))
        self.wait(1)

        # Normalize u2 to get q2
        self.play(Write(step3_text))
        q2_coords = u2_coords / np.linalg.norm(u2_coords)
        q2 = Arrow(axes.c2p(0,0), axes.c2p(*q2_coords[:2]), buff=0, color=ORANGE)
        q2_label = MathTex("q_2", color=ORANGE).next_to(q2.get_end(), LEFT)
        
        self.play(Transform(u2, q2), Transform(u2_label, q2_label))
        self.wait(1)

        # --- 5. Show the Final Orthonormal Basis ---
        final_text = Tex("Final Orthonormal Basis: $\{q_1, q_2\}$", color=GOLD).next_to(title, DOWN)
        right_angle = RightAngle(Line(axes.c2p(0,0), q1.get_end()), Line(axes.c2p(0,0), q2.get_end()), length=0.4, color=WHITE)

        self.play(
            FadeOut(v2, v2_label, proj_vec, proj_label, dashed_line, subtraction_arrow, step1_text, step2_text, step3_text),
            FadeIn(final_text),
            Create(right_angle)
        )
        self.wait(3)