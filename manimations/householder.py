from manim import *

class HouseholderVectorConstruction(Scene):
    """
    A Manim scene to visualize the construction of the vector 'v'
    used in a Householder transformation.
    v = x - ||x|| * e_1
    """
    def construct(self):
        # Set up the coordinate system
        axes = Axes(
            x_range=[-1, 6, 1],
            y_range=[-1, 4, 1],
            axis_config={"color": GRAY},
            x_length=7,
            y_length=5,
        )
        axes_labels = axes.get_axis_labels(x_label="x_1", y_label="x_2")

        # Define the initial vector 'x'
        x_vec = np.array([4, 3, 0])
        x_arrow = Arrow(
            start=axes.c2p(0, 0),
            end=axes.c2p(4, 3),
            buff=0,
            color=BLUE,
            stroke_width=7
        )
        x_label = MathTex("x", color=BLUE).next_to(x_arrow.get_end(), UR, buff=0.2)

        # Define the target vector ||x|| * e_1
        norm_x = np.linalg.norm(x_vec)
        target_vec = np.array([norm_x, 0, 0])
        target_arrow = Arrow(
            start=axes.c2p(0, 0),
            end=axes.c2p(norm_x, 0),
            buff=0,
            color=GREEN,
            stroke_width=7
        )
        target_label = MathTex(r"\|x\|_2 e_1", color=GREEN).next_to(target_arrow.get_end(), DOWN, buff=0.2)

        # Define the Householder vector v = x - ||x||*e_1
        v_vec = x_vec - target_vec
        v_arrow = Arrow(
            start=axes.c2p(norm_x, 0),
            end=axes.c2p(4, 3),
            buff=0,
            color=RED,
            stroke_width=7
        )
        v_label = MathTex("v = x - \|x\|_2 e_1", color=RED).next_to(v_arrow, RIGHT, buff=0.2)
        v_label.shift(DOWN * 0.5)
        
        # Define the reflection line (hyperplane)
        # This line is perpendicular to v and passes through the midpoint of x and target.
        midpoint = (x_vec + target_vec) / 2
        # The direction of the line is orthogonal to v
        line_direction = np.array([-v_vec[1], v_vec[0], 0])
        
        reflection_line = DashedLine(
            start=axes.c2p(*(midpoint - line_direction * 3)),
            end=axes.c2p(*(midpoint + line_direction * 3)),
            color=YELLOW,
            stroke_width=5
        )
        reflection_label = Text("Reflection Hyperplane", color=YELLOW, font_size=24).next_to(reflection_line, UL, buff=0.1)

        # --- Animation Sequence ---
        self.play(Write(axes), Write(axes_labels))
        
        # 1. Introduce the vector x
        self.play(
            LaggedStart(
                GrowArrow(x_arrow),
                Write(x_label),
                lag_ratio=0.5
            )
        )
        self.wait(1)

        # 2. Show its projection onto the e_1 axis
        self.play(
            LaggedStart(
                GrowArrow(target_arrow),
                Write(target_label),
                lag_ratio=0.5
            )
        )
        self.wait(1)

        # 3. Construct the vector v
        self.play(
            LaggedStart(
                GrowArrow(v_arrow),
                Write(v_label),
                lag_ratio=0.5
            )
        )
        self.wait(1)

        # 4. Show the line of reflection
        self.play(
            Create(reflection_line),
            Write(reflection_label)
        )
        self.wait(3)