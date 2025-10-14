from manim import *
import numpy as np

class GivensRotation(Scene):
    """
    A Manim scene to visualize a Givens rotation applied to a 2D vector.
    The rotation aligns the vector with the positive x-axis.
    """
    def construct(self):
        # Set up the coordinate system
        axes = Axes(
            x_range=[-1, 6, 1],
            y_range=[-1, 5, 1],
            axis_config={"color": GRAY},
            x_length=7,
            y_length=6,
        )
        axes_labels = axes.get_axis_labels(x_label="x_1", y_label="x_2")
        self.add(axes, axes_labels)

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

        # Define the target vector (rotated vector)
        norm_x = np.linalg.norm(x_vec)
        target_vec_coords = np.array([norm_x, 0, 0])
        target_label = MathTex(r"Gx = \begin{pmatrix} r \\ 0 \end{pmatrix}", color=GREEN)
        target_label.next_to(axes.c2p(norm_x, 0), DOWN, buff=0.4)
        
        # Calculate the rotation angle
        # The angle is the negative of the vector's angle with the x-axis
        rotation_angle = -np.arctan2(x_vec[1], x_vec[0])

        # Create the rotation arc, ensuring it's centered on the axes' origin
        arc = Arc(
            radius=1.5,
            angle=rotation_angle,
            start_angle=np.arctan2(x_vec[1], x_vec[0]),
            color=RED,
            arc_center=axes.c2p(0, 0)
        )
        
        # Position the label more robustly along the arc
        label_pos = arc.point_from_proportion(0.5)
        theta_label = MathTex(r"-\theta", color=RED).move_to(label_pos).shift(0.2*UP + 0.4*RIGHT)
        
        # --- Animation Sequence ---
        
        # 1. Introduce the initial vector x
        self.play(
            GrowArrow(x_arrow),
            Write(x_label)
        )
        self.wait(1)

        # 2. Show the rotation arc and angle
        self.play(
            Create(arc),
            Write(theta_label)
        )
        self.wait(1)

        # 3. Perform the rotation
        self.play(
            Rotate(
                x_arrow,
                angle=rotation_angle,
                about_point=axes.c2p(0, 0)
            ),
            # Move the label with the arrow tip
            x_label.animate.next_to(axes.c2p(norm_x, 0), UR, buff=0.2)
        )
        
        # 4. Show the label for the final vector
        self.play(Write(target_label))
        self.wait(3)