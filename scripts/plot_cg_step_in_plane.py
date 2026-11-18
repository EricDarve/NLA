"""Render the CG energy slice with LaTeX, TikZ, and pdftocairo."""
from pathlib import Path
from subprocess import run
from tempfile import TemporaryDirectory

OUTPUT = Path(__file__).resolve().parents[1] / "_static" / "cg_step_in_plane.svg"

FIGURE = r"""
\documentclass[tikz,border=4pt]{standalone}
\usepackage{amsmath,amssymb}
\usetikzlibrary{arrows.meta}
\definecolor{cgblue}{HTML}{205FBA}
\definecolor{ink}{HTML}{252525}
\begin{document}
\begin{tikzpicture}[x=2.35cm,y=2.35cm,text=ink,font=\sffamily,
  vec/.style={-{Stealth[length=2.5mm]},line width=1.25pt},
  labelbg/.style={fill=white,inner sep=1.2pt}]
\node at (1.72,2.09) {\large Energy restricted to $x_k+S,\quad S=\operatorname{span}\{p_k,r_k\}$};
% The slice has quadratic part with eigenvalues 1,3 and center (2,1).
\begin{scope}
\clip (-.43,-.48) rectangle (3.88,1.88);
\foreach \level in {.18,.75,1.5,3,4.5}{
  \pgfmathsetmacro{\major}{sqrt(2*\level)}
  \pgfmathsetmacro{\minor}{sqrt(2*\level/3)}
  \draw[gray!46,rotate around={45:(2,1)}] (2,1) ellipse
    [x radius=\major,y radius=\minor];
}
\draw[gray,rotate around={45:(2,1)}] (2,1) ellipse
  [x radius={sqrt(1.5)},y radius={sqrt(.5)}];
\end{scope}
\draw[gray,dashed] (0,-.33)--(0,.44) (-.2,0)--(3.25,0);
\fill[cgblue!19,opacity=.8] (1.5,0)--(1.5,1.5)--(2.25,1.5)--cycle;
\draw[ink,thin] (1.43,0)--(1.43,.07)--(1.5,.07);
\draw[ink,thin] (.07,0)--(.07,.07)--(0,.07);
\draw[vec,cgblue] (0,0)--(1.5,0) node[midway,below=4pt] {$\mu_k p_k$};
% The short black arrow indicates the direction of the projected residual.
\draw[vec,ink] (0,0)--(.50,0) node[midway,above=4pt,labelbg] {$\operatorname{proj}_S r_{k-1}$};
\draw[vec,ink] (1.5,0)--(1.5,1.5) node[midway,left] {$r_k$};
\draw[vec,cgblue] (1.5,1.5)--(2.25,1.5) node[midway,above] {$\tau_k p_k$};
\draw[vec,cgblue] (1.5,0)--(2.25,1.5) node[right=2pt,labelbg] {$p_{k+1}=r_k+\tau_kp_k$};
\foreach \x/\y in {0/0,1.5/0,2/1}{\fill[ink] (\x,\y) circle[radius=1.6pt];}
\node[below left] at (0,0) {$x_{k-1}$};
\node[below] at (1.5,0) {$x_k$};
\node[right=3pt,labelbg] at (2,1) {$x_{k+1}$};
\node at (1.72,-.57) {The center $x_{k+1}=x_k+\mu_{k+1}p_{k+1}$ minimizes the energy in this plane.};
\draw[gray!40] (-.42,-.72)--(3.86,-.72);
\node at (1.72,-.91) {Residual orthogonality:\quad $r_i^Tr_j=0\ (i\ne j),\qquad r_i^Tp_j=0\ (1\le j\le i)$.};
\node at (1.72,-1.17) {Error orthogonality:\quad $(x-x_i)^TAp_j=0\ (1\le j\le i)$.};
\node at (1.72,-1.43) {Conjugate directions:\quad $p_i^TAp_j=0\ (i\ne j)$.};
\end{tikzpicture}
\end{document}
"""

if __name__ == "__main__":
    with TemporaryDirectory() as directory:
        tex = Path(directory) / "figure.tex"
        tex.write_text(FIGURE)
        run(["pdflatex", "-halt-on-error", "-interaction=nonstopmode", tex.name],
            cwd=directory, check=True)
        run(["pdftocairo", "-svg", "figure.pdf", str(OUTPUT)],
            cwd=directory, check=True)
