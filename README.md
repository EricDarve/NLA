# CME 302 Numerical Linear Algebra 2026

This GitHub repo contains notes for the course CME302 Numerical Linear Algebra at Stanford University, Fall 2026. The course is taught by [Eric Darve](https://profiles.stanford.edu/eric-darve).

[Read the course notes](https://ericdarve.github.io/NLA/)

These notes are written in MyST Markdown and built with [Jupyter Book 1.x](https://jupyterbook.org/v1/intro.html). The course textbook is [Numerical Linear Algebra with Julia](https://epubs.siam.org/doi/book/10.1137/1.9781611976557), by Eric Darve and Mary Wootters, also available through [Amazon](https://www.amazon.com/Numerical-Linear-Algebra-Julia-Darve/dp/1611976545) and [Google Books](https://play.google.com/store/books/details/Numerical_Linear_Algebra_with_Julia?id=lt9BEAAAQBAJ).

## Setup

From the repository root, create a local Python build environment:

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
```

The Makefile uses the tools in `.venv/bin` directly; activating the environment is optional. This avoids accidentally using an older global `jupyter-book` installation. The book requires Jupyter Book 1.x, as pinned in `requirements.txt`. To use an environment elsewhere, pass `VENV=/path/to/environment` to `make`.

## Building and publishing

All commands are wrapped in the `Makefile`:

| Command | What it does |
| --- | --- |
| `make book` | Incremental HTML build into `_build/html` |
| `make build` | Full rebuild, ignoring caches |
| `make clean` | Delete everything under `_build/` |
| `make pdf` | Single-file PDF into `_build/pdf/book.pdf` |
| `make site` | Full rebuild, then publish to the `gh-pages` branch |

To edit and preview a page, run `make book` and open `_build/html/index.html`. Use `make build` instead of `make book` after editing `_toc.yml` or `_config.yml`, or when cross-references change, since the incremental build reuses cached pages and can leave stale output.

Some Markdown pages contain executable Python cells. The current `_config.yml` uses automatic notebook execution with a 600-second timeout per cell and stops the build if a cell raises an error. The required Python packages are included in `requirements.txt`.

To publish the site:

```
$ make site
```

This rebuilds from scratch and force-pushes `_build/html` to the `gh-pages` branch, which GitHub Pages serves. Two things to keep in mind:

- It publishes the *built directory*, not your last commit, so commit your source changes first to avoid `main` drifting from the live site.
- The `-n` flag in the `ghp-import` call writes a `.nojekyll` file. This is required: without it GitHub Pages runs Jekyll, which ignores the `_static/` and `_sources/` directories and the site loses all of its styling.

`make pdf` uses the `pdfhtml` builder, which renders the book to HTML and then prints it from a headless Chromium driven by [Playwright](https://playwright.dev/python/). Playwright is installed by `requirements.txt`, but its browser has to be downloaded once:

```
$ .venv/bin/playwright install chromium
```

Equations occasionally render poorly through this path. `.venv/bin/jupyter-book build ./ --builder pdflatex` typesets the math properly but needs a LaTeX installation.

## Repo layout

| Path | Contents |
| --- | --- |
| `content/` | Book pages in MyST Markdown, including executable Python examples and images in `content/images/` |
| `content/LICENSE.md` | License for the course notes |
| `_toc.yml` | Table of contents, which sets the chapter order |
| `_config.yml` | Book title, execution, and Sphinx settings |
| `requirements.txt` | Python build dependencies, including the Jupyter Book 1.x constraint |
| `Makefile` | Build, clean, PDF, and publishing commands |
| `_static/` | Files copied verbatim into the built site |
| `manimations/` | Manim sources for the animations, rendered separately; see `manimations/README.md` |
| `_bibliography/` | BibTeX references, not currently cited by any page |
