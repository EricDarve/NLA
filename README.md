# CME 302 Numerical Linear Algebra 2026

This GitHub repo contains the class notes for the Numerical Linear Algebra course at Stanford University, Spring 2026. The course is taught by [Eric Darve](https://profiles.stanford.edu/eric-darve).

[Site location](https://ericdarve.github.io/NLA/)

The book was written using Markdown and [Jupyter Book](https://jupyterbook.org/en/stable/intro.html). The book is available online at [Numerical Linear Algebra with Julia](https://epubs.siam.org/doi/book/10.1137/1.9781611976557), and can be purchased on [Amazon](https://www.amazon.com/Numerical-Linear-Algebra-Julia-Darve/dp/1611976545) or [Google Books](https://play.google.com/store/books/details/Numerical_Linear_Algebra_with_Julia?id=lt9BEAAAQBAJ).

## Setup

Create a local build environment with Python 3.12:

```sh
python3.12 -m venv .venv
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
| `content/` | All book pages, as Markdown |
| `_toc.yml` | Table of contents, which sets the chapter order |
| `_config.yml` | Book title, execution, and Sphinx settings |
| `_static/` | Files copied verbatim into the built site |
| `manimations/` | Manim sources for the animations, rendered separately; see `manimations/README.md` |
| `addl_material/` | Notes not currently listed in `_toc.yml` |
| `_bibliography/` | BibTeX references, not currently cited by any page |
