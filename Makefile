VENV ?= .venv
JUPYTER_BOOK ?= $(VENV)/bin/jupyter-book
GHP_IMPORT ?= $(VENV)/bin/ghp-import
export PATH := $(abspath $(VENV))/bin:$(PATH)

# Incremental HTML build into _build/html (fast).
# Use `make build` instead after editing _toc.yml or _config.yml, or when
# cross-references change, since this target reuses cached pages.
book: check-book
	"$(JUPYTER_BOOK)" build ./

# Full rebuild, ignoring all caches. This is what `make site` publishes.
build: check-book
	"$(JUPYTER_BOOK)" build --all ./

# Delete everything under _build/.
clean: check-book
	"$(JUPYTER_BOOK)" clean --all ./

# Single-file PDF into _build/pdf/book.pdf.
# Requires playwright plus its browser: playwright install chromium
pdf: check-book
	"$(JUPYTER_BOOK)" build ./ --builder pdfhtml

# Full rebuild, then publish _build/html to the gh-pages branch on origin.
# Note: this publishes the built directory, not your last commit. Commit first.
site: build
	"$(GHP_IMPORT)" -n -p -f _build/html

check-book:
	@test -x "$(JUPYTER_BOOK)" || { echo 'Build environment missing. Follow the Setup instructions in README.md.'; exit 1; }

.PHONY: book build clean pdf site check-book
