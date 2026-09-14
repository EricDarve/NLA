# Incremental HTML build into _build/html (fast).
# Use `make build` instead after editing _toc.yml or _config.yml, or when
# cross-references change, since this target reuses cached pages.
book:
	jupyter-book build ./

# Full rebuild, ignoring all caches. This is what `make site` publishes.
build:
	jupyter-book build --all ./

# Delete everything under _build/.
clean:
	jupyter-book clean --all ./

# Single-file PDF into _build/pdf/book.pdf.
# Requires playwright plus its browser: playwright install chromium
pdf:
	jupyter-book build ./ --builder pdfhtml

# Full rebuild, then publish _build/html to the gh-pages branch on origin.
# Note: this publishes the built directory, not your last commit. Commit first.
site: build
	ghp-import -n -p -f _build/html

.PHONY: book build clean pdf site
