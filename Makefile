book:
	jupyter-book build ./

clean:
	python scripts/clean.py

build:
	jupyter-book build --all ./

site: build
	ghp-import -n -p -f _build/html