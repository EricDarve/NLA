# CME 302 Class Notes 2024

- [GitHub Pages web site](https://ericdarve.github.io/NLA/)
- [CME 302](https://canvas.stanford.edu/courses/178134)
- Built using [Obsidian](https://obsidian.md/) and [Quartz](https://quartz.jzhao.xyz/).

### How to generate the static web pages

```
npx quartz build --serve
```

This will start a local web server to run your Quartz on your computer. Open a web browser and visit `http://localhost:8080/` to view it.

Note: update
```
export const QUARTZ_SOURCE_BRANCH = "v4"
```

### Sync with GitHub

```
npx quartz sync
```
