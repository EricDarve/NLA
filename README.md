# CME 302 Class Notes 2024

- [GitHub Pages web site](https://ericdarve.github.io/NLA/)
- Built using [Obsidian](https://obsidian.md/) and [Quartz](https://quartz.jzhao.xyz/), v4.3.1, Quartz [github](https://github.com/jackyzha0/quartz) repo.

### How to generate the static web pages

```
npx quartz build --serve
```

This will start a local web server to run your Quartz on your computer. Open a web browser and visit `http://localhost:8080/` to view it.

Note: update file `constants.js` (`v4` &#8594; `main`):
```
export const QUARTZ_SOURCE_BRANCH = "v4"
```

### Sync with GitHub

```
npx quartz sync
```
