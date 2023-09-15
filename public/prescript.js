(function () {// quartz/components/scripts/quartz/components/scripts/darkmode.inline.ts
var userPref = window.matchMedia("(prefers-color-scheme: light)").matches ? "light" : "dark";
var currentTheme = localStorage.getItem("theme") ?? userPref;
document.documentElement.setAttribute("saved-theme", currentTheme);
document.addEventListener("nav", () => {
  const switchTheme = (e) => {
    if (e.target.checked) {
      document.documentElement.setAttribute("saved-theme", "dark");
      localStorage.setItem("theme", "dark");
    } else {
      document.documentElement.setAttribute("saved-theme", "light");
      localStorage.setItem("theme", "light");
    }
  };
  const toggleSwitch = document.querySelector("#darkmode-toggle");
  toggleSwitch.removeEventListener("change", switchTheme);
  toggleSwitch.addEventListener("change", switchTheme);
  if (currentTheme === "dark") {
    toggleSwitch.checked = true;
  }
});
})();