// Inline the engine and data into app.html to make a single shareable file.
//   node interactive/build.js [fragment-output.html]
// Writes docs/nursery_trap.html (a complete standalone page that works
// offline and can be emailed; docs/ is served by GitHub Pages). If a second
// path is given, also writes the page body alone there, for hosts that
// supply their own <html> skeleton.
"use strict";
const fs = require("fs");
const path = require("path");

const dir = __dirname;
const src = fs.readFileSync(path.join(dir, "app.html"), "utf8");
const inlined = src.replace(/<script src="([^"]+)"><\/script>/g, (_, file) => {
  const code = fs.readFileSync(path.join(dir, file), "utf8");
  if (code.includes("</script")) throw new Error(`${file} contains </script>`);
  return `<script>\n${code}\n</script>`;
});

const standalone = `<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">
</head>
<body>
${inlined}
</body>
</html>
`;
fs.writeFileSync(path.join(dir, "..", "docs", "nursery_trap.html"), standalone);
console.log(`Wrote docs/nursery_trap.html (${(standalone.length / 1024).toFixed(0)} KB)`);
if (process.argv[2]) {
  fs.writeFileSync(process.argv[2], inlined);
  console.log(`Wrote ${process.argv[2]}`);
}
