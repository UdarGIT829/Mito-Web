# Standalone mitochondrial genome visualizer

Open `index.html` in a browser. The visualizer has no package, server, or
application dependencies.

The ring is generated from pinned rCRS / NC_012920.1 coordinates in
`mitochondrial-genome.js`. Features are SVG paths, so later selection and
mutation-position behavior can be added without changing the base rendering.

To embed the prototype elsewhere:

```html
<link rel="stylesheet" href="./mitochondrial-genome.css">
<div id="mitochondrial-genome"></div>
<script src="./mitochondrial-genome.js"></script>
<script>
  MitochondrialGenomeVisualizer.render(
    document.querySelector("#mitochondrial-genome")
  );
</script>
```

Pass `{ showLabels: false }` as the second argument to omit gene labels.

The standalone page enables an optional checkbox panel with
`{ debugControls: true }`. Leave that option out (or set it to `false`) when
embedding the visualizer in the application. A starting selection can be
provided with:

```js
MitochondrialGenomeVisualizer.render(container, {
  selectedGenes: ["MT-ND1", "MT-ND2"]
});
```

Selections can also be updated without rerendering:

```js
MitochondrialGenomeVisualizer.applySelection(container, ["MT-CO1"]);
```

The mode box in the center of the ring is part of the visualizer and is always
rendered. Its two direct buttons switch between the expected position-driven
workflow and a genes-only workflow. Genes-only mode hides mutation markers
while retaining the current gene selection. `debugControls` affects only the
development panel below the visualizer.

Exact rCRS positions can be supplied during rendering or updated later:

```js
MitochondrialGenomeVisualizer.render(container, {
  mutations: [73, 9207, 10760]
});

MitochondrialGenomeVisualizer.applyMutations(container, [73, 9207]);
```

Markers are drawn at their proportional 1–16,569 coordinate. If a position
overlaps genes on both visual tracks, both tracks receive a marker.

All mitochondrial tRNA features participate in labeling, direct selection, and
position-driven selection. Their display names use the conventional `TRN`
prefix (for example, `MT-TP` is displayed as `TRNP`).

Labels use four ordered billboard lanes around a reduced-radius genome ring.
Leader lines preserve the exact feature association while keeping dense tRNA
clusters and adjacent selected labels from occupying the same coordinates.

The expected application workflow can provide selected mutation positions
directly. The visualizer places their markers and derives the selected genes:

```js
MitochondrialGenomeVisualizer.render(container, {
  selectedMutations: [9207, 10760]
});

const result = MitochondrialGenomeVisualizer.applySelectedMutations(
  container,
  [9207, 10760]
);
// result.genes includes every feature containing either coordinate.
```
