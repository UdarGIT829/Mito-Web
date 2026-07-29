(function (global) {
  "use strict";

  const NS = "http://www.w3.org/2000/svg";
  const GENOME_LENGTH = 16569;

  // Pinned rCRS / NC_012920.1 coordinates. This local copy is intentional:
  // the reference visualizer must remain usable outside the larger application.
  const FEATURES = [
    { name: "D-loop", start: 16024, end: 576, type: "control" },
    { name: "MT-TF", start: 577, end: 647, type: "trna" },
    { name: "MT-RNR1", start: 648, end: 1601, type: "rrna" },
    { name: "MT-TV", start: 1602, end: 1670, type: "trna" },
    { name: "MT-RNR2", start: 1671, end: 3229, type: "rrna" },
    { name: "MT-TL1", start: 3230, end: 3304, type: "trna" },
    { name: "MT-ND1", start: 3307, end: 4262, type: "complex-i" },
    { name: "MT-TI", start: 4263, end: 4331, type: "trna" },
    { name: "MT-TQ", start: 4329, end: 4400, type: "trna" },
    { name: "MT-TM", start: 4402, end: 4469, type: "trna" },
    { name: "MT-ND2", start: 4470, end: 5511, type: "complex-i" },
    { name: "MT-TW", start: 5512, end: 5579, type: "trna" },
    { name: "MT-TA", start: 5587, end: 5655, type: "trna" },
    { name: "MT-TN", start: 5657, end: 5729, type: "trna" },
    { name: "MT-TC", start: 5761, end: 5826, type: "trna" },
    { name: "MT-TY", start: 5826, end: 5891, type: "trna" },
    { name: "MT-CO1", start: 5904, end: 7445, type: "complex-iv" },
    { name: "MT-TS1", start: 7446, end: 7514, type: "trna" },
    { name: "MT-TD", start: 7518, end: 7585, type: "trna" },
    { name: "MT-CO2", start: 7586, end: 8269, type: "complex-iv" },
    { name: "MT-TK", start: 8295, end: 8364, type: "trna" },
    { name: "MT-ATP8", start: 8366, end: 8572, type: "complex-v", track: "inner" },
    { name: "MT-ATP6", start: 8527, end: 9207, type: "complex-v" },
    { name: "MT-CO3", start: 9207, end: 9990, type: "complex-iv" },
    { name: "MT-TG", start: 9991, end: 10058, type: "trna" },
    { name: "MT-ND3", start: 10059, end: 10404, type: "complex-i" },
    { name: "MT-TR", start: 10405, end: 10469, type: "trna" },
    { name: "MT-ND4L", start: 10470, end: 10766, type: "complex-i", track: "inner" },
    { name: "MT-ND4", start: 10760, end: 12137, type: "complex-i" },
    { name: "MT-TH", start: 12138, end: 12206, type: "trna" },
    { name: "MT-TS2", start: 12207, end: 12265, type: "trna" },
    { name: "MT-TL2", start: 12266, end: 12336, type: "trna" },
    { name: "MT-ND5", start: 12337, end: 14148, type: "complex-i" },
    { name: "MT-ND6", start: 14149, end: 14673, type: "complex-i" },
    { name: "MT-TE", start: 14674, end: 14742, type: "trna" },
    { name: "MT-CYB", start: 14747, end: 15887, type: "complex-iii" },
    { name: "MT-TT", start: 15888, end: 15953, type: "trna" },
    { name: "MT-TP", start: 15956, end: 16023, type: "trna" }
  ];

  const COLORS = {
    control: "#bcbcbc",
    rrna: "#19f527",
    trna: "#fff",
    "complex-i": "#fff500",
    "complex-iii": "#ff7b4d",
    "complex-iv": "#77b8e5",
    "complex-v": "#ff9f05"
  };

  function svgElement(name, attributes) {
    const element = document.createElementNS(NS, name);
    Object.entries(attributes || {}).forEach(([key, value]) => {
      element.setAttribute(key, value);
    });
    return element;
  }

  function positionToAngle(position) {
    return ((position - 1) / GENOME_LENGTH) * 360 - 90;
  }

  function point(cx, cy, radius, degrees) {
    const radians = degrees * Math.PI / 180;
    return {
      x: cx + radius * Math.cos(radians),
      y: cy + radius * Math.sin(radians)
    };
  }

  function annularArcPath(start, end, outerRadius, innerRadius) {
    const cx = 500;
    const cy = 500;
    let startAngle = positionToAngle(start);
    let endAngle = positionToAngle(end + 1);
    if (end < start) endAngle += 360;

    const span = endAngle - startAngle;
    const largeArc = span > 180 ? 1 : 0;
    const a = point(cx, cy, outerRadius, startAngle);
    const b = point(cx, cy, outerRadius, endAngle);
    const c = point(cx, cy, innerRadius, endAngle);
    const d = point(cx, cy, innerRadius, startAngle);

    return [
      `M ${a.x} ${a.y}`,
      `A ${outerRadius} ${outerRadius} 0 ${largeArc} 1 ${b.x} ${b.y}`,
      `L ${c.x} ${c.y}`,
      `A ${innerRadius} ${innerRadius} 0 ${largeArc} 0 ${d.x} ${d.y}`,
      "Z"
    ].join(" ");
  }

  function midpoint(feature) {
    let end = feature.end;
    if (end < feature.start) end += GENOME_LENGTH;
    return ((feature.start + end) / 2 - 1) % GENOME_LENGTH + 1;
  }

  function featureContainsPosition(feature, position) {
    if (feature.end >= feature.start) {
      return position >= feature.start && position <= feature.end;
    }
    return position >= feature.start || position <= feature.end;
  }

  function addFeature(svg, feature, showLabels) {
    const isInner = feature.track === "inner";
    const outerRadius = isInner ? 395 : 442;
    const innerRadius = isInner ? 366 : 402;
    const path = svgElement("path", {
      class: "mitochondrial-genome__segment",
      d: annularArcPath(feature.start, feature.end, outerRadius, innerRadius),
      fill: COLORS[feature.type],
      "data-feature": feature.name,
      "data-start": feature.start,
      "data-end": feature.end
    });
    path.appendChild(svgElement("title"));
    path.lastChild.textContent =
      `${feature.name}: ${feature.start.toLocaleString()}–${feature.end.toLocaleString()} bp`;
    svg.appendChild(path);

    if (!showLabels || feature.type === "trna") return;

    const angle = positionToAngle(midpoint(feature));
    const labelRadius = isInner ? 337 : 366;
    const labelPoint = point(500, 500, labelRadius, angle);
    const text = svgElement("text", {
      class: `mitochondrial-genome__label${feature.name.length > 7 ? " mitochondrial-genome__label--small" : ""}`,
      x: labelPoint.x,
      y: labelPoint.y,
      dy: ".35em",
      "data-label-for": feature.name
    });
    text.textContent = feature.name.replace("MT-", "");
    svg.appendChild(text);
  }

  function applySelection(container, selectedNames) {
    const selected = new Set(selectedNames);
    container.classList.toggle(
      "mitochondrial-genome--filtering",
      selected.size > 0
    );
    container.querySelectorAll("[data-feature]").forEach((segment) => {
      segment.classList.toggle("is-selected", selected.has(segment.dataset.feature));
    });
    container.querySelectorAll("[data-label-for]").forEach((label) => {
      label.classList.toggle("is-selected", selected.has(label.dataset.labelFor));
    });
  }

  function normalizeMutationPositions(positions) {
    return Array.from(new Set(positions.map(Number))).map((position) => {
      if (!Number.isInteger(position) || position < 1 || position > GENOME_LENGTH) {
        throw new RangeError(
          `Mutation positions must be whole numbers from 1 to ${GENOME_LENGTH}.`
        );
      }
      return position;
    }).sort((a, b) => a - b);
  }

  function applyMutations(container, positions) {
    const svg = container.querySelector("svg");
    if (!svg) throw new Error("Render the visualizer before adding mutations.");
    const mutationPositions = normalizeMutationPositions(positions);
    svg.querySelector(".mitochondrial-genome__mutations")?.remove();
    const layer = svgElement("g", {
      class: "mitochondrial-genome__mutations",
      "aria-label": "Mutation positions"
    });

    mutationPositions.forEach((position) => {
      const affectedFeatures = FEATURES.filter((feature) =>
        featureContainsPosition(feature, position)
      );
      // Intergenic variants still receive a marker on the primary ring.
      const tracks = affectedFeatures.length
        ? Array.from(new Set(affectedFeatures.map((feature) => feature.track || "outer")))
        : ["outer"];

      tracks.forEach((track) => {
        const isInner = track === "inner";
        const innerRadius = isInner ? 360 : 396;
        const outerRadius = isInner ? 401 : 449;
        const angle = positionToAngle(position);
        const inner = point(500, 500, innerRadius, angle);
        const outer = point(500, 500, outerRadius, angle);
        const dot = point(500, 500, outerRadius + 8, angle);
        const line = svgElement("line", {
          class: "mitochondrial-genome__mutation-line",
          x1: inner.x,
          y1: inner.y,
          x2: outer.x,
          y2: outer.y,
          "data-mutation-position": position
        });
        const marker = svgElement("circle", {
          class: "mitochondrial-genome__mutation-dot",
          cx: dot.x,
          cy: dot.y,
          r: isInner ? 5 : 7,
          "data-mutation-position": position
        });
        const title = svgElement("title");
        const featureNames = affectedFeatures
          .filter((feature) => (feature.track || "outer") === track)
          .map((feature) => feature.name);
        title.textContent =
          `Position ${position.toLocaleString()}` +
          (featureNames.length ? ` · ${featureNames.join(", ")}` : " · intergenic");
        marker.appendChild(title);
        layer.append(line, marker);
      });
    });
    svg.appendChild(layer);
    return mutationPositions;
  }

  function genesForMutationPositions(positions) {
    const mutationPositions = normalizeMutationPositions(positions);
    const genes = new Set();
    mutationPositions.forEach((position) => {
      FEATURES.filter((feature) => featureContainsPosition(feature, position))
        .forEach((feature) => genes.add(feature.name));
    });
    return {
      positions: mutationPositions,
      genes: Array.from(genes)
    };
  }

  function applySelectedMutations(container, positions) {
    const result = genesForMutationPositions(positions);
    applyMutations(container, result.positions);
    applySelection(container, result.genes);
    return result;
  }

  function addDebugControls(
    container,
    selectableFeatures,
    selectedNames,
    mutationPositions,
    selectedMutationPositions
  ) {
    const panel = document.createElement("section");
    panel.className = "mitochondrial-genome__debug";
    panel.setAttribute("aria-label", "Visualizer debug controls");
    panel.innerHTML = `
      <div class="mitochondrial-genome__debug-header">
        <h2 class="mitochondrial-genome__debug-title">Debug: highlighted genes</h2>
        <div class="mitochondrial-genome__debug-actions">
          <button type="button" data-action="all">All</button>
          <button type="button" data-action="clear">Clear</button>
        </div>
      </div>
      <div class="mitochondrial-genome__debug-mutations">
        <label for="mitochondrial-mutation-positions">Mutation positions</label>
        <input
          id="mitochondrial-mutation-positions"
          type="text"
          inputmode="numeric"
          placeholder="e.g. 73, 9207, 10760"
          aria-describedby="mitochondrial-mutation-error"
        >
        <button type="button" data-action="mutations">Apply</button>
      </div>
      <div
        id="mitochondrial-mutation-error"
        class="mitochondrial-genome__debug-error"
        aria-live="polite"
      ></div>
      <div class="mitochondrial-genome__debug-mutations mitochondrial-genome__debug-mutations--selected">
        <label for="mitochondrial-selected-mutation-positions">Selected mutation positions</label>
        <input
          id="mitochondrial-selected-mutation-positions"
          type="text"
          inputmode="numeric"
          placeholder="e.g. 9207, 10760"
          aria-describedby="mitochondrial-selected-mutation-error"
        >
        <button type="button" data-action="selected-mutations">Apply selection</button>
      </div>
      <div
        id="mitochondrial-selected-mutation-error"
        class="mitochondrial-genome__debug-error"
        aria-live="polite"
      ></div>
      <div class="mitochondrial-genome__debug-options"></div>
    `;

    const mutationInput = panel.querySelector("#mitochondrial-mutation-positions");
    const mutationError = panel.querySelector("#mitochondrial-mutation-error");
    const selectedMutationInput = panel.querySelector(
      "#mitochondrial-selected-mutation-positions"
    );
    const selectedMutationError = panel.querySelector(
      "#mitochondrial-selected-mutation-error"
    );
    mutationInput.value = mutationPositions.join(", ");
    selectedMutationInput.value = selectedMutationPositions.join(", ");
    const options = panel.querySelector(".mitochondrial-genome__debug-options");
    selectableFeatures.forEach((feature) => {
      const label = document.createElement("label");
      label.className = "mitochondrial-genome__debug-option";
      const input = document.createElement("input");
      input.type = "checkbox";
      input.value = feature.name;
      input.checked = selectedNames.has(feature.name);
      label.append(input, feature.name.replace("MT-", ""));
      options.appendChild(label);
    });

    function updateFromInputs() {
      const checked = Array.from(
        options.querySelectorAll("input:checked"),
        (input) => input.value
      );
      applySelection(container, checked);
    }

    function syncGeneInputs(selectedGenes) {
      const genes = new Set(selectedGenes);
      options.querySelectorAll("input").forEach((input) => {
        input.checked = genes.has(input.value);
      });
    }

    options.addEventListener("change", updateFromInputs);
    panel.querySelector('[data-action="all"]').addEventListener("click", () => {
      options.querySelectorAll("input").forEach((input) => {
        input.checked = true;
      });
      updateFromInputs();
    });
    panel.querySelector('[data-action="clear"]').addEventListener("click", () => {
      options.querySelectorAll("input").forEach((input) => {
        input.checked = false;
      });
      updateFromInputs();
    });
    function updateMutations() {
      try {
        const values = mutationInput.value.trim()
          ? mutationInput.value.split(",").map((value) => value.trim())
          : [];
        const normalized = applyMutations(container, values);
        mutationInput.value = normalized.join(", ");
        mutationError.textContent = "";
      } catch (error) {
        mutationError.textContent = error.message;
      }
    }
    panel.querySelector('[data-action="mutations"]')
      .addEventListener("click", updateMutations);
    mutationInput.addEventListener("keydown", (event) => {
      if (event.key === "Enter") updateMutations();
    });
    function updateSelectedMutations() {
      try {
        const values = selectedMutationInput.value.trim()
          ? selectedMutationInput.value.split(",").map((value) => value.trim())
          : [];
        const result = applySelectedMutations(container, values);
        selectedMutationInput.value = result.positions.join(", ");
        mutationInput.value = result.positions.join(", ");
        syncGeneInputs(result.genes);
        selectedMutationError.textContent = "";
        mutationError.textContent = "";
      } catch (error) {
        selectedMutationError.textContent = error.message;
      }
    }
    panel.querySelector('[data-action="selected-mutations"]')
      .addEventListener("click", updateSelectedMutations);
    selectedMutationInput.addEventListener("keydown", (event) => {
      if (event.key === "Enter") updateSelectedMutations();
    });
    container.appendChild(panel);
  }

  function render(container, options) {
    if (!container) throw new Error("A container element is required.");
    const settings = Object.assign({
      showLabels: true,
      debugControls: false,
      selectedGenes: [],
      mutations: [],
      selectedMutations: []
    }, options);
    const selectedNames = new Set(settings.selectedGenes);
    const mutationPositions = normalizeMutationPositions(settings.mutations);
    const selectedMutationPositions = normalizeMutationPositions(
      settings.selectedMutations
    );
    const svg = svgElement("svg", {
      class: "mitochondrial-genome__canvas",
      viewBox: "0 0 1000 1000",
      role: "img",
      "aria-labelledby": "mitochondrial-genome-title mitochondrial-genome-description"
    });
    const title = svgElement("title", { id: "mitochondrial-genome-title" });
    title.textContent = "Human mitochondrial genome";
    const description = svgElement("desc", { id: "mitochondrial-genome-description" });
    description.textContent =
      "Circular map of the 16,569 base pair rCRS mitochondrial genome, with genes drawn to scale.";
    svg.append(title, description);

    FEATURES.filter((feature) => feature.track !== "inner")
      .forEach((feature) => addFeature(svg, feature, settings.showLabels));
    FEATURES.filter((feature) => feature.track === "inner")
      .forEach((feature) => addFeature(svg, feature, settings.showLabels));

    container.replaceChildren(svg);
    if (selectedMutationPositions.length) {
      const mutationSelection = genesForMutationPositions(selectedMutationPositions);
      mutationSelection.genes.forEach((gene) => selectedNames.add(gene));
      applyMutations(container, selectedMutationPositions);
    } else {
      applyMutations(container, mutationPositions);
    }
    applySelection(container, selectedNames);
    if (settings.debugControls) {
      addDebugControls(
        container,
        FEATURES.filter((feature) => feature.type !== "trna"),
        selectedNames,
        selectedMutationPositions.length
          ? selectedMutationPositions
          : mutationPositions,
        selectedMutationPositions
      );
    }
    return svg;
  }

  global.MitochondrialGenomeVisualizer = Object.freeze({
    GENOME_LENGTH,
    FEATURES,
    applyMutations,
    applySelectedMutations,
    applySelection,
    genesForMutationPositions,
    render
  });
}(window));
