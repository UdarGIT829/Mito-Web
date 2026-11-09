// Faint ball-and-stick "molecule" network behind the page.
//
// Atoms sit on a jittered grid (a few large, most small); each one bonds to
// one or two of its nearest neighbours, so the bonds form loose chains
// rather than a mesh. Drawn once (per load and resize) onto a canvas taller
// than the window, which drifts up a little as the page scrolls.
(function () {
  var canvas = document.getElementById("molecules");
  if (!canvas || !canvas.getContext) return;
  var ctx = canvas.getContext("2d");
  var PARALLAX = 0.35;         // canvas extra height, as a fraction of the window
  var reduceMotion = window.matchMedia("(prefers-reduced-motion: reduce)").matches;

  function draw() {
    var dpr = Math.min(window.devicePixelRatio || 1, 2);
    var w = window.innerWidth, h = Math.round(window.innerHeight * (1 + PARALLAX));
    canvas.width = Math.round(w * dpr); canvas.height = Math.round(h * dpr);
    canvas.style.width = w + "px"; canvas.style.height = h + "px";
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, w, h);

    var cell = Math.max(150, Math.min(230, w / 7));
    var atoms = [];
    for (var y = -cell / 2; y < h + cell; y += cell) {
      for (var x = -cell / 2; x < w + cell; x += cell) {
        if (Math.random() < 0.18) continue;            // leave some holes
        var big = Math.random() < 0.12;
        atoms.push({
          x: x + (Math.random() - 0.5) * cell * 0.9,
          y: y + (Math.random() - 0.5) * cell * 0.9,
          r: big ? 22 + Math.random() * 26 : 5 + Math.random() * 6,
        });
      }
    }
    // Bonds: each atom to its 1-2 nearest neighbours within reach.
    var bonds = [], seen = {};
    atoms.forEach(function (a, i) {
      var near = atoms.map(function (b, j) { return { j: j, d: Math.hypot(a.x - b.x, a.y - b.y) }; })
        .filter(function (n) { return n.j !== i && n.d < cell * 1.6; })
        .sort(function (p, q) { return p.d - q.d; })
        .slice(0, Math.random() < 0.5 ? 1 : 2);
      near.forEach(function (n) {
        var k = Math.min(i, n.j) + ":" + Math.max(i, n.j);
        if (!seen[k]) { seen[k] = 1; bonds.push([a, atoms[n.j]]); }
      });
    });

    ctx.lineCap = "round";
    bonds.forEach(function (b) {
      ctx.strokeStyle = "rgba(255,255,255,0.55)";
      ctx.lineWidth = 3 + Math.random() * 4;
      ctx.beginPath(); ctx.moveTo(b[0].x, b[0].y); ctx.lineTo(b[1].x, b[1].y); ctx.stroke();
    });
    atoms.forEach(function (a) {
      ctx.fillStyle = a.r > 20 ? "rgba(255,255,255,0.42)" : "rgba(255,255,255,0.62)";
      ctx.beginPath(); ctx.arc(a.x, a.y, a.r, 0, Math.PI * 2); ctx.fill();
    });
  }

  function scroll() {
    if (reduceMotion) return;
    var max = document.documentElement.scrollHeight - window.innerHeight;
    var p = max > 0 ? window.scrollY / max : 0;
    canvas.style.transform = "translateY(" + (-p * window.innerHeight * PARALLAX).toFixed(1) + "px)";
  }

  var resizeTimer;
  window.addEventListener("resize", function () { clearTimeout(resizeTimer); resizeTimer = setTimeout(function () { draw(); scroll(); }, 150); });
  window.addEventListener("scroll", scroll, { passive: true });
  draw(); scroll();
})();
