// Live 3D DNA double helix on the right of the project page.
//
// - Two navy backbones on opposite sides (half a turn apart), so the
//   strands cross evenly. A full turn is about three times as tall as the
//   helix is wide.
// - BASES_PER_TURN rungs per turn, each a diameter split into two rounded
//   halves that meet at the axis, with a white knob at each inner end and a
//   grey outline (drawn on the silhouette rim), like a flat vector
//   illustration.
// - Each half shades from white (facing the viewer) to grey (facing away),
//   recomputed every frame as the helix turns.
// - Scrolling turns the helix and slides it upward like a screw. The helix
//   repeats every PITCH units, so the slide wraps without a visible jump.
//   It also turns slowly on its own unless the visitor prefers reduced
//   motion. #helix is removed if WebGL is unavailable.

import * as THREE from "./vendor/three.module.min.js";

const RADIUS = 1.0;
const PITCH = 6.0;              // rise per full turn (3 x the diameter)
const BASES_PER_TURN = 10;
const TURNS = 4;                // covers the view plus one wrap
const BACKBONE_R = 0.17;
const RUNG_R = 0.12;
const KNOB_R = 0.058;
const RIM = 0.42;               // silhouette band drawn as outline (|n·v| below this)
const JOINT_GAP = 0.05;         // gap between the two halves at the axis
const IDLE_SPIN = 0.07;         // rad/s
const SCROLL_SPIN = 0.0018;     // rad per scrolled pixel
const SCROLL_RISE = 0.0035;     // units per scrolled pixel
const DEPTH_PUSH = 0.24;        // rung depth bias (> BACKBONE_R), so strands cover rung ends cleanly

const FRONT = new THREE.Color("#f3f5f8");
const BACK = new THREE.Color("#c3c9d1");
const OUTLINE_COLOR = new THREE.Color("#a9b0b9");
const BACKBONE = new THREE.Color("#12213a");

const reduceMotion = window.matchMedia("(prefers-reduced-motion: reduce)").matches;
const container = document.getElementById("helix");

let renderer;
try {
  renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true, powerPreference: "low-power" });
} catch (e) {
  container.remove();
  throw e;
}
renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 1.75));
renderer.setClearColor(0x000000, 0);
container.appendChild(renderer.domElement);

// Two-step ramp: flat shading with one soft shadow band.
const ramp = new THREE.DataTexture(new Uint8Array([170, 170, 170, 255, 255, 255, 255, 255]), 2, 1);
ramp.minFilter = ramp.magFilter = THREE.NearestFilter;
ramp.needsUpdate = true;

const scene = new THREE.Scene();
scene.add(new THREE.AmbientLight(0xffffff, 1.35));
const key = new THREE.DirectionalLight(0xffffff, 1.3);
key.position.set(-2, 5, 6);
scene.add(key);

const helix = new THREE.Group();          // spins and slides
const frame = new THREE.Group();          // placement on screen
frame.add(helix);
scene.add(frame);

// ── Backbones ─────────────────────────────────────────────────────────────
const height = TURNS * PITCH;
const y0 = -height / 2;
const point = (angle, y, r = RADIUS) => new THREE.Vector3(r * Math.cos(angle), y, r * Math.sin(angle));
const backboneMat = new THREE.MeshToonMaterial({ color: BACKBONE, gradientMap: ramp });
for (const phase of [0, Math.PI]) {
  const pts = [];
  const steps = TURNS * 64;
  for (let i = 0; i <= steps; i++) {
    const t = i / steps;
    pts.push(point(phase + t * TURNS * Math.PI * 2, y0 + t * height));
  }
  const tube = new THREE.TubeGeometry(new THREE.CatmullRomCurve3(pts), steps * 2, BACKBONE_R, 16, false);
  helix.add(new THREE.Mesh(tube, backboneMat));
}

// ── Rungs ─────────────────────────────────────────────────────────────────
// Where a rung enters its strand, the two shapes would intersect along a
// ragged curve. Instead every rung part is drawn as if DEPTH_PUSH further
// from the camera (depth only; its position on screen is unchanged), so near
// a strand the strand always wins and covers the rung end with a clean edge,
// like a flat illustration. Elsewhere the depth order is unaffected.
// The grey outline is drawn in the same shader, by colouring the rim of
// the silhouette (where the surface turns away from the viewer), so it is
// part of the rung and can't stick out past a strand the way a separate,
// larger outline shell would.
// The shader's output is already sRGB at this point, so use sRGB components.
const rim = OUTLINE_COLOR.clone().convertLinearToSRGB();
const outlineGLSL = `vec3(${rim.r.toFixed(4)}, ${rim.g.toFixed(4)}, ${rim.b.toFixed(4)})`;
function pushDepth(material, rimOutline = false) {
  material.onBeforeCompile = (shader) => {
    shader.vertexShader = shader.vertexShader.replace("#include <project_vertex>", `#include <project_vertex>
      {
        vec4 pushed = projectionMatrix * (mvPosition + vec4(0.0, 0.0, -${DEPTH_PUSH.toFixed(3)}, 0.0));
        gl_Position.z = pushed.z / pushed.w * gl_Position.w;
      }`);
    if (rimOutline) {
      shader.fragmentShader = shader.fragmentShader.replace("#include <dithering_fragment>", `
        float facing = abs(dot(normalize(normal), normalize(vViewPosition)));
        gl_FragColor.rgb = mix(gl_FragColor.rgb, ${outlineGLSL}, step(facing, ${RIM.toFixed(3)}));
        #include <dithering_fragment>`);
    }
  };
  material.customProgramCacheKey = () => (rimOutline ? "rung-rim" : "rung");
  return material;
}

// One half: a capsule from just short of the axis to the backbone's centre
// line, so its rounded outer end (RUNG_R < BACKBONE_R)
// is buried inside the strand and never pokes through it.
const halfLen = RADIUS - JOINT_GAP / 2 - RUNG_R;          // between the two cap centres
const capsule = new THREE.CapsuleGeometry(RUNG_R, halfLen, 6, 14);
const knob = new THREE.SphereGeometry(KNOB_R, 12, 10);
const knobMat = pushDepth(new THREE.MeshBasicMaterial({ color: 0xffffff }));
const up = new THREE.Vector3(0, 1, 0);
const halves = [];   // { mesh, material } for per-frame shading
const total = TURNS * BASES_PER_TURN;
for (let i = 0; i < total; i++) {
  const a = (i / BASES_PER_TURN) * Math.PI * 2;
  const y = y0 + (i / BASES_PER_TURN) * PITCH;            // exactly where both strands pass at angle a
  for (const side of [0, Math.PI]) {
    const outward = point(a + side, 0, 1);                 // unit direction from the axis
    const center = outward.clone().multiplyScalar(JOINT_GAP / 2 + RUNG_R + halfLen / 2).setY(y);
    const q = new THREE.Quaternion().setFromUnitVectors(up, outward);
    const material = pushDepth(new THREE.MeshToonMaterial({ color: FRONT.clone(), gradientMap: ramp }), true);
    const mesh = new THREE.Mesh(capsule, material);
    mesh.position.copy(center); mesh.quaternion.copy(q);
    // The knob sits at the inner end and is nudged toward the viewer every
    // frame, so it reads as a bright dot on the end of the bar.
    const dot = new THREE.Mesh(knob, knobMat);
    const dotBase = outward.clone().multiplyScalar(JOINT_GAP / 2 + RUNG_R * 1.2).setY(y);
    helix.add(mesh, dot);
    halves.push({ mesh, material, outward, dot, dotBase });
  }
}

// ── Camera, layout, scroll ────────────────────────────────────────────────
const camera = new THREE.PerspectiveCamera(30, 1, 0.1, 100);
camera.position.set(0, 0, 12);

function layout() {
  const w = window.innerWidth, h = window.innerHeight;
  renderer.setSize(w, h, false);
  camera.aspect = w / h;
  camera.updateProjectionMatrix();
  const halfW = Math.tan(THREE.MathUtils.degToRad(15)) * 12 * camera.aspect;
  // Wide screens: beside the content column. Narrow: at the edge, and CSS
  // fades it behind the text.
  frame.position.x = w >= 1100 ? halfW * 0.6 : halfW * 0.8;
}
window.addEventListener("resize", layout);
layout();
frame.rotation.z = -0.04;

let scrollTarget = window.scrollY, scrollNow = scrollTarget;
window.addEventListener("scroll", () => { scrollTarget = window.scrollY; }, { passive: true });

const tmp = new THREE.Vector3(), toViewer = new THREE.Vector3();
let idle = 0, last = performance.now();
function tick(now) {
  const dt = Math.min((now - last) / 1000, 0.1);
  last = now;
  if (!reduceMotion) idle += dt * IDLE_SPIN;
  scrollNow += (scrollTarget - scrollNow) * Math.min(1, dt * 5);   // ease toward the scroll position

  helix.rotation.y = idle + scrollNow * SCROLL_SPIN;
  helix.position.y = (scrollNow * SCROLL_RISE) % PITCH;            // periodic, so the wrap is invisible

  // Shade each half by how much it faces the viewer (+z after rotation).
  toViewer.set(0, 0, RUNG_R * 0.72).applyAxisAngle(up, -helix.rotation.y);   // +z in helix space
  for (const h of halves) {
    tmp.copy(h.outward).applyAxisAngle(up, helix.rotation.y);
    h.material.color.copy(BACK).lerp(FRONT, (tmp.z + 1) / 2);
    h.dot.position.copy(h.dotBase).add(toViewer);
  }
  renderer.render(scene, camera);
  requestAnimationFrame(tick);
}
requestAnimationFrame(tick);
container.classList.add("live");
