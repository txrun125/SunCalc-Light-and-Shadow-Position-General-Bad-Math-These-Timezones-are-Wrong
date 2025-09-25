/* ===== vector helpers ===== */
const deg = (r) => (r * 180) / Math.PI,
	rad = (d) => (d * Math.PI) / 180;
const vadd = (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
const vsub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
const vmul = (a, s) => [a[0] * s, a[1] * s, a[2] * s];
const vdot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const vcross = (a, b) => [
	a[1] * b[2] - a[2] * b[1],
	a[2] * b[0] - a[0] * b[2],
	a[0] * b[1] - a[1] * b[0]
];
const vlen = (a) => Math.hypot(a[0], a[1], a[2]);
const vnorm = (a) => {
	const L = vlen(a) || 1;
	return [a[0] / L, a[1] / L, a[2] / L];
};

function normalToEuler([x, y, z]) {
	const yaw = Math.atan2(x, z);
	const hyp = Math.hypot(x, z);
	const pitch = Math.atan2(y, hyp);
	return [pitch, yaw];
}
function tangentFrame(n) {
	// stable ONB on face plane
	const ref = Math.abs(n[1]) < 0.9 ? [0, 1, 0] : [1, 0, 0];
	const u = vnorm(vcross(ref, n));
	const v = vnorm(vcross(n, u));
	return { u, v, n };
}

/* ===== base icosahedron ===== */
function icosahedron() {
	const φ = (1 + Math.sqrt(5)) / 2;
	const V = [
		[-1, φ, 0],
		[1, φ, 0],
		[-1, -φ, 0],
		[1, -φ, 0],
		[0, -1, φ],
		[0, 1, φ],
		[0, -1, -φ],
		[0, 1, -φ],
		[φ, 0, -1],
		[φ, 0, 1],
		[-φ, 0, -1],
		[-φ, 0, 1]
	].map(vnorm);
	const F = [
		[0, 11, 5],
		[0, 5, 1],
		[0, 1, 7],
		[0, 7, 10],
		[0, 10, 11],
		[1, 5, 9],
		[5, 11, 4],
		[11, 10, 2],
		[10, 7, 6],
		[7, 1, 8],
		[3, 9, 4],
		[3, 4, 2],
		[3, 2, 6],
		[3, 6, 8],
		[3, 8, 9],
		[4, 9, 5],
		[2, 4, 11],
		[6, 2, 10],
		[8, 6, 7],
		[9, 8, 1]
	];
	return { V, F };
}

/* build unique edges + vertex neighbors (for pent ordering) */
function edgesAndNeighbors(F, nVerts) {
	const Eset = new Set();
	const nbr = Array.from({ length: nVerts }, () => new Set());
	for (const [a, b, c] of F) {
		const add = (i, j) => {
			const k = i < j ? `${i},${j}` : `${j},${i}`;
			Eset.add(k);
			nbr[i].add(j);
			nbr[j].add(i);
		};
		add(a, b);
		add(b, c);
		add(c, a);
	}
	const E = Array.from(Eset).map((k) => k.split(",").map((x) => +x));
	const neighbors = nbr.map((s) => Array.from(s));
	return { E, neighbors };
}

/* sort 5 neighbors around vertex i in a stable CCW order */
function orderNeighborsCCW(i, neighbors, V) {
	const n = V[i];
	const { u, v } = tangentFrame(n);
	const ang = (j) => {
		const w = vsub(V[j], n); // direction from vertex to neighbor
		const proj = vsub(w, vmul(n, vdot(w, n)));
		const x = vdot(proj, u),
			y = vdot(proj, v);
		return Math.atan2(y, x);
	};
	return neighbors[i].slice().sort((a, b) => ang(a) - ang(b));
}

/* truncation point P(i->j) at fraction s from i toward j (no normalization) */
const truncPoint = (Vi, Vj, s) => vadd(vmul(Vi, 1 - s), vmul(Vj, s));

/* === build truncated icosahedron faces (pent + hex) ===
   - Pentagon at each icosa vertex i: the 5 points P(i->j_k), neighbors ordered CCW.
   - Hexagon at each icosa triangle [a,b,c]: [P(a->b), P(b->a), P(b->c), P(c->b), P(c->a), P(a->c)].
*/
function buildTruncatedIcosahedron(s) {
	const { V, F } = icosahedron();
	const { neighbors } = edgesAndNeighbors(F, V.length);

	const faces = [];

	// pentagons (12)
	for (let i = 0; i < V.length; i++) {
		const ring = orderNeighborsCCW(i, neighbors, V);
		const verts = ring.map((j) => truncPoint(V[i], V[j], s));
		faces.push({ kind: "pent", verts });
	}

	// hexagons (20)
	for (const [a, b, c] of F) {
		const verts = [
			truncPoint(V[a], V[b], s),
			truncPoint(V[b], V[a], s),
			truncPoint(V[b], V[c], s),
			truncPoint(V[c], V[b], s),
			truncPoint(V[c], V[a], s),
			truncPoint(V[a], V[c], s)
		];
		faces.push({ kind: "hex", verts });
	}

	// ensure outward normals + plane offsets
	for (const f of faces) {
		// compute normal from first three vertices
		let n = vnorm(
			vcross(vsub(f.verts[1], f.verts[0]), vsub(f.verts[2], f.verts[0]))
		);
		// make it outward (pointing roughly like average vertex position)
		const avg = f.verts
			.reduce((a, b) => vadd(a, b), [0, 0, 0])
			.map((x) => x / f.verts.length);
		if (vdot(n, avg) < 0) {
			n = vmul(n, -1);
			f.verts.reverse();
		}
		f.n = n;
		f.d = vdot(n, f.verts[0]); // plane offset from origin
	}
	return faces;
}

/* project each face’s 3D vertices to its local (u,v) plane coords, centered */
function face2D(f) {
	const { u, v, n } = tangentFrame(f.n);
	const p0 = vmul(n, f.d); // plane origin (closest point to origin)
	let pts = f.verts.map((p) => {
		const w = vsub(p, p0);
		return [vdot(w, u), vdot(w, v)];
	});
	// center around centroid
	const c = pts
		.reduce((a, b) => [a[0] + b[0], a[1] + b[1]], [0, 0])
		.map((x) => x / pts.length);
	pts = pts.map(([x, y]) => [x - c[0], y - c[1]]);
	// bounding box
	let minX = Infinity,
		maxX = -Infinity,
		minY = Infinity,
		maxY = -Infinity;
	for (const [x, y] of pts) {
		if (x < minX) minX = x;
		if (x > maxX) maxX = x;
		if (y < minY) minY = y;
		if (y > maxY) maxY = y;
	}
	const w = maxX - minX,
		h = maxY - minY;
	// remap to 0..1 box
	const pts01 = pts.map(([x, y]) => [(x - minX) / w, (y - minY) / h]);
	return { pts01, w, h };
}

/* create DOM for one face using exact projected polygon */
function addFace(ball, f, scale) {
	const { pts01, w, h } = face2D(f);
	const sizeW = w * scale,
		sizeH = h * scale;

	const [pitch, yaw] = normalToEuler(f.n);
	const wrap = document.createElement("div");
	wrap.className = "face";
	wrap.style.width = `${sizeW}px`;
	wrap.style.height = `${sizeH}px`;
	// Put plane at its true offset d (no shared radius!)
	wrap.style.transform = `rotateY(${deg(yaw)}deg) rotateX(${deg(
		-pitch
	)}deg) translateZ(${f.d * scale}px)`;

	const tile = document.createElement("div");
	tile.className = `tile ${f.kind}`;
	if (f.kind === "hex" && f.n[0] === 0) {
		console.log(tile, f);
		wrap.style.transform = `rotateY(${deg(yaw)}deg) rotateX(${deg(
			-pitch
		)}deg) rotateZ(${f.n[2] > 0 ? 30 : -30}deg) 
  translateZ(${f.d * scale}px)`;
	}

	// build clip-path from 0..1 points
	const poly = pts01
		.map(([x, y]) => `${(x * 100).toFixed(3)}% ${(y * 100).toFixed(3)}%`)
		.join(",");
	tile.style.clipPath = `polygon(${poly})`;

	wrap.appendChild(tile);
	ball.appendChild(wrap);
}

/* ===== render ===== */
const ball = document.getElementById("ball");
const sEl = document.getElementById("s");
const scaleEl = document.getElementById("scale");
const spinEl = document.getElementById("spin");

function render() {
	const s = parseFloat(sEl.value); // truncation fraction (≈0.333 looks “soccer”)
	const scale = parseFloat(scaleEl.value);
	document.documentElement.style.setProperty("--scale", scale + "px");
	ball.classList.add("spin");

	ball.innerHTML = "";
	const faces = buildTruncatedIcosahedron(s);
	// draw hexes first, then pents (purely visual layering)
	faces.sort((a, b) => (a.kind === b.kind ? 0 : a.kind === "hex" ? -1 : 1));
	for (const f of faces) addFace(ball, f, scale);
}

sEl.addEventListener("input", render);
scaleEl.addEventListener("input", render);

render();
