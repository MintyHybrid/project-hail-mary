/* ===========================================================================
   network-nav.js — the landing reasoning graph and the Data & Code package
   graph. Both are the same force-directed simulation with different tuning.

   Differences from the design prototype, deliberately:
   - the simulation parks itself once kinetic energy drops below a threshold
     (the prototype ran requestAnimationFrame forever, which pins a core and
     drains laptop battery on a page people leave open while reading);
   - prefers-reduced-motion gets a settled static layout and no animation;
   - under 820px the graph is replaced by a plain ordered chapter list, since
     a 10-node force graph is unreadable and undraggable on a phone.
   =========================================================================== */

(function () {
  "use strict";

  var D = window.HailMary;
  if (!D) return;

  var REDUCED = window.matchMedia &&
    window.matchMedia("(prefers-reduced-motion: reduce)").matches;

  /* ------------------------------------------------------------ audience -- */
  // "expert" | "plain", shared with workbench.js through localStorage so the
  // choice survives navigation between the landing page and a chapter.
  var AUDIENCE_KEY = "hm-audience";

  function getAudience() {
    try {
      return localStorage.getItem(AUDIENCE_KEY) === "plain" ? "plain" : "expert";
    } catch (e) {
      return "expert";
    }
  }

  function setAudience(value) {
    try { localStorage.setItem(AUDIENCE_KEY, value); } catch (e) { /* private mode */ }
    document.dispatchEvent(new CustomEvent("hm:audience", { detail: value }));
  }

  window.HailMary.getAudience = getAudience;
  window.HailMary.setAudience = setAudience;

  /* ------------------------------------------------------ force simulation -- */

  function Sim(opts) {
    this.nodes = opts.nodes;          // [{id, x, y, vx, vy}]
    this.links = opts.links;          // [[idA, idB]]
    this.bounds = opts.bounds;        // fn -> {minX, maxX, minY, maxY}
    this.cfg = opts.cfg;
    this.onTick = opts.onTick;
    this.index = {};
    for (var i = 0; i < this.nodes.length; i++) this.index[this.nodes[i].id] = i;
    this.drag = null;
    this._raf = null;
    this._idle = 0;
  }

  Sim.prototype.step = function () {
    var ns = this.nodes, c = this.cfg, b = this.bounds();
    var f = new Array(ns.length);
    var i, j;
    for (i = 0; i < ns.length; i++) f[i] = { x: 0, y: 0 };

    // repulsion (all pairs — n is tiny, so O(n²) is fine)
    for (i = 0; i < ns.length; i++) {
      for (j = i + 1; j < ns.length; j++) {
        var dx = ns[i].x - ns[j].x, dy = ns[i].y - ns[j].y;
        var d2 = dx * dx + dy * dy || 1, d = Math.sqrt(d2);
        var rep = c.rep / d2, ux = dx / d, uy = dy / d;
        f[i].x += ux * rep; f[i].y += uy * rep;
        f[j].x -= ux * rep; f[j].y -= uy * rep;
      }
    }

    // spring attraction along links
    for (i = 0; i < this.links.length; i++) {
      var a = this.index[this.links[i][0]], z = this.index[this.links[i][1]];
      if (a == null || z == null) continue;
      var lx = ns[z].x - ns[a].x, ly = ns[z].y - ns[a].y;
      var ld = Math.hypot(lx, ly) || 1;
      var s = c.spr * (ld - c.rest), sx = lx / ld, sy = ly / ld;
      f[a].x += sx * s; f[a].y += sy * s;
      f[z].x -= sx * s; f[z].y -= sy * s;
    }

    var cx = (b.minX + b.maxX) / 2, cy = (b.minY + b.maxY) / 2;
    var energy = 0;

    for (i = 0; i < ns.length; i++) {
      var n = ns[i];
      if (this.drag && this.drag.moved && this.drag.id === n.id) {
        n.x = Math.max(b.minX, Math.min(b.maxX, this.drag.cx));
        n.y = Math.max(b.minY, Math.min(b.maxY, this.drag.cy));
        n.vx = n.vy = 0;
        energy += 100;              // dragging always counts as active
        continue;
      }
      f[i].x += (cx - n.x) * c.cen;
      f[i].y += (cy - n.y) * c.cen;
      n.vx = (n.vx + f[i].x) * c.damp;
      n.vy = (n.vy + f[i].y) * c.damp;
      n.x += n.vx; n.y += n.vy;
      if (n.x < b.minX) { n.x = b.minX; n.vx = -n.vx * 0.6; }
      if (n.x > b.maxX) { n.x = b.maxX; n.vx = -n.vx * 0.6; }
      if (n.y < b.minY) { n.y = b.minY; n.vy = -n.vy * 0.6; }
      if (n.y > b.maxY) { n.y = b.maxY; n.vy = -n.vy * 0.6; }
      energy += n.vx * n.vx + n.vy * n.vy;
    }
    return energy;
  };

  Sim.prototype.start = function () {
    var self = this;
    if (this._raf) return;
    var tick = function () {
      var energy = self.step();
      self.onTick();
      // park once everything has settled; any interaction calls start() again
      if (energy < 0.05) {
        self._idle++;
        if (self._idle > 30) { self._raf = null; return; }
      } else {
        self._idle = 0;
      }
      self._raf = requestAnimationFrame(tick);
    };
    this._raf = requestAnimationFrame(tick);
  };

  Sim.prototype.stop = function () {
    if (this._raf) cancelAnimationFrame(this._raf);
    this._raf = null;
  };

  Sim.prototype.settle = function (iterations) {
    for (var i = 0; i < iterations; i++) this.step();
    this.onTick();
  };

  /* --------------------------------------------------------------- utils -- */

  function el(tag, cls, parent) {
    var node = document.createElement(tag);
    if (cls) node.className = cls;
    if (parent) parent.appendChild(node);
    return node;
  }

  function svgEl(tag) {
    return document.createElementNS("http://www.w3.org/2000/svg", tag);
  }

  // deterministic pseudo-random so the initial layout is identical every load
  function rng(seed) {
    return function () {
      seed |= 0; seed = seed + 0x6D2B79F5 | 0;
      var t = Math.imul(seed ^ seed >>> 15, 1 | seed);
      t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t;
      return ((t ^ t >>> 14) >>> 0) / 4294967296;
    };
  }

  /* ================================================== LANDING REASONING MAP */

  function initLanding() {
    var hero = document.querySelector("[data-hm-hero]");
    if (!hero) return;

    document.documentElement.classList.add("hm-landing");

    var graph = hero.querySelector(".hm-graph");
    var svg = svgEl("svg");
    svg.setAttribute("width", "100%");
    svg.setAttribute("height", "100%");
    graph.appendChild(svg);

    // clicking empty canvas collapses any open card
    graph.addEventListener("pointerdown", function (e) {
      if (e.target === graph || e.target === svg) setExpanded(null);
    });

    var LEFT_GUTTER = 460;   // keeps nodes clear of the title overlay
    var TOP = 150, PAD = 70, BOTTOM = 110;

    var state = {
      expandedId: null,
      hoveredId: null,
      expandAll: false,
      audience: getAudience()
    };

    // --- build DOM ------------------------------------------------------
    var lines = [], views = {};

    D.EDGES.forEach(function () {
      var ln = svgEl("line");
      ln.setAttribute("stroke", "var(--hm-edge)");
      ln.setAttribute("stroke-width", "1.6");
      svg.appendChild(ln);
      lines.push(ln);
    });

    var simNodes = [];
    var rand = rng(7);

    D.CHAPTERS.forEach(function (meta) {
      var wrap = el("div", "hm-node", graph);
      if (meta.role === "origin") wrap.classList.add("is-origin");
      if (meta.role === "dest") wrap.classList.add("is-dest");

      // satellites
      var satCount = 3 + (meta.id.charCodeAt(0) % 3);
      var radius = meta.role === "origin" ? 31 : meta.role === "dest" ? 28 : 22;
      for (var s = 0; s < satCount; s++) {
        var sat = el("div", "sat", wrap);
        var ang = (Math.PI * 2 * s / satCount) - Math.PI / 2 +
                  (meta.id.charCodeAt(0) % 7) * 0.35;
        var rad = radius + 15 + (s % 2) * 7;
        sat.style.left = (Math.cos(ang) * rad) + "px";
        sat.style.top = (Math.sin(ang) * rad) + "px";
      }

      var circle = el("div", "circle", wrap);
      circle.textContent = meta.n;
      circle.setAttribute("role", "button");
      circle.setAttribute("tabindex", "0");
      circle.setAttribute("aria-label", "Chapter " + meta.n + ": " + meta.title);

      var label = el("div", "label", wrap);
      label.textContent = meta.short;
      label.style.top = (radius + 7) + "px";

      var summary = el("div", "summary", wrap);
      summary.style.top = (radius + 30) + "px";
      summary.hidden = true;

      var card = el("div", "card", wrap);
      card.hidden = true;
      card.style.transform = "translate(" + (radius + 16) + "px, -50%)";
      var kicker = el("div", "kicker", card);
      kicker.textContent = meta.id === "A" ? "APPENDIX" : "CHAPTER " + meta.n;
      var title = el("div", "title", card);
      title.textContent = meta.title;
      var body = el("div", "body", card);
      var enter = el("a", "enter", card);
      enter.href = meta.href;
      enter.textContent = (meta.id === "A"
        ? "explore the package network"
        : "enter chapter") + " →";

      views[meta.id] = {
        meta: meta, wrap: wrap, circle: circle, card: card,
        body: body, summary: summary, radius: radius
      };

      simNodes.push({
        id: meta.id,
        x: LEFT_GUTTER + rand() * 400,
        y: TOP + rand() * 380,
        vx: 0, vy: 0
      });

      // --- interaction ---
      circle.addEventListener("pointerenter", function () {
        state.hoveredId = meta.id; render();
      });
      circle.addEventListener("pointerleave", function () {
        if (state.hoveredId === meta.id) { state.hoveredId = null; render(); }
      });
      circle.addEventListener("keydown", function (e) {
        if (e.key === "Enter" || e.key === " ") {
          e.preventDefault();
          setExpanded(state.expandedId === meta.id ? null : meta.id);
        }
      });
      circle.addEventListener("pointerdown", function (e) {
        e.preventDefault();
        var r = graph.getBoundingClientRect();
        sim.drag = {
          id: meta.id, moved: false,
          sx: e.clientX - r.left, sy: e.clientY - r.top,
          cx: e.clientX - r.left, cy: e.clientY - r.top
        };
        sim.start();
      });
    });

    function setExpanded(id) {
      state.expandedId = id;
      render();
    }

    var sim = new Sim({
      nodes: simNodes,
      links: D.EDGES,
      cfg: { rep: 9500, spr: 0.016, rest: 155, cen: 0.0018, damp: 0.94 },
      bounds: function () {
        var r = graph.getBoundingClientRect();
        var W = Math.max(720, r.width), H = Math.max(560, r.height);
        return {
          minX: LEFT_GUTTER, maxX: W - PAD,
          // BOTTOM (not PAD) keeps node labels and summaries, which hang up to
          // ~45px below a node's centre, clear of the scroll cue at the foot
          // of the hero.
          minY: TOP, maxY: Math.min(H - BOTTOM, TOP + 470)
        };
      },
      onTick: paint
    });

    // --- painting -------------------------------------------------------
    function paint() {
      var byId = {};
      simNodes.forEach(function (n) { byId[n.id] = n; });

      D.EDGES.forEach(function (e, i) {
        var a = byId[e[0]], b = byId[e[1]];
        if (!a || !b) return;
        var hot = state.expandedId === e[0] || state.expandedId === e[1];
        var ln = lines[i];
        ln.setAttribute("x1", a.x); ln.setAttribute("y1", a.y);
        ln.setAttribute("x2", b.x); ln.setAttribute("y2", b.y);
        ln.setAttribute("stroke", hot ? "var(--hm-gold)" : "var(--hm-edge)");
        ln.setAttribute("stroke-width", hot ? "2.4" : "1.6");
      });

      simNodes.forEach(function (n) {
        var v = views[n.id];
        v.wrap.style.transform = "translate(" + n.x + "px," + n.y + "px)";
      });
    }

    function render() {
      var graphW = graph.getBoundingClientRect().width;

      Object.keys(views).forEach(function (id) {
        var v = views[id], meta = v.meta;
        var expanded = state.expandedId === id && !state.expandAll;
        var text = state.audience === "expert" ? meta.expert : meta.plain;

        v.wrap.classList.toggle("is-hovered", state.hoveredId === id);
        v.wrap.classList.toggle("is-expanded", state.expandedId === id);

        v.summary.hidden = !state.expandAll;
        if (state.expandAll) v.summary.textContent = text;

        v.card.hidden = !expanded;
        if (expanded) {
          v.body.textContent = text;
          // flip the card to the left when it would overflow the viewport
          var n = simNodes.filter(function (s) { return s.id === id; })[0];
          var flip = n && (n.x > graphW - 300) && (n.x - v.radius - 276 > LEFT_GUTTER + 8);
          v.card.style.transform = flip
            ? "translate(-" + (276 + v.radius) + "px, -50%)"
            : "translate(" + (v.radius + 16) + "px, -50%)";
        }
        v.circle.setAttribute("aria-expanded", String(state.expandedId === id));
      });

      paint();
    }

    // --- global pointer handling (drag + click discrimination) ----------
    document.addEventListener("pointermove", function (e) {
      if (!sim.drag) return;
      var r = graph.getBoundingClientRect();
      var x = e.clientX - r.left, y = e.clientY - r.top;
      if (Math.hypot(x - sim.drag.sx, y - sim.drag.sy) > 4) sim.drag.moved = true;
      sim.drag.cx = x; sim.drag.cy = y;
    });

    document.addEventListener("pointerup", function () {
      if (sim.drag && !sim.drag.moved) {
        setExpanded(state.expandedId === sim.drag.id ? null : sim.drag.id);
      }
      sim.drag = null;
    });

    // --- header controls -------------------------------------------------
    var segButtons = hero.querySelectorAll("[data-hm-audience]");
    Array.prototype.forEach.call(segButtons, function (btn) {
      btn.addEventListener("click", function () {
        setAudience(btn.getAttribute("data-hm-audience"));
      });
    });

    document.addEventListener("hm:audience", function (e) {
      state.audience = e.detail;
      Array.prototype.forEach.call(segButtons, function (btn) {
        btn.setAttribute("aria-pressed",
          String(btn.getAttribute("data-hm-audience") === e.detail));
      });
      render();
    });

    var toggle = hero.querySelector("[data-hm-expand-all]");
    if (toggle) {
      toggle.addEventListener("click", function () {
        state.expandAll = !state.expandAll;
        state.expandedId = null;
        toggle.setAttribute("aria-pressed", String(state.expandAll));
        render();
      });
    }

    // reflect the stored audience into the buttons on first paint
    Array.prototype.forEach.call(segButtons, function (btn) {
      btn.setAttribute("aria-pressed",
        String(btn.getAttribute("data-hm-audience") === state.audience));
    });

    // --- resize ----------------------------------------------------------
    var rt;
    window.addEventListener("resize", function () {
      clearTimeout(rt);
      rt = setTimeout(function () { sim.start(); render(); }, 160);
    });

    // --- go --------------------------------------------------------------
    sim.settle(REDUCED ? 400 : 60);
    render();
    if (!REDUCED) sim.start();

    buildFallbackList(hero, state);
  }

  /* Plain list for narrow screens and for anyone the graph doesn't serve. */
  function buildFallbackList(hero, state) {
    var host = hero.querySelector(".hm-fallback-list");
    if (!host) return;
    var ol = el("ol", null, host);
    D.CHAPTERS.forEach(function (m) {
      var li = el("li", null, ol);
      var a = el("a", null, li);
      a.href = m.href;
      var n = el("span", "n", a); n.textContent = m.n;
      var t = el("span", "t", a); t.textContent = m.title;
      var s = el("span", "s", a);
      s.textContent = state.audience === "expert" ? m.expert : m.plain;
      document.addEventListener("hm:audience", function (e) {
        s.textContent = e.detail === "expert" ? m.expert : m.plain;
      });
    });
  }

  /* ==================================================== PACKAGE NETWORK ==== */

  function initPackages() {
    var host = document.querySelector("[data-hm-packages]");
    if (!host) return;

    var svg = svgEl("svg");
    svg.setAttribute("width", "100%");
    svg.setAttribute("height", "100%");
    host.appendChild(svg);

    var tooltip = el("div", "hm-pkg-tooltip", host);
    tooltip.hidden = true;
    var tMeta = el("div", "meta", tooltip);
    var tName = el("div", "name", tooltip);
    var tDesc = el("div", "desc", tooltip);

    // links: package -> each chapter that uses it
    var links = [];
    D.PACKAGES.forEach(function (pk) {
      pk.chapters.forEach(function (cid) { links.push([pk.id, cid]); });
    });

    var rect = host.getBoundingClientRect();
    var W = Math.max(600, rect.width), H = Math.max(420, rect.height);
    var cx = W / 2, cy = H / 2;
    var R1 = Math.min(W, H) * 0.24, R2 = Math.min(W, H) * 0.46;

    var simNodes = [], views = {};
    var hoveredId = null;

    D.CHAPTERS.forEach(function (m, i) {
      var ang = (Math.PI * 2 * i / D.CHAPTERS.length) - Math.PI / 2;
      simNodes.push({
        id: m.id, kind: "chapter",
        x: cx + Math.cos(ang) * R1, y: cy + Math.sin(ang) * R1, vx: 0, vy: 0
      });
    });

    var sorted = D.PACKAGES.slice().sort(function (a, b) {
      return a.topic.localeCompare(b.topic) || a.id.localeCompare(b.id);
    });
    sorted.forEach(function (pk, i) {
      var ang = (Math.PI * 2 * i / sorted.length) - Math.PI / 2;
      simNodes.push({
        id: pk.id, kind: "pkg",
        x: cx + Math.cos(ang) * R2, y: cy + Math.sin(ang) * R2, vx: 0, vy: 0
      });
    });

    var lines = links.map(function () {
      var ln = svgEl("line");
      svg.appendChild(ln);
      return ln;
    });

    simNodes.forEach(function (n) {
      var isChapter = n.kind === "chapter";
      var meta = isChapter
        ? D.CHAPTERS.filter(function (m) { return m.id === n.id; })[0]
        : D.PACKAGES.filter(function (p) { return p.id === n.id; })[0];

      var wrap = el("div", "hm-pkg-node", host);
      wrap.classList.add(isChapter ? "kind-chapter" : "kind-pkg");
      if (!isChapter) {
        wrap.classList.add(meta.lang === "Python" ? "lang-python" : "lang-r");
      }
      var chip = el("div", "chip", wrap);
      chip.textContent = isChapter ? meta.short : meta.id;

      if (!isChapter) {
        chip.addEventListener("pointerenter", function () {
          hoveredId = n.id; render();
        });
        chip.addEventListener("pointerleave", function () {
          if (hoveredId === n.id) { hoveredId = null; render(); }
        });
      }
      chip.addEventListener("pointerdown", function (e) {
        e.preventDefault();
        var r = host.getBoundingClientRect();
        sim.drag = {
          id: n.id, moved: false,
          sx: e.clientX - r.left, sy: e.clientY - r.top,
          cx: e.clientX - r.left, cy: e.clientY - r.top
        };
        sim.start();
      });

      views[n.id] = { wrap: wrap, chip: chip, meta: meta, kind: n.kind };
    });

    var sim = new Sim({
      nodes: simNodes,
      links: links,
      cfg: { rep: 7000, spr: 0.02, rest: 118, cen: 0.006, damp: 0.9 },
      bounds: function () {
        var r = host.getBoundingClientRect();
        return {
          minX: 32, maxX: Math.max(600, r.width) - 32,
          minY: 32, maxY: Math.max(420, r.height) - 56
        };
      },
      onTick: paint
    });

    function paint() {
      var byId = {};
      simNodes.forEach(function (n) { byId[n.id] = n; });

      links.forEach(function (l, i) {
        var p = byId[l[0]], c = byId[l[1]];
        if (!p || !c) return;
        var pk = D.PACKAGES.filter(function (x) { return x.id === l[0]; })[0];
        var hot = hoveredId === l[0];
        var col = pk && pk.lang === "Python"
          ? (hot ? "var(--hm-gold)" : "var(--hm-gold-dim)")
          : (hot ? "var(--hm-spark)" : "var(--hm-spark-dim)");
        var ln = lines[i];
        ln.setAttribute("x1", p.x); ln.setAttribute("y1", p.y);
        ln.setAttribute("x2", c.x); ln.setAttribute("y2", c.y);
        ln.setAttribute("stroke", col);
        ln.setAttribute("stroke-width", hot ? "2.2" : "1.1");
        ln.setAttribute("opacity", hot ? "1" : "0.65");
      });

      simNodes.forEach(function (n) {
        views[n.id].wrap.style.transform =
          "translate(" + n.x + "px," + n.y + "px)";
      });
    }

    function render() {
      Object.keys(views).forEach(function (id) {
        views[id].wrap.classList.toggle("is-hovered", hoveredId === id);
      });
      if (hoveredId) {
        var pk = D.PACKAGES.filter(function (p) { return p.id === hoveredId; })[0];
        var pos = simNodes.filter(function (n) { return n.id === hoveredId; })[0];
        if (pk && pos) {
          var r = host.getBoundingClientRect();
          tMeta.textContent = pk.lang + " · " + pk.topic;
          tMeta.style.color = pk.lang === "Python"
            ? "var(--hm-gold)" : "var(--hm-spark)";
          tName.textContent = pk.id;
          tDesc.textContent = pk.desc;
          tooltip.hidden = false;
          tooltip.style.left = Math.min(pos.x + 20, r.width - 240) + "px";
          tooltip.style.top = Math.max(pos.y - 70, 10) + "px";
        }
      } else {
        tooltip.hidden = true;
      }
      paint();
    }

    document.addEventListener("pointermove", function (e) {
      if (!sim.drag) return;
      var r = host.getBoundingClientRect();
      var x = e.clientX - r.left, y = e.clientY - r.top;
      if (Math.hypot(x - sim.drag.sx, y - sim.drag.sy) > 4) sim.drag.moved = true;
      sim.drag.cx = x; sim.drag.cy = y;
    });
    document.addEventListener("pointerup", function () { sim.drag = null; });

    var rt;
    window.addEventListener("resize", function () {
      clearTimeout(rt);
      rt = setTimeout(function () { sim.start(); }, 160);
    });

    sim.settle(REDUCED ? 400 : 80);
    render();
    if (!REDUCED) sim.start();
  }

  /* ------------------------------------------------------------------ go -- */
  function boot() {
    initLanding();
    initPackages();
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", boot);
  } else {
    boot();
  }
})();
