/* ===========================================================================
   workbench.js — the per-chapter reasoning stepper, download bar and audience
   toggle.

   These are injected client-side rather than written into the .qmd sources on
   purpose: the chapters are `execute: freeze: auto`, and editing them would
   invalidate _freeze/ and force a full re-render — which needs mafft, iqtree2,
   hmmscan and the untracked data/ directory. Injecting keeps the frozen
   renders valid, and keeps prose files free of navigation markup.

   Step content comes from assets/hail-mary-data.js (CHAPTER_STEPS), which is
   seed text pending author review — see the warning at the top of that file.
   =========================================================================== */

(function () {
  "use strict";

  var D = window.HailMary;
  if (!D) return;

  var REPO = "https://github.com/MintyHybrid/project-hail-mary";

  // Set to true once CHAPTER_STEPS has been reviewed and rewritten, to drop
  // the "draft summary" notice under each step panel. Reviewed by the author
  // against the chapter text: the empirical figures match index.qmd and the
  // open questions (Ch.2 GARD screen, Ch.6 HGT origin, Ch.8/9 model-not-result)
  // are stated as open, consistent with the preprint's hedging.
  var STEPS_REVIEWED = true;

  function el(tag, cls, parent) {
    var node = document.createElement(tag);
    if (cls) node.className = cls;
    if (parent) parent.appendChild(node);
    return node;
  }

  /* Which chapter is this page? Matched on the rendered filename. */
  function currentChapter() {
    var path = window.location.pathname;
    var file = path.substring(path.lastIndexOf("/") + 1) || "index.html";
    for (var i = 0; i < D.CHAPTERS.length; i++) {
      var href = D.CHAPTERS[i].href;
      var base = href.substring(href.lastIndexOf("/") + 1);
      if (base === file) return D.CHAPTERS[i];
    }
    return null;
  }

  function boot() {
    var chapter = currentChapter();
    // the landing page and the appendix have no stepper
    if (!chapter || chapter.id === "A") return;

    var steps = D.CHAPTER_STEPS[chapter.id];
    if (!steps) return;

    var main = document.querySelector("#quarto-document-content") ||
               document.querySelector("main.content") ||
               document.querySelector("main");
    if (!main) return;

    var accent = chapter.role === "dest"
      ? "var(--hm-gold)" : "var(--hm-spark)";

    /* ------------------------------------------------------- workbench -- */
    var wb = el("div", "hm-workbench");
    wb.style.setProperty("--step-accent", accent);

    // breadcrumb + audience control
    var bar = el("div", "hm-workbench-bar", wb);
    var crumb = el("div", "hm-crumb", bar);
    var home = el("a", null, crumb);
    home.href = "../index.html";
    home.textContent = "← HOME";
    var sep = el("span", "sep", crumb); sep.textContent = "/";
    var chLabel = el("span", "ch", crumb);
    chLabel.textContent = "CH. " + chapter.n;

    var seg = el("div", "hm-segmented", bar);
    ["expert", "plain"].forEach(function (mode) {
      var b = el("button", null, seg);
      b.type = "button";
      b.textContent = mode;
      b.setAttribute("data-hm-audience", mode);
      b.addEventListener("click", function () { D.setAudience(mode); });
    });

    // stepper
    var stepper = el("div", "hm-stepper", wb);
    var stepEls = [];
    D.STEP_LABELS.forEach(function (label, i) {
      var step = el("div", "step", stepper);
      var btn = el("button", null, step);
      btn.type = "button";
      var dot = el("span", "dot", btn);
      dot.textContent = String(i + 1);
      var lbl = el("span", "lbl", btn);
      lbl.textContent = label;
      if (i < D.STEP_LABELS.length - 1) el("span", "line", step);
      btn.addEventListener("click", function () { setStep(i); });
      stepEls.push(step);
    });

    // panel
    var panel = el("div", "hm-step-panel", wb);
    var kicker = el("div", "step-kicker", panel);
    var headline = el("p", "step-headline", panel);
    var detail = el("p", "step-detail", panel);
    if (!STEPS_REVIEWED) {
      var note = el("div", "hm-step-review", panel);
      note.textContent =
        "⚑ draft summary — this stepper paraphrases the chapter below; " +
        "the chapter text is authoritative.";
    }

    var active = 0;

    function setStep(i) {
      active = i;
      render();
    }

    function render() {
      var audience = D.getAudience();
      var pair = steps[active] || ["", ""];

      stepEls.forEach(function (s, i) {
        s.classList.toggle("is-active", i === active);
      });
      kicker.textContent = D.STEP_LABELS[active];
      headline.textContent = pair[0] || "";

      // plain mode shows the headline alone; expert adds the detail
      var hasDetail = audience === "expert" && pair[1];
      detail.textContent = pair[1] || "";
      detail.classList.toggle("is-hidden", !hasDetail);

      Array.prototype.forEach.call(
        seg.querySelectorAll("[data-hm-audience]"), function (b) {
          b.setAttribute("aria-pressed",
            String(b.getAttribute("data-hm-audience") === audience));
        });
    }

    document.addEventListener("hm:audience", render);
    setStep(0);

    // insert directly under the title block
    var header = main.querySelector("#title-block-header");
    if (header && header.nextSibling) {
      main.insertBefore(wb, header.nextSibling);
    } else {
      main.insertBefore(wb, main.firstChild);
    }

    /* ------------------------------------------------------- downloads -- */
    var srcName = chapter.href
      .replace(/^notebooks\//, "")
      .replace(/\.html$/, ".qmd");

    var dl = el("div", "hm-downloads");

    var srcLink = el("a", "hm-btn", dl);
    srcLink.href = REPO + "/blob/main/notebooks/" + srcName;
    srcLink.rel = "noopener";
    srcLink.textContent = "⭳ chapter source (.qmd)";

    var dataLink = el("a", "hm-btn", dl);
    dataLink.href = "../data-code.html";
    dataLink.textContent = "⭳ data & code for this chapter";

    var pdfLink = el("button", "hm-btn", dl);
    pdfLink.type = "button";
    pdfLink.textContent = "⭳ export chapter as PDF";
    pdfLink.addEventListener("click", function () { window.print(); });

    main.appendChild(dl);
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", boot);
  } else {
    boot();
  }
})();

/* ===========================================================================
   r3dmol background theme-sync.

   The 3D structure viewers (r3dmol / 3Dmol.js) are static htmlwidgets: their
   background is baked at render time (viewer_spec backgroundColor in the .qmd),
   so a frozen widget can't follow the light/dark toggle the way CSS can, and
   the value baked from the old palette (#EFF1EC) shows as a light box on the
   dark theme.

   The binding exposes each viewer's methods on the widget element as
   `el.widget`, so we can recolour at runtime instead: read the active theme's
   --hm-bg-deep token and push it into every viewer, on load and on every
   theme switch. This keeps the 3D canvas on the same inset colour as the
   other figure plates in both themes, with no re-render of the frozen chapter.
   =========================================================================== */
(function () {
  "use strict";

  function tokenHex(name, fallback) {
    var v = getComputedStyle(document.documentElement)
      .getPropertyValue(name).trim();
    var m = /#?([0-9a-f]{6})/i.exec(v);
    return m ? parseInt(m[1], 16) : fallback;
  }

  function syncViewers() {
    var hex = tokenHex("--hm-bg-deep", 0x181c13);
    var els = document.querySelectorAll(".r3dmol");
    for (var i = 0; i < els.length; i++) {
      var w = els[i].widget;
      if (w && typeof w.setBackgroundColor === "function") {
        try { w.setBackgroundColor({ hex: hex, alpha: 1 }); } catch (e) { /* not ready yet */ }
      }
    }
  }

  // Nothing to do on pages without a viewer.
  if (!document.querySelector(".r3dmol")) return;

  // htmlwidgets render on window load; retry a couple of times in case the
  // viewer binds a beat later, then keep in step with the theme toggle.
  function scheduleSync() {
    syncViewers();
    setTimeout(syncViewers, 300);
    setTimeout(syncViewers, 1200);
  }
  if (document.readyState === "complete") scheduleSync();
  else window.addEventListener("load", scheduleSync);

  // Quarto's scheme toggle swaps a class on <body>; re-read the token after it
  // settles so the viewer follows light <-> dark.
  var mo = new MutationObserver(function () { setTimeout(syncViewers, 50); });
  mo.observe(document.body, { attributes: true, attributeFilter: ["class"] });
})();
