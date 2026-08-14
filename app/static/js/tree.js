// Minimal dependency-free Newick parser + SVG phylogram renderer.
// Adapted from AI_morpho2's phylogeny tree viewer, generalized (no per-tip
// fragment/marker coloring — this app fetches a single marker per job).

function parseNewick(s) {
  s = s.trim().replace(/;$/, '');
  s = s.replace(/\[[^\]]*\]/g, ''); // strip NHX/FigTree comments
  let i = 0;
  function parse() {
    const node = {children: []};
    if (s[i] === '(') {
      i++;
      node.children.push(parse());
      while (s[i] === ',') { i++; node.children.push(parse()); }
      i++; // ')'
    }
    let name = '';
    while (i < s.length && s[i] !== ',' && s[i] !== ')' && s[i] !== ':' && s[i] !== ';')
      name += s[i++];
    node.name = name.trim().replace(/^['"]|['"]$/g, '');
    if (s[i] === ':') {
      i++;
      let len = '';
      while (i < s.length && s[i] !== ',' && s[i] !== ')' && s[i] !== ';') len += s[i++];
      node.length = Math.max(0, parseFloat(len) || 0);
    } else {
      node.length = 0;
    }
    return node;
  }
  try { return parse(); } catch (e) { return null; }
}

function getLeaves(node) {
  if (!node.children || !node.children.length) return [node];
  return node.children.flatMap(getLeaves);
}

function assignDepths(node, d) {
  node.x = d;
  (node.children || []).forEach(c => assignDepths(c, d + (c.length || 1)));
}

function formatTipLabel(name) {
  if (!name || name.includes('!') || name.includes('=')) return '';
  if (name.includes('|')) name = name.split('|').pop();
  name = name.replace(/^_R_/i, '');
  return name.replace(/_/g, ' ').trim();
}

function renderPhylogram(newick, svgEl) {
  while (svgEl.firstChild) svgEl.removeChild(svgEl.firstChild);
  const tree = parseNewick(newick);
  if (!tree) {
    const t = document.createElementNS('http://www.w3.org/2000/svg', 'text');
    t.setAttribute('x', '10'); t.setAttribute('y', '20'); t.setAttribute('fill', 'red');
    t.textContent = 'Could not parse Newick string.';
    svgEl.appendChild(t);
    return;
  }

  const leaves = getLeaves(tree);
  const N = leaves.length;

  const ROW_H = 20, PAD_TOP = 24, PAD_LEFT = 24, LABEL_GAP = 8, MAX_LABEL = 240;
  const containerW = (svgEl.parentElement && svgEl.parentElement.clientWidth) || 920;
  const SVG_W = Math.max(containerW - 4, PAD_LEFT + LABEL_GAP + MAX_LABEL + 120);
  const BRANCH_W = SVG_W - PAD_LEFT - LABEL_GAP - MAX_LABEL - 10;
  const SVG_H = N * ROW_H + PAD_TOP * 2;

  svgEl.setAttribute('width', SVG_W);
  svgEl.setAttribute('height', SVG_H);
  svgEl.setAttribute('viewBox', `0 0 ${SVG_W} ${SVG_H}`);

  leaves.forEach((leaf, idx) => { leaf.leafY = PAD_TOP + idx * ROW_H + ROW_H / 2; });
  assignDepths(tree, 0);
  const maxDepth = Math.max(...leaves.map(l => l.x)) || 1;
  const xScale = v => PAD_LEFT + (v / maxDepth) * BRANCH_W;

  function assignY(node) {
    if (!node.children || !node.children.length) return;
    node.children.forEach(assignY);
    node.leafY = node.children.reduce((s, c) => s + c.leafY, 0) / node.children.length;
  }
  assignY(tree);

  const NS = 'http://www.w3.org/2000/svg';
  function svgLine(x1, y1, x2, y2, stroke, sw, dash) {
    const el = document.createElementNS(NS, 'line');
    el.setAttribute('x1', x1.toFixed(1)); el.setAttribute('y1', y1.toFixed(1));
    el.setAttribute('x2', x2.toFixed(1)); el.setAttribute('y2', y2.toFixed(1));
    el.setAttribute('stroke', stroke || '#444');
    el.setAttribute('stroke-width', sw || 1.2);
    if (dash) el.setAttribute('stroke-dasharray', dash);
    svgEl.appendChild(el);
  }
  function svgText(x, y, content, opts) {
    if (!content) return;
    const el = document.createElementNS(NS, 'text');
    el.setAttribute('x', x.toFixed(1)); el.setAttribute('y', y.toFixed(1));
    el.setAttribute('dominant-baseline', 'middle');
    el.setAttribute('font-size', (opts && opts.size) || '11');
    el.setAttribute('font-family', 'Inter, sans-serif');
    el.setAttribute('fill', (opts && opts.fill) || '#222');
    if (opts && opts.italic) el.setAttribute('font-style', 'italic');
    if (opts && opts.anchor) el.setAttribute('text-anchor', opts.anchor);
    el.textContent = content;
    svgEl.appendChild(el);
  }

  function supportNum(node) {
    if (!node.name) return null;
    const m = String(node.name).match(/^\s*-?\d+(\.\d+)?\s*$/);
    return m ? parseFloat(node.name) : null;
  }
  function supColor(s) {
    const pct = (s <= 1) ? s * 100 : s;
    return pct >= 70 ? '#2f6f4e' : (pct >= 40 ? '#a1751f' : '#a13a3a');
  }
  function supLabel(s) { return (s <= 1) ? s.toFixed(2) : String(Math.round(s)); }

  function drawNode(node) {
    const nx = xScale(node.x);
    if (node.children && node.children.length) {
      const yTop = Math.min(...node.children.map(c => c.leafY));
      const yBot = Math.max(...node.children.map(c => c.leafY));
      svgLine(nx, yTop, nx, yBot);
      node.children.forEach(c => {
        const s = (c.children && c.children.length) ? supportNum(c) : null;
        const col = (s !== null) ? supColor(s) : '#8e8b82';
        svgLine(nx, c.leafY, xScale(c.x), c.leafY, col, (s !== null) ? 2 : undefined);
        if (s !== null) {
          const midX = (nx + xScale(c.x)) / 2;
          svgText(midX, c.leafY - 4, supLabel(s), {size: '9.5', fill: col, anchor: 'middle'});
        }
        drawNode(c);
      });
    } else {
      const labelX = PAD_LEFT + BRANCH_W + LABEL_GAP;
      svgLine(nx, node.leafY, labelX - 2, node.leafY, '#d8d2c6', 0.8, '2,3');
      const label = formatTipLabel(node.name);
      svgText(labelX, node.leafY, label, {italic: true, size: '11.5', fill: '#252523'});
    }
  }
  drawNode(tree);

  const scaleLen = maxDepth * 0.1;
  const scaleBarW = (scaleLen / maxDepth) * BRANCH_W;
  const sbY = SVG_H - 12;
  svgLine(PAD_LEFT, sbY, PAD_LEFT + scaleBarW, sbY, '#6c6a64', 1.5);
  svgLine(PAD_LEFT, sbY - 4, PAD_LEFT, sbY + 4, '#6c6a64', 1.5);
  svgLine(PAD_LEFT + scaleBarW, sbY - 4, PAD_LEFT + scaleBarW, sbY + 4, '#6c6a64', 1.5);
  svgText(PAD_LEFT + scaleBarW / 2, sbY - 8, scaleLen.toPrecision(2),
          {size: '9', fill: '#6c6a64', anchor: 'middle'});

  svgText(PAD_LEFT, 12, 'branch support:', {size: '9', fill: '#6c6a64'});
  svgText(PAD_LEFT + 90, 12, '\u226570 high', {size: '9', fill: '#2f6f4e'});
  svgText(PAD_LEFT + 152, 12, '40\u201369 med', {size: '9', fill: '#a1751f'});
  svgText(PAD_LEFT + 216, 12, '<40 low', {size: '9', fill: '#a13a3a'});
}
