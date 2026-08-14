"""Generalized phylogeny pipeline: NCBI fetch -> Galaxy MAFFT/trimAl align+trim ->
local NJ preview tree -> optional Galaxy RAxML-NG bootstrap ML tree -> optional
outgroup rerooting.

Every NCBI/Galaxy credential is an explicit function argument, supplied by the
caller for each call — nothing here reads a server-side secret and nothing
writes a credential to disk. Adapted from AI_morpho2's app/routes/phylogeny.py,
stripped of every piece specific to a Specimen/Project morphometrics database
(candidate-sequence review, "learned decisions", multi-fragment concatenation,
hardcoded Gyrodactylidae outgroup families).
"""
import json
import os
import re
import threading
import time

import requests

# ── Marker presets ────────────────────────────────────────────────────────────
# code -> (display label, default NCBI query, default minimum sequence length)
MARKER_PRESETS = {
    'COI': {
        'label': 'COI — cytochrome c oxidase I (animal DNA barcode)',
        'query': '(cytochrome c oxidase subunit 1[All Fields] OR cytochrome oxidase '
                 'subunit I[All Fields] OR COI[All Fields] OR COX1[All Fields])',
        'min_length': 500,
    },
    '16S': {
        'label': '16S ribosomal RNA (bacteria/archaea, animal mitochondrial)',
        'query': '16S ribosomal RNA[All Fields]',
        'min_length': 400,
    },
    '18S': {
        'label': '18S ribosomal RNA (broad eukaryote phylogeny)',
        'query': '(18S ribosomal RNA[All Fields] OR small subunit ribosomal '
                 'RNA[All Fields]) NOT (internal transcribed spacer[All Fields])',
        'min_length': 400,
    },
    '28S': {
        'label': '28S ribosomal RNA (broad eukaryote phylogeny)',
        'query': '(28S ribosomal RNA[All Fields] OR large subunit ribosomal '
                 'RNA[All Fields]) NOT (internal transcribed spacer[All Fields])',
        'min_length': 400,
    },
    'ITS': {
        'label': 'ITS — internal transcribed spacer (fungi, plants, close animal relatives)',
        'query': '(internal transcribed spacer[All Fields] OR ITS[All Fields]) '
                 'NOT (18S[All Fields] OR 28S[All Fields])',
        'min_length': 300,
    },
    'cytb': {
        'label': 'cytb — cytochrome b (animal mitochondrial phylogeny)',
        'query': '(cytochrome b[All Fields] OR cytb[All Fields] OR CYTB[All Fields])',
        'min_length': 400,
    },
    'matK': {
        'label': 'matK — maturase K (plant DNA barcode)',
        'query': '(maturase K[All Fields] OR matK[All Fields])',
        'min_length': 400,
    },
    'rbcL': {
        'label': 'rbcL — RuBisCO large subunit (plant DNA barcode)',
        'query': '(ribulose-1,5-bisphosphate carboxylase[All Fields] OR rbcL[All Fields])',
        'min_length': 400,
    },
}


# ── NCBI helpers ──────────────────────────────────────────────────────────────
# Bio.Entrez keeps its credentials as module-global state, so concurrent jobs
# from different users (different NCBI email/key) must not interleave calls —
# serialize with a lock rather than risk one user's request using another's key.
_NCBI_LOCK = threading.Lock()


def _entrez_setup(email, api_key):
    from Bio import Entrez
    Entrez.email = email
    Entrez.api_key = api_key or None
    Entrez.tool = 'phylogen'
    return 0.11 if api_key else 0.4  # NCBI: 10 req/s with a key, 3 req/s without


def ncbi_search(term, email, api_key, retmax=5000):
    from Bio import Entrez
    with _NCBI_LOCK:
        delay = _entrez_setup(email, api_key)
        last_exc = None
        for attempt in range(5):
            try:
                h = Entrez.esearch(db='nuccore', term=term, retmax=retmax)
                result = Entrez.read(h)
                h.close()
                time.sleep(delay)
                return result['IdList'], int(result['Count'])
            except Exception as exc:
                last_exc = exc
                is_429 = '429' in str(exc) or 'Too Many Requests' in str(exc)
                if attempt < 4:
                    time.sleep(5 * (attempt + 1) if is_429 else 2)
        raise RuntimeError(f'NCBI search failed: {last_exc}')


def ncbi_fetch_batch(ids, email, api_key, batch_size=200):
    from Bio import Entrez, SeqIO
    records = {}
    with _NCBI_LOCK:
        delay = _entrez_setup(email, api_key)
        for i in range(0, len(ids), batch_size):
            batch = ids[i:i + batch_size]
            last_exc = None
            for attempt in range(5):
                try:
                    h = Entrez.efetch(db='nuccore', id=','.join(batch),
                                       rettype='fasta', retmode='text')
                    for rec in SeqIO.parse(h, 'fasta'):
                        records[rec.id] = rec
                    h.close()
                    break
                except Exception as exc:
                    last_exc = exc
                    is_429 = '429' in str(exc) or 'Too Many Requests' in str(exc)
                    if attempt < 4:
                        time.sleep(5 * (attempt + 1) if is_429 else 2)
                    else:
                        raise RuntimeError(f'NCBI fetch failed: {last_exc}')
            time.sleep(delay)
    return records


def _parse_species_name(description):
    """'ACC.1 Genus species strain X' -> 'Genus_species'."""
    parts = description.split()
    if len(parts) >= 3:
        return f'{parts[1]}_{parts[2]}'
    if len(parts) >= 2:
        return parts[1]
    return parts[0] if parts else 'unknown'


def process_records(records, min_length=400, max_length_factor=2.0):
    """One record per species: the longest sequence not wildly longer than the
    set's mean (an outlier guard against e.g. a whole mitogenome dwarfing single-
    gene amplicons), never dropping a species outright. Renames to
    'accession|Genus_species'."""
    from Bio.SeqRecord import SeqRecord

    records = {k: v for k, v in records.items() if len(v.seq) >= min_length}
    if not records:
        return []

    by_sp_seq = {}
    for rec in records.values():
        sp = _parse_species_name(rec.description)
        key = (sp, str(rec.seq).upper())
        if key not in by_sp_seq or len(rec.seq) > len(by_sp_seq[key].seq):
            by_sp_seq[key] = rec
    records = {r.id: r for r in by_sp_seq.values()}

    lengths = [len(r.seq) for r in records.values()]
    max_allowed = max_length_factor * (sum(lengths) / len(lengths))

    by_species = {}
    for rec in records.values():
        sp = _parse_species_name(rec.description)
        by_species.setdefault(sp, []).append(rec)
    for sp in by_species:
        by_species[sp].sort(key=lambda r: len(r.seq), reverse=True)

    result = []
    for sp, recs in by_species.items():
        chosen = next((r for r in recs if len(r.seq) <= max_allowed), None)
        if chosen is None:
            chosen = min(recs, key=lambda r: len(r.seq))
        result.append(SeqRecord(chosen.seq, id=f'{chosen.id}|{sp}', name='', description=''))
    return result


def fetch_species_list(species_names, gene_query, email, api_key, min_length, progress=None):
    """One targeted NCBI search per species name. Returns (records, found, missing)."""
    all_records = {}
    found, missing = [], []
    names = [s.strip() for s in species_names if s.strip()]
    for i, sp in enumerate(names):
        if progress:
            progress(f'Searching NCBI for {sp} ({i + 1}/{len(names)})…')
        query = f'"{sp}"[Organism] AND ({gene_query})'
        try:
            ids, _ = ncbi_search(query, email, api_key, retmax=50)
        except Exception:
            ids = []
        if not ids:
            missing.append(sp)
            continue
        recs = ncbi_fetch_batch(ids, email, api_key)
        all_records.update(recs)
        found.append(sp)
    processed = process_records(all_records, min_length)
    return processed, found, missing


def fetch_taxon(taxon, gene_query, email, api_key, min_length, max_species=40, progress=None):
    """Whole-taxon discovery: search a higher taxon name, cap to max_species."""
    relaxed_q = re.sub(r'\s*NOT\s*\([^)]*\)', '', gene_query or '').strip()
    query = f'"{taxon}"[Organism] AND ({relaxed_q})'
    if progress:
        progress(f'Searching NCBI: {query}')
    ids, count = ncbi_search(query, email, api_key, retmax=2000)
    if progress:
        progress(f'Found {count} matching records; downloading {len(ids)}…')
    records = ncbi_fetch_batch(ids, email, api_key)
    processed = process_records(records, min_length)
    n_total_species = len(processed)
    processed.sort(key=lambda r: r.id.split('|', 1)[1] if '|' in r.id else r.id)
    capped = processed[:max_species]
    species_found = [r.id.split('|', 1)[1].replace('_', ' ') for r in capped if '|' in r.id]
    return capped, species_found, n_total_species


def summarize_fetch(records):
    """Review-table summary for the pre-align confirmation screen: one row per
    recovered sequence (species, accession, length)."""
    found = []
    for r in records:
        if '|' in r.id:
            acc, sp = r.id.split('|', 1)
        else:
            acc, sp = r.id, r.id
        found.append({'species': sp.replace('_', ' '), 'accession': acc, 'length': len(r.seq)})
    found.sort(key=lambda row: row['species'])
    return found


def filter_fasta_exclude(input_path, exclude_accessions, output_path):
    """Rewrite input_path to output_path, dropping any record whose header
    contains one of exclude_accessions (user-rejected sequences from the
    review screen). Returns the list of removed header lines."""
    exclude = {a.strip() for a in (exclude_accessions or []) if a.strip()}
    if not exclude:
        import shutil
        shutil.copyfile(input_path, output_path)
        return []
    removed, skip, buf = [], False, []
    with open(input_path) as fin:
        for line in fin:
            if line.startswith('>'):
                skip = any(acc in line for acc in exclude)
                if skip:
                    removed.append(line.strip())
            if not skip:
                buf.append(line)
    with open(output_path, 'w') as fout:
        fout.writelines(buf)
    return removed


def write_fasta(records, path):
    from Bio import SeqIO
    with open(path, 'w') as fh:
        SeqIO.write(records, fh, 'fasta')


def count_fasta(path):
    if not os.path.exists(path):
        return 0
    with open(path) as fh:
        return sum(1 for line in fh if line.startswith('>'))


def orient_fasta_by_reference(in_path, k=8):
    """Reverse-complement sequences whose k-mers match the longest reference
    sequence better in RC orientation. Galaxy's MAFFT wrapper has no
    --adjustdirection, so this runs before upload. Rewrites in_path in place
    only if something flipped. Returns list of flipped ids (best-effort)."""
    flipped_ids = []
    try:
        from Bio import SeqIO
        from Bio.Seq import Seq

        def kmers(s):
            return {s[i:i + k] for i in range(len(s) - k + 1)}

        records = list(SeqIO.parse(in_path, 'fasta'))
        if len(records) < 2:
            return flipped_ids
        ref = max(records, key=lambda r: len(r.seq))
        ref_k = kmers(str(ref.seq).upper())
        if not ref_k:
            return flipped_ids
        for rec in records:
            if rec is ref:
                continue
            s = str(rec.seq).upper()
            if len(s) < k:
                continue
            fwd = sum(1 for km in kmers(s) if km in ref_k)
            rc = sum(1 for km in kmers(str(Seq(s).reverse_complement())) if km in ref_k)
            if rc > fwd:
                rec.seq = rec.seq.reverse_complement()
                rec.description = ''
                flipped_ids.append(rec.id)
        if flipped_ids:
            SeqIO.write(records, in_path, 'fasta')
    except Exception:
        return []
    return flipped_ids


# ── Local NJ preview tree ─────────────────────────────────────────────────────

def local_nj_tree(trimmed_path):
    """Neighbor-Joining tree from the trimmed alignment (identity distance,
    pure Biopython — no external service). Returns a Newick string."""
    from Bio import AlignIO, Phylo
    from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
    import io

    alignment = AlignIO.read(trimmed_path, 'fasta')
    calc = DistanceCalculator('identity')
    dm = calc.get_distance(alignment)
    nj_tree = DistanceTreeConstructor().nj(dm)
    buf = io.StringIO()
    Phylo.write(nj_tree, buf, 'newick')
    return buf.getvalue().strip()


# ── Substitution model estimation ─────────────────────────────────────────────
# No local ModelTest-NG binary is available on Railway (same constraint as
# MAFFT/trimAl), so this is a fast, always-available heuristic — base
# composition homogeneity + transition/transversion ratio — rather than a full
# BIC search across candidate models. Good enough to pick a sensible RAxML-NG
# model string per fragment; NOT a substitute for a real ModelTest-NG run.
_PURINES = frozenset('AG')


def estimate_model(trimmed_path, max_pairs=200):
    """Returns a model string ('JC+G' | 'K80+G' | 'HKY+G' | 'GTR+G') estimated
    from base-frequency skew and the observed ts/tv ratio in the alignment."""
    from Bio import AlignIO
    try:
        aln = AlignIO.read(trimmed_path, 'fasta')
    except Exception:
        return 'GTR+G'
    seqs = [str(r.seq).upper() for r in aln]
    if len(seqs) < 2:
        return 'GTR+G'

    counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    for s in seqs:
        for b in s:
            if b in counts:
                counts[b] += 1
    total = sum(counts.values())
    if total == 0:
        return 'GTR+G'
    max_dev = max(abs(counts[b] / total - 0.25) for b in counts)

    ts = tv = n_pairs = 0
    for i in range(len(seqs)):
        if n_pairs > max_pairs:
            break
        for j in range(i + 1, len(seqs)):
            n_pairs += 1
            for a, b in zip(seqs[i], seqs[j]):
                if a not in 'ACGT' or b not in 'ACGT' or a == b:
                    continue
                if (a in _PURINES) == (b in _PURINES):
                    ts += 1
                else:
                    tv += 1
            if n_pairs > max_pairs:
                break
    ts_tv = (ts / tv) if tv else (2.0 if ts else 1.0)

    even_freqs = max_dev < 0.04
    strong_titv = ts_tv > 1.3
    if even_freqs and not strong_titv:
        base = 'JC'
    elif even_freqs and strong_titv:
        base = 'K80'
    elif not even_freqs and strong_titv:
        base = 'HKY'
    else:
        base = 'GTR'
    return f'{base}+G'


# ── Fragment concatenation (multi-marker) ─────────────────────────────────────

def concatenate_fragments(frag_paths, out_path):
    """frag_paths: list of (trimmed_path, fragment_code). Concatenates trimmed
    alignments; a taxon missing one fragment gets all-gap columns for it.
    Returns (n_taxa, partition_spec, species_coverage): partition_spec is
    [{'name','start','end'}] (1-based inclusive column ranges), species_coverage
    maps 'genus species' (lowercase) -> '+'.join(fragment codes present)."""
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    loaded = []
    for path, code in frag_paths:
        by_sp = {}
        for r in SeqIO.parse(path, 'fasta'):
            sp = r.id.split('|', 1)[1] if '|' in r.id else r.id
            by_sp[sp] = r
        width = len(next(iter(by_sp.values())).seq) if by_sp else 0
        loaded.append((code, by_sp, width))

    all_sp = sorted({sp for _, by_sp, _ in loaded for sp in by_sp})
    records, coverage = [], {}
    for sp in all_sp:
        seq_parts, marks, rec_id = [], [], None
        for code, by_sp, width in loaded:
            r = by_sp.get(sp)
            if r is not None:
                seq_parts.append(str(r.seq))
                marks.append(code)
                rec_id = rec_id or r.id
            else:
                seq_parts.append('-' * width)
        records.append(SeqRecord(Seq(''.join(seq_parts)), id=rec_id or sp, name='', description=''))
        coverage[sp.replace('_', ' ').lower()] = '+'.join(marks)

    write_fasta(records, out_path)
    spec, cursor = [], 1
    for code, _by_sp, width in loaded:
        spec.append({'name': code, 'start': cursor, 'end': cursor + width - 1})
        cursor += width
    return len(records), spec, coverage


# ── Model-corrected partitioned NJ ────────────────────────────────────────────

def _p_distance_correction(model):
    """Map an estimated/RAxML-NG model string to a distance-correction family:
    'JC' (Jukes-Cantor) for simple models, 'K80' (Kimura-2P) when the model
    distinguishes transitions/transversions (K80/HKY/GTR/...)."""
    if not model:
        return 'JC'
    m = model.upper()
    if any(t in m for t in ('K80', 'K2P', 'HKY', 'TN', 'TIM', 'TVM', 'GTR', 'SYM')):
        return 'K80'
    return 'JC'


def _corrected_pair_distance(si, sj, fam):
    import math
    sites = ts = tv = 0
    for a, b in zip(si, sj):
        if a not in 'ACGT' or b not in 'ACGT':
            continue
        sites += 1
        if a == b:
            continue
        if (a in _PURINES) == (b in _PURINES):
            ts += 1
        else:
            tv += 1
    if sites == 0:
        return 0.0, 0
    if fam == 'K80':
        P, Q = ts / sites, tv / sites
        try:
            d = -0.5 * math.log(1 - 2 * P - Q) - 0.25 * math.log(1 - 2 * Q)
        except ValueError:
            d = 2.0
    else:
        p = (ts + tv) / sites
        try:
            d = -0.75 * math.log(1 - 4 / 3 * p)
        except ValueError:
            d = 2.0
    if not (d == d) or d < 0:
        d = 2.0
    return d, sites


def local_nj_tree_partitioned(concat_path, partition_spec, models):
    """NJ tree from a concatenated alignment using model-corrected pairwise
    distances per partition (models: {fragment_code: model_string}), combined
    weighted by comparable sites per partition."""
    from Bio import AlignIO, Phylo
    from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
    import io

    aln = AlignIO.read(concat_path, 'fasta')
    ids = [rec.id for rec in aln]
    seqs = [str(rec.seq).upper() for rec in aln]
    n = len(ids)
    parts = [(p['start'] - 1, p['end'], _p_distance_correction(models.get(p['name'])))
             for p in (partition_spec or [])]
    if not parts:
        parts = [(0, len(seqs[0]) if seqs else 0, 'JC')]

    matrix = [[0.0] * (i + 1) for i in range(n)]
    for i in range(n):
        for j in range(i):
            dacc = wacc = 0.0
            for s0, s1, fam in parts:
                d, w = _corrected_pair_distance(seqs[i][s0:s1], seqs[j][s0:s1], fam)
                dacc += d * w
                wacc += w
            matrix[i][j] = (dacc / wacc) if wacc else 1.0
    dm = DistanceMatrix(ids, matrix)
    nj_tree = DistanceTreeConstructor().nj(dm)
    buf = io.StringIO()
    Phylo.write(nj_tree, buf, 'newick')
    return buf.getvalue().strip()


# ── Rerooting ─────────────────────────────────────────────────────────────────

def reroot_newick(newick, names):
    """Reroot at the MRCA of tips matching `names` (substring, case-insensitive
    on 'genus species' form). Falls back to midpoint rooting if none match or
    `names` is empty. Returns (new_newick, rooted_on, not_found)."""
    from Bio import Phylo
    from io import StringIO

    tree = Phylo.read(StringIO(newick), 'newick')
    want = [n for n in (names or []) if n and n.strip()]
    matched, not_found = [], []
    if want:
        terminals = tree.get_terminals()
        for name in want:
            t = next((t for t in terminals if t.name and (
                name.lower() in t.name.lower().replace('_', ' ') or
                t.name.lower().replace('_', ' ') in name.lower())), None)
            if t:
                matched.append(t)
            else:
                not_found.append(name)
    if not matched:
        tree.root_at_midpoint()
        rooted_on = 'midpoint (no outgroup match)' if want else 'midpoint'
    elif len(matched) == 1:
        tree.root_with_outgroup(matched[0])
        rooted_on = matched[0].name
    else:
        tree.root_with_outgroup(tree.common_ancestor(matched))
        rooted_on = ', '.join(m.name for m in matched)
    buf = StringIO()
    Phylo.write(tree, buf, 'newick')
    return buf.getvalue().strip(), rooted_on, not_found


# ── Galaxy (usegalaxy.eu) REST helpers ────────────────────────────────────────

def _galaxy_base():
    from flask import current_app
    return current_app.config['GALAXY_BASE_URL']


def _galaxy_headers(api_key):
    return {'x-api-key': api_key, 'Accept': 'application/json'}


def _galaxy_raise_for_status(r, what):
    if r.status_code < 400:
        return
    detail = ''
    try:
        body = r.json()
        detail = body.get('err_msg') or body.get('message') or json.dumps(body)
    except ValueError:
        detail = (r.text or '')[:500]
    raise RuntimeError(f'Galaxy {what} failed ({r.status_code}): {detail}')


def galaxy_create_history(api_key, name='phylogen'):
    base = _galaxy_base()
    r = requests.post(f'{base}/api/histories',
                       headers={**_galaxy_headers(api_key), 'Content-Type': 'application/json'},
                       json={'name': name}, timeout=30)
    r.raise_for_status()
    return r.json()['id']


def galaxy_upload_file(api_key, history_id, file_path, file_type='fasta'):
    base = _galaxy_base()
    with open(file_path, 'rb') as fh:
        r = requests.post(
            f'{base}/api/tools', headers=_galaxy_headers(api_key),
            data={
                'tool_id': 'upload1',
                'history_id': history_id,
                'inputs': json.dumps({
                    'files_0|NAME': os.path.basename(file_path),
                    'file_count': '1',
                    'file_type': file_type,
                    'dbkey': '?',
                }),
            },
            files={'files_0|file_data': (os.path.basename(file_path), fh)},
            timeout=180,
        )
    _galaxy_raise_for_status(r, f'upload ({file_type})')
    data = r.json()
    dataset_id = data['outputs'][0]['id']
    upload_job_id = data['jobs'][0]['id'] if data.get('jobs') else None
    return dataset_id, upload_job_id


def galaxy_wait_for_job(api_key, job_id, max_wait=300):
    base = _galaxy_base()
    deadline = time.time() + max_wait
    while time.time() < deadline:
        r = requests.get(f'{base}/api/jobs/{job_id}', headers=_galaxy_headers(api_key), timeout=30)
        r.raise_for_status()
        state = r.json().get('state', 'running')
        if state in ('ok', 'error', 'deleted', 'paused'):
            return state
        time.sleep(5)
    return 'timeout'


def galaxy_run_tool(api_key, history_id, tool_id, inputs):
    base = _galaxy_base()
    r = requests.post(
        f'{base}/api/tools',
        headers={**_galaxy_headers(api_key), 'Content-Type': 'application/json'},
        json={'tool_id': tool_id, 'history_id': history_id, 'inputs': inputs}, timeout=60)
    _galaxy_raise_for_status(r, f'tool run ({tool_id})')
    data = r.json()
    if data.get('err_msg'):
        raise RuntimeError(f'Galaxy tool error: {data["err_msg"]}')
    jobs = data.get('jobs') or []
    if not jobs:
        raise RuntimeError(f'Galaxy returned no jobs: {str(data)[:400]}')
    return jobs[0]['id']


def galaxy_download_dataset(api_key, ds_id, dest_path):
    base = _galaxy_base()
    dl = requests.get(f'{base}/api/datasets/{ds_id}/display',
                       headers=_galaxy_headers(api_key), stream=True, timeout=600)
    dl.raise_for_status()
    with open(dest_path, 'wb') as fh:
        for chunk in dl.iter_content(65536):
            fh.write(chunk)
    return dest_path


def galaxy_run_chain(api_key, history_id, tool_id, input_key, dataset_id,
                      extra_params_json, label, max_wait=3600):
    inputs = {input_key: {'src': 'hda', 'id': dataset_id}}
    try:
        extra = json.loads(extra_params_json or '{}')
        if isinstance(extra, dict):
            inputs.update(extra)
    except Exception:
        pass
    gjob = galaxy_run_tool(api_key, history_id, tool_id, inputs)
    state = galaxy_wait_for_job(api_key, gjob, max_wait=max_wait)
    if state != 'ok':
        raise RuntimeError(f'Galaxy {label} job ended in state "{state}".')
    return gjob


def galaxy_pick_fasta_output(api_key, job_id, dest_path):
    base = _galaxy_base()
    r = requests.get(f'{base}/api/jobs/{job_id}/outputs', headers=_galaxy_headers(api_key), timeout=30)
    r.raise_for_status()
    cands = []
    for out in r.json():
        ds = out.get('dataset') or {}
        ds_id = ds.get('id') or out.get('id')
        if ds_id:
            cands.append((ds_id, (out.get('name') or '').lower()))
    cands.sort(key=lambda c: any(k in c[1] for k in ('html', 'report', 'log', 'summary')))
    tmp = dest_path + '.cand'
    for ds_id, _name in cands:
        try:
            galaxy_download_dataset(api_key, ds_id, tmp)
        except Exception:
            continue
        if os.path.exists(tmp) and count_fasta(tmp) > 0:
            os.replace(tmp, dest_path)
            return ds_id, dest_path
    if os.path.exists(tmp):
        try:
            os.remove(tmp)
        except OSError:
            pass
    return None, None


def galaxy_align_trim(api_key, in_path, aligned_out, trimmed_out):
    """MAFFT-align then trimAl-trim on Galaxy. Trimming is best-effort: if
    trimAl drops sequences or produces nothing usable, the untrimmed alignment
    is used instead so no sequence is ever silently lost. Returns
    (aligned_out, trimmed_out, flipped_ids, trim_skipped)."""
    import shutil
    from flask import current_app
    cfg = current_app.config

    flipped = orient_fasta_by_reference(in_path)
    hist = galaxy_create_history(api_key, 'phylogen_AlignTrim')
    ds_id, up_job = galaxy_upload_file(api_key, hist, in_path, 'fasta')
    if up_job:
        st = galaxy_wait_for_job(api_key, up_job, max_wait=600)
        if st != 'ok':
            raise RuntimeError(f'Galaxy upload failed (state: {st}).')

    mafft_job = galaxy_run_chain(
        api_key, hist, cfg['GALAXY_MAFFT_TOOL_ID'], cfg['GALAXY_MAFFT_INPUT_KEY'],
        ds_id, cfg['GALAXY_MAFFT_PARAMS'], 'MAFFT')
    aln_ds, _ = galaxy_pick_fasta_output(api_key, mafft_job, aligned_out)
    if not aln_ds:
        raise RuntimeError('Galaxy MAFFT produced no FASTA output.')

    n_aligned = count_fasta(aligned_out)
    trim_skipped = False
    try:
        trim_job = galaxy_run_chain(
            api_key, hist, cfg['GALAXY_TRIMAL_TOOL_ID'], cfg['GALAXY_TRIMAL_INPUT_KEY'],
            aln_ds, cfg['GALAXY_TRIMAL_PARAMS'], 'trimAl')
        trm_ds, _ = galaxy_pick_fasta_output(api_key, trim_job, trimmed_out)
    except Exception:
        trm_ds = None
    if trm_ds and count_fasta(trimmed_out) < n_aligned:
        trm_ds = None
    if not trm_ds:
        shutil.copyfile(aligned_out, trimmed_out)
        trim_skipped = True
    return aligned_out, trimmed_out, flipped, trim_skipped


# ── Galaxy RAxML-NG (bootstrap ML tree) ───────────────────────────────────────

def _galaxy_find_tool_id(api_key, needle):
    base = _galaxy_base()
    r = requests.get(f'{base}/api/tools', params={'q': needle, 'in_panel': 'false'},
                      headers=_galaxy_headers(api_key), timeout=30)
    r.raise_for_status()
    ids = []

    def _collect(o):
        if isinstance(o, dict):
            tid = o.get('id')
            if isinstance(tid, str):
                ids.append(tid)
            for v in o.values():
                _collect(v)
        elif isinstance(o, list):
            for v in o:
                _collect(v)
    _collect(r.json())
    n = needle.lower()
    matches = [t for t in ids if n in t.lower()]
    seg = [t for t in matches if f'/{n}/' in t.lower()] or matches
    seg.sort(key=lambda t: t.split('/')[-1], reverse=True)
    return seg[0] if seg else None


def _write_raxmlng_partition_file(partition_spec, partition_models, out_path):
    """Write a RAxML-NG partition file: one '<model>, <name> = <start>-<end>'
    line per fragment. Returns out_path, or None if there is nothing to write."""
    lines = []
    for p in (partition_spec or []):
        if p.get('end', 0) < p.get('start', 1):
            continue
        model = (partition_models or {}).get(p['name']) or 'GTR+G'
        lines.append(f"{model}, {p['name']} = {p['start']}-{p['end']}")
    if not lines:
        return None
    with open(out_path, 'w') as fh:
        fh.write('\n'.join(lines) + '\n')
    return out_path


def submit_raxml(fasta_path, api_key, n_bootstraps=1000, model='GTR+G',
                  partition_spec=None, partition_models=None):
    """Upload the (possibly concatenated) alignment and run RAxML-NG (--all: ML
    search + bootstrap + support) on Galaxy. Returns (history_id, galaxy_job_id).

    With more than one partition, a RAxML-NG partition file (one model per
    fragment, from `partition_models`) is uploaded and passed via
    model_type=multi_file; otherwise a single model string is used."""
    from flask import current_app
    tool_id = current_app.config.get('GALAXY_RAXMLNG_TOOL_ID') or None
    if not tool_id:
        tool_id = (_galaxy_find_tool_id(api_key, 'raxmlng')
                   or _galaxy_find_tool_id(api_key, 'raxml_ng')
                   or 'toolshed.g2.bx.psu.edu/repos/iuc/raxmlng/raxmlng/2.0.2+galaxy0')

    history_id = galaxy_create_history(api_key, 'phylogen_RAxML-NG')
    ds_id, up_job = galaxy_upload_file(api_key, history_id, fasta_path, 'fasta')
    if up_job:
        state = galaxy_wait_for_job(api_key, up_job, max_wait=300)
        if state != 'ok':
            raise RuntimeError(f'Galaxy upload job failed (state: {state})')

    inputs = {
        'general_opts|cmdtype|infile': {'src': 'hda', 'id': ds_id},
        'general_opts|cmdtype|command': '--all',
        'general_opts|cmdtype|bs_metric': 'fbp',
        'bootstrap_opts|bs_reps': int(n_bootstraps),
        'bootstrap_opts|bs_mre': 'true',
        'random_seed': 1234567890,
    }

    part_path = os.path.join(os.path.dirname(fasta_path), 'raxmlng_partitions.txt')
    part_file = (_write_raxmlng_partition_file(partition_spec, partition_models, part_path)
                 if partition_spec and len(partition_spec) > 1 else None)
    if part_file:
        part_ds_id, part_job = galaxy_upload_file(api_key, history_id, part_file, 'txt')
        if part_job:
            state = galaxy_wait_for_job(api_key, part_job, max_wait=300)
            if state != 'ok':
                raise RuntimeError(f'Galaxy partition upload failed (state: {state})')
        inputs['general_opts|cmdtype|model|model_type'] = 'multi_file'
        inputs['general_opts|cmdtype|model|model_file'] = {'src': 'hda', 'id': part_ds_id}
        inputs['general_opts|cmdtype|model|brlen_linkage'] = 'scaled'
        inputs['general_opts|cmdtype|model|model_file_auto'] = 'false'
    else:
        single_model = model
        if partition_spec and len(partition_spec) == 1:
            single_model = (partition_models or {}).get(partition_spec[0]['name']) or model
        inputs['general_opts|cmdtype|model|model_type'] = 'single_string'
        inputs['general_opts|cmdtype|model|model_string'] = single_model or 'GTR+G'

    job_id = galaxy_run_tool(api_key, history_id, tool_id, inputs)
    return history_id, job_id


def galaxy_check_status(api_key, job_id):
    """Returns (stage, message) — stage in RUNNING/COMPLETED/FAILED/SUSPENDED."""
    base = _galaxy_base()
    r = requests.get(f'{base}/api/jobs/{job_id}', headers=_galaxy_headers(api_key), timeout=30)
    r.raise_for_status()
    data = r.json()
    state = data.get('state', 'unknown')
    msg = data.get('stderr', '') or data.get('stdout', '') or state
    if state == 'ok':
        return 'COMPLETED', msg
    if state in ('error', 'deleted'):
        return 'FAILED', msg
    if state == 'paused':
        return 'SUSPENDED', msg
    return 'RUNNING', msg


def galaxy_download_results(api_key, job_id, dest_dir):
    base = _galaxy_base()
    os.makedirs(dest_dir, exist_ok=True)
    r = requests.get(f'{base}/api/jobs/{job_id}/outputs', headers=_galaxy_headers(api_key), timeout=30)
    r.raise_for_status()
    downloaded = []
    for out in r.json():
        ds = out.get('dataset') or {}
        ds_id = ds.get('id') or out.get('id')
        name = (out.get('name') or 'output').replace(' ', '_').replace('/', '_')
        if not ds_id:
            continue
        dest = os.path.join(dest_dir, f'{name}.dat')
        try:
            dl = requests.get(f'{base}/api/datasets/{ds_id}/display',
                               headers=_galaxy_headers(api_key), stream=True, timeout=300)
            dl.raise_for_status()
            with open(dest, 'wb') as fh:
                for chunk in dl.iter_content(65536):
                    fh.write(chunk)
            downloaded.append(os.path.basename(dest))
        except Exception:
            pass
    return downloaded


_TREE_NAME_PRIORITY = (
    'bipartitions', 'branchsupportvalues', 'support', 'bipartitionsbranchlabels',
    'bestscoringmltree', 'besttree', 'mltree', 'result', 'bestmodel',
)
_TREE_NAMES_WITH_SUPPORT = {'bipartitions', 'branchsupportvalues', 'support',
                             'bipartitionsbranchlabels'}


def _norm_ds_name(name):
    stem = os.path.splitext(name)[0]
    return re.sub(r'[^a-z0-9]', '', stem.lower())


def _looks_like_newick(path):
    try:
        with open(path, 'r', errors='ignore') as fh:
            head = fh.read(8192)
    except OSError:
        return False
    head = head.lstrip()
    return head.startswith('(') and ')' in head


def find_best_tree(results_dir):
    """(path, has_support) for the best RAxML-NG tree file, preferring one that
    carries bootstrap support as node labels. (None, False) if nothing found."""
    try:
        names = os.listdir(results_dir)
    except OSError:
        return None, False
    norm = {name: _norm_ds_name(name) for name in names}
    for key in _TREE_NAME_PRIORITY:
        for match in (str.endswith, lambda n, k: k in n):
            for name in names:
                if match(norm[name], key):
                    path = os.path.join(results_dir, name)
                    if os.path.isfile(path) and _looks_like_newick(path):
                        return path, key in _TREE_NAMES_WITH_SUPPORT
    return None, False


# ── Divergence-time dating (IQ-TREE + LSD2 on Galaxy) ─────────────────────────
# LSD2 (least-squares dating) is bundled into IQ-TREE and hosted on
# usegalaxy.eu as toolshed.g2.bx.psu.edu/repos/iuc/iqtree/iqtree. Verified
# parameter wiring (fetched live from the Galaxy tool model, not guessed):
#   s                       = alignment dataset
#   t  + n=0                = starting tree + zero search iterations, the
#                             closest equivalent this wrapper exposes to
#                             IQ-TREE's own '-te' (fix topology, don't search)
#   m                       = substitution model string
#   o                       = outgroup taxa (comma-separated), for rooting
#   date_source|select_source = 'dataset', date_source|date = calibration file
#   date_ci, clock_sd       = LSD2 confidence-interval and relaxed-clock knobs
# v1 limitation: one model string is used for the whole (possibly
# concatenated) alignment during dating -- the RAxML-NG step still properly
# partitions per fragment, but LSD2's own topology/rate re-optimization here
# does not; stated plainly in the UI rather than silently approximated.

IQTREE_TOOL_ID = 'toolshed.g2.bx.psu.edu/repos/iuc/iqtree/iqtree/3.1.3+galaxy0'


def _match_tip_name(tip_names, query):
    """Case-insensitive substring match (either direction) of `query` against
    tip_names (each normalized: '_' -> ' '). Returns the matching raw tip name,
    or None."""
    q = (query or '').strip().lower()
    if not q:
        return None
    for t in tip_names:
        norm = t.lower().replace('_', ' ')
        if q in norm or norm in q:
            return t
    return None


def newick_tip_names(newick):
    from Bio import Phylo
    from io import StringIO
    tree = Phylo.read(StringIO(newick), 'newick')
    return [t.name for t in tree.get_terminals() if t.name]


def _fmt_past_age(age_mya):
    """Format a positive Mya value as an LSD2 past-time number: 0 stays '0'
    (avoids a cosmetic '-0.0'), everything else becomes negative."""
    v = abs(float(age_mya))
    return '0' if v == 0 else f'-{v}'


def build_calibration_file(calibrations, tip_names, out_path):
    """calibrations: list of {'taxa': [species name, ...], 'min_age': float|None,
    'max_age': float|None} (millions of years before present). Writes an LSD2
    ancestral-date file: 'tipA,tipB<TAB><lower>:<upper>' on the LSD2 numeric
    time axis (0 = present, more negative = further in the past; the more
    negative bound is written first, matching LSD2's documented ascending-
    order convention). Returns (path_or_None, matched_rows, unmatched_names)."""
    lines, matched_rows, unmatched = [], [], []
    for cal in calibrations or []:
        taxa_in, min_age, max_age = cal.get('taxa') or [], cal.get('min_age'), cal.get('max_age')
        if min_age is None and max_age is None:
            continue
        matched_taxa = []
        for name in taxa_in:
            tip = _match_tip_name(tip_names, name)
            if tip:
                matched_taxa.append(tip)
            else:
                unmatched.append(name)
        if not matched_taxa:
            continue
        if min_age is not None and max_age is not None and float(min_age) == float(max_age):
            date_field = _fmt_past_age(min_age)
        else:
            lo = _fmt_past_age(max_age) if max_age is not None else 'NA'
            hi = _fmt_past_age(min_age) if min_age is not None else 'NA'
            date_field = f'{lo}:{hi}'
        lines.append(f"{','.join(matched_taxa)}\t{date_field}")
        matched_rows.append({'taxa': matched_taxa, 'min_age': min_age, 'max_age': max_age})
    if not lines:
        return None, [], unmatched
    with open(out_path, 'w') as fh:
        fh.write('\n'.join(lines) + '\n')
    return out_path, matched_rows, unmatched


def submit_dating(concat_path, ml_tree_newick, date_file_path, api_key,
                   model='GTR+G', outgroup_names=None, date_ci=100, clock_sd=0.2):
    """Submit an IQ-TREE + LSD2 divergence-dating run on Galaxy. Returns
    (history_id, galaxy_job_id)."""
    history_id = galaxy_create_history(api_key, 'phylogen_Dating')

    aln_ds, up_job = galaxy_upload_file(api_key, history_id, concat_path, 'fasta')
    if up_job:
        state = galaxy_wait_for_job(api_key, up_job, max_wait=300)
        if state != 'ok':
            raise RuntimeError(f'Galaxy alignment upload failed (state: {state})')

    tree_path = concat_path + '.dating_input_tree.nwk'
    with open(tree_path, 'w') as fh:
        fh.write(ml_tree_newick.strip() + '\n')
    tree_ds, tree_job = galaxy_upload_file(api_key, history_id, tree_path, 'nhx')
    if tree_job:
        state = galaxy_wait_for_job(api_key, tree_job, max_wait=300)
        if state != 'ok':
            raise RuntimeError(f'Galaxy tree upload failed (state: {state})')

    date_ds, date_job = galaxy_upload_file(api_key, history_id, date_file_path, 'txt')
    if date_job:
        state = galaxy_wait_for_job(api_key, date_job, max_wait=300)
        if state != 'ok':
            raise RuntimeError(f'Galaxy date-file upload failed (state: {state})')

    inputs = {
        's': {'src': 'hda', 'id': aln_ds},
        't': {'src': 'hda', 'id': tree_ds},
        'n': 0,
        'm': model or 'GTR+G',
        'date_source|select_source': 'dataset',
        'date_source|date': {'src': 'hda', 'id': date_ds},
        'date_ci': int(date_ci),
        'clock_sd': float(clock_sd),
    }
    if outgroup_names:
        inputs['o'] = ','.join(outgroup_names)

    job_id = galaxy_run_tool(api_key, history_id, IQTREE_TOOL_ID, inputs)
    return history_id, job_id


def find_dating_outputs(results_dir):
    """Locate the LSD2 dated-tree (NEXUS) and text report among a downloaded
    IQ-TREE dating job's outputs. Matches by normalized substring against the
    tool's declared output labels ('Tree labeled with dates', 'LSD Report')
    rather than exact filenames, which vary with Galaxy dataset numbering.
    Returns (timetree_path_or_None, report_path_or_None)."""
    try:
        names = os.listdir(results_dir)
    except OSError:
        return None, None
    timetree = report = None
    for name in names:
        n = re.sub(r'[^a-z0-9]', '', name.lower())
        path = os.path.join(results_dir, name)
        if not os.path.isfile(path):
            continue
        if 'labeledwithdates' in n or 'timetree' in n:
            timetree = path
        elif 'lsdreport' in n:
            report = path
    return timetree, report


def parse_timetree_newick(nexus_or_newick_path):
    """Read the LSD2 dated-tree output (NEXUS, node labels = ages) and return a
    plain Newick string with time-scaled branch lengths and tip labels intact,
    but internal node labels stripped. Keeping them would make the shared tree
    renderer misread ages as bootstrap-support percentages and color-code them
    accordingly, which would be actively misleading."""
    from Bio import Phylo
    import io
    try:
        tree = Phylo.read(nexus_or_newick_path, 'nexus')
    except Exception:
        tree = Phylo.read(nexus_or_newick_path, 'newick')
    for clade in tree.find_clades():
        clade.comment = None
        if not clade.is_terminal():
            clade.name = None
            clade.confidence = None
    buf = io.StringIO()
    Phylo.write(tree, buf, 'newick')
    return buf.getvalue().strip()
