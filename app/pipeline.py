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


def submit_raxml(fasta_path, api_key, n_bootstraps=1000, model='GTR+G'):
    """Upload the trimmed alignment and run RAxML-NG (--all: ML search +
    bootstrap + support) on Galaxy. Returns (history_id, galaxy_job_id)."""
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
        'general_opts|cmdtype|model|model_type': 'single_string',
        'general_opts|cmdtype|model|model_string': model or 'GTR+G',
        'bootstrap_opts|bs_reps': int(n_bootstraps),
        'bootstrap_opts|bs_mre': 'true',
        'random_seed': 1234567890,
    }
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
