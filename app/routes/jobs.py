import io
import json
import os
import re
import shutil
import threading
import time
import zipfile
from datetime import datetime, timezone

from flask import (Blueprint, abort, current_app, g, jsonify, redirect,
                    render_template, request, send_file, url_for)

from app import db
from app import creds as secret_cache
from app import pipeline
from app.models import Job

jobs_bp = Blueprint('jobs', __name__)


# ── Helpers ───────────────────────────────────────────────────────────────────

def _lines(text):
    return [ln.strip() for ln in (text or '').splitlines() if ln.strip()]


def _fail(job_id, message):
    job = db.session.get(Job, job_id)
    if not job:
        return
    job.status = 'error'
    job.error_message = message
    job.status_message = message
    db.session.commit()


def _owns(job):
    return job.session_id == g.session_id


def _parse_fragments(form):
    """Build the job's fragment list from checked presets plus any repeated
    custom-fragment rows. Returns (fragments, errors)."""
    fragments, errors, seen_codes = [], [], set()

    for code in form.getlist('fragment_codes'):
        if code in pipeline.MARKER_PRESETS and code not in seen_codes:
            info = pipeline.MARKER_PRESETS[code]
            fragments.append({'code': code, 'label': f'{code} — {info["label"]}',
                               'query': info['query'], 'min_length': info['min_length']})
            seen_codes.add(code)

    names = form.getlist('custom_fragment_name')
    queries = form.getlist('custom_fragment_query')
    minlens = form.getlist('custom_fragment_minlen')
    for i, (name, query) in enumerate(zip(names, queries)):
        name, query = name.strip(), query.strip()
        if not name and not query:
            continue
        if not name or not query:
            errors.append(f'Custom fragment #{i + 1}: give it both a name and a search query.')
            continue
        code = re.sub(r'[^A-Za-z0-9_]', '', name.replace(' ', '_'))[:20] or f'FRAG{i + 1}'
        while code in seen_codes:
            code = f'{code}_'
        seen_codes.add(code)
        try:
            min_len = int(minlens[i]) if i < len(minlens) and minlens[i] else 300
        except ValueError:
            min_len = 300
        fragments.append({'code': code, 'label': name, 'query': query, 'min_length': min_len})

    if not fragments:
        errors.append('Select at least one marker, or add a custom fragment.')
    return fragments, errors


# ── Create ───────────────────────────────────────────────────────────────────

@jobs_bp.route('/jobs', methods=['POST'])
def create_job():
    f = request.form

    ncbi_email = (f.get('ncbi_email') or '').strip()
    ncbi_api_key = (f.get('ncbi_api_key') or '').strip()
    galaxy_api_key = (f.get('galaxy_api_key') or '').strip()
    input_mode = f.get('input_mode', 'manual')
    outgroup_text = f.get('outgroup_text', '')

    errors = []
    if '@' not in ncbi_email:
        errors.append('A valid NCBI contact email is required (Entrez requires one).')
    if not galaxy_api_key:
        errors.append('A usegalaxy.eu API key is required — align/trim runs on Galaxy.')

    if input_mode == 'manual':
        species = _lines(f.get('species_text', ''))
        if len(species) < 2:
            errors.append('List at least 2 species (one per line) to build a tree.')
        taxon_query = None
        species_text = '\n'.join(species)
    elif input_mode == 'taxon':
        taxon_query = (f.get('taxon_query') or '').strip()
        if not taxon_query:
            errors.append('Enter a higher taxon name for taxon discovery.')
        species_text = None
    else:
        errors.append('Unknown input mode.')
        taxon_query = None
        species_text = None

    fragments, frag_errors = _parse_fragments(f)
    errors.extend(frag_errors)

    try:
        taxon_max_species = int(f.get('taxon_max_species') or 40)
    except ValueError:
        taxon_max_species = 40
    taxon_max_species = max(3, min(taxon_max_species, 150))

    if errors:
        jobs = (Job.query.filter_by(session_id=g.session_id)
                .order_by(Job.created_at.desc()).limit(20).all())
        return render_template('index.html', jobs=jobs, markers=pipeline.MARKER_PRESETS,
                                errors=errors, form=f), 400

    job = Job(
        session_id=g.session_id,
        input_mode=input_mode,
        species_text=species_text,
        taxon_query=taxon_query,
        taxon_max_species=taxon_max_species,
        fragments=fragments,
        outgroup_text='\n'.join(_lines(outgroup_text)) or None,
        status='created',
        status_message='Job created — starting…',
    )
    db.session.add(job)
    db.session.commit()

    secret_cache.put(job.id, ncbi_email=ncbi_email, ncbi_api_key=ncbi_api_key or None,
                      galaxy_api_key=galaxy_api_key)

    app_obj = current_app._get_current_object()
    threading.Thread(target=_fetch_thread, args=(app_obj, job.id), daemon=True).start()

    return redirect(url_for('jobs.job_detail', job_id=job.id))


# ── Stage 1: fetch (pauses at 'fetched' for user review) ────────────────────

def _fetch_thread(app, job_id):
    with app.app_context():
        job = db.session.get(Job, job_id)
        if not job:
            return
        creds = secret_cache.get(job_id)
        if not creds:
            _fail(job_id, 'Your NCBI/Galaxy credentials were not available when the '
                           'pipeline started (they expire after inactivity) — please retry.')
            return
        try:
            _run_fetch(job, creds)
        except Exception as exc:
            _fail(job_id, str(exc))


def _run_fetch(job, creds):
    email = creds['ncbi_email']
    ncbi_key = creds.get('ncbi_api_key')

    def progress(msg):
        job.status_message = msg
        db.session.commit()

    job.status = 'fetching'
    db.session.commit()

    outgroup_names = _lines(job.outgroup_text) if job.outgroup_text else []
    fetch_results = {}
    fragment_files = {}
    all_species_seen = set()

    for frag in job.fragments:
        code, query, min_len = frag['code'], frag['query'], frag['min_length']
        progress(f'[{code}] Searching NCBI…')

        if job.input_mode == 'manual':
            species = _lines(job.species_text)
            records, _found, missing = pipeline.fetch_species_list(
                species, query, email, ncbi_key, min_len, progress)
        else:
            records, _found, n_total = pipeline.fetch_taxon(
                job.taxon_query, query, email, ncbi_key, min_len, job.taxon_max_species, progress)
            missing = []
            if n_total > job.taxon_max_species:
                progress(f'[{code}] {n_total} species had hits; capped to {job.taxon_max_species}.')

        og_records, og_missing = [], []
        if outgroup_names:
            og_records, _f, og_missing = pipeline.fetch_species_list(
                outgroup_names, query, email, ncbi_key, min_len, progress)

        combined = list(records) + list(og_records)
        fetch_results[code] = {
            'label': frag['label'],
            'found': pipeline.summarize_fetch(combined),
            'missing': missing,
            'outgroup_missing': og_missing,
            'n_sequences': len(combined),
        }

        raw_path = os.path.join(job.result_dir, f'{code}_raw.fasta')
        pipeline.write_fasta(combined, raw_path)
        fragment_files[code] = {'raw': raw_path}
        for r in combined:
            all_species_seen.add(r.id.split('|', 1)[1] if '|' in r.id else r.id)

    job.fetch_results = fetch_results
    job.fragment_files = fragment_files
    if not all_species_seen:
        raise RuntimeError('No sequences were found for any fragment/species combination. '
                            'Try different markers, or broaden the species/taxon.')
    job.n_sequences = sum(v['n_sequences'] for v in fetch_results.values())
    job.status = 'fetched'
    job.status_message = (
        f'Fetched {len(job.fragments)} fragment(s) across {len(all_species_seen)} species. '
        f'Review the sequences below, then approve to align.')
    db.session.commit()


# ── Stage 2: user approves -> align + trim + concatenate + NJ ───────────────

def _all_found_accessions(job):
    accs = set()
    for frag_result in (job.fetch_results or {}).values():
        for row in frag_result.get('found', []):
            accs.add(row['accession'])
    return accs


@jobs_bp.route('/api/job/<job_id>/approve_and_align', methods=['POST'])
def approve_and_align(job_id):
    job = db.session.get(Job, job_id) or abort(404)
    if job.status != 'fetched':
        return jsonify({'error': f'Job is not awaiting approval (status: {job.status}).'}), 400

    rejected_raw = request.form.get('rejected_accessions', '')
    rejected = sorted({a.strip() for a in rejected_raw.split(',') if a.strip()})
    all_found = _all_found_accessions(job)
    approved_count = len(all_found - set(rejected))

    if all_found and approved_count == 0:
        # Every fetched sequence was rejected -- nothing to align. Abandon the
        # job rather than run a pipeline on zero sequences.
        result_dir = job.result_dir
        db.session.delete(job)
        db.session.commit()
        secret_cache.forget(job_id)
        shutil.rmtree(result_dir, ignore_errors=True)
        return jsonify({'status': 'rejected_all', 'redirect': url_for('main.index')})

    creds = secret_cache.get(job_id) or {}
    galaxy_api_key = (request.form.get('galaxy_api_key') or '').strip() or creds.get('galaxy_api_key')
    if not galaxy_api_key:
        return jsonify({'error': 'A usegalaxy.eu API key is required to align/trim.'}), 400
    secret_cache.put(job_id, **{**creds, 'galaxy_api_key': galaxy_api_key})

    job.rejected_accessions = rejected
    job.status = 'aligning'
    job.status_message = ('Approved. Aligning fragment(s) on Galaxy…' if not rejected else
                           f'Approved ({len(rejected)} sequence(s) excluded). Aligning fragment(s) on Galaxy…')
    db.session.commit()

    app_obj = current_app._get_current_object()
    threading.Thread(target=_align_thread, args=(app_obj, job_id), daemon=True).start()
    return jsonify({'status': 'aligning'})


def _align_thread(app, job_id):
    with app.app_context():
        job = db.session.get(Job, job_id)
        if not job:
            return
        creds = secret_cache.get(job_id) or {}
        galaxy_key = creds.get('galaxy_api_key')
        if not galaxy_key:
            _fail(job_id, 'Galaxy API key is no longer available — please retry from the review step.')
            return
        try:
            _run_align_concat_nj(job, galaxy_key)
        except Exception as exc:
            _fail(job_id, str(exc))


def _run_align_concat_nj(job, galaxy_key):
    def progress(msg):
        job.status_message = msg
        db.session.commit()

    fragment_files = dict(job.fragment_files or {})
    rejected = job.rejected_accessions or []
    frag_paths = []
    for frag in job.fragments:
        code = frag['code']
        raw_path = fragment_files.get(code, {}).get('raw')
        if raw_path and rejected:
            approved_path = os.path.join(job.result_dir, f'{code}_raw_approved.fasta')
            removed = pipeline.filter_fasta_exclude(raw_path, rejected, approved_path)
            if removed:
                progress(f'[{code}] Excluded {len(removed)} rejected sequence(s).')
            raw_path = approved_path
        if not raw_path or pipeline.count_fasta(raw_path) < 2:
            progress(f'[{code}] Skipped — fewer than 2 approved sequences for this fragment.')
            continue

        progress(f'[{code}] Aligning on Galaxy (MAFFT)…')
        aligned_path = os.path.join(job.result_dir, f'{code}_aligned.fasta')
        trimmed_path = os.path.join(job.result_dir, f'{code}_trimmed.fasta')
        _, _, flipped, trim_skipped = pipeline.galaxy_align_trim(
            galaxy_key, raw_path, aligned_path, trimmed_path)

        model = pipeline.estimate_model(trimmed_path)
        n_seq = pipeline.count_fasta(trimmed_path)
        fragment_files[code] = {**fragment_files.get(code, {}), 'aligned': aligned_path,
                                 'trimmed': trimmed_path, 'model': model, 'n_sequences': n_seq}
        frag_paths.append((trimmed_path, code))

        note = ''
        if flipped:
            note += f' {len(flipped)} sequence(s) reverse-complemented.'
        if trim_skipped:
            note += ' trimAl produced no usable output; using the untrimmed alignment.'
        progress(f'[{code}] Aligned & trimmed ({n_seq} sequences). Estimated model: {model}.{note}')

    job.fragment_files = fragment_files
    db.session.commit()

    if not frag_paths:
        raise RuntimeError('No fragment had enough sequences to align. '
                            'Try different markers, or broaden the species/taxon.')

    concat_path = os.path.join(job.result_dir, 'concat_trimmed.fasta')
    n_taxa, partition_spec, coverage = pipeline.concatenate_fragments(frag_paths, concat_path)
    job.concat_path = concat_path
    job.partition_spec = partition_spec
    job.species_coverage = coverage
    job.n_sequences = n_taxa
    db.session.commit()

    job.status = 'nj_running'
    progress(f'Concatenated {len(frag_paths)} fragment(s) across {n_taxa} taxa. Computing NJ tree…')
    try:
        models = {code: v.get('model') for code, v in fragment_files.items()}
        nj = pipeline.local_nj_tree_partitioned(concat_path, partition_spec, models)
        job.nj_newick = nj
        if job.outgroup_text:
            rooted, rooted_on, not_found = pipeline.reroot_newick(nj, _lines(job.outgroup_text))
            job.nj_rooted_newick = rooted
            job.nj_root_note = f'Rooted on {rooted_on}.' + (
                f' Not found in tree: {", ".join(not_found)}.' if not_found else '')
        job.status = 'nj_ready'
        job.status_message = (
            f'Tree ready ({n_taxa} taxa, {len(frag_paths)} fragment(s)). Review the sequences '
            f'included below, then confirm to enable the bootstrapped ML tree build.')
    except Exception as exc:
        job.status = 'trimmed'
        job.status_message = f'Alignment/trim/concatenation complete, but the NJ preview failed: {exc}'
    db.session.commit()


@jobs_bp.route('/api/job/<job_id>/confirm_nj', methods=['POST'])
def confirm_nj(job_id):
    job = db.session.get(Job, job_id) or abort(404)
    if job.status not in ('nj_ready', 'trimmed'):
        return jsonify({'error': 'Nothing to confirm yet.'}), 400
    job.nj_confirmed = True
    db.session.commit()
    return jsonify({'status': 'ok'})


# ── ML tree (Galaxy RAxML-NG) ────────────────────────────────────────────────

@jobs_bp.route('/api/job/<job_id>/build_ml', methods=['POST'])
def build_ml(job_id):
    job = db.session.get(Job, job_id) or abort(404)
    if not job.nj_confirmed:
        return jsonify({'error': 'Confirm the sequences included in the NJ tree first.'}), 400
    if not job.concat_path or not os.path.exists(job.concat_path):
        return jsonify({'error': 'No concatenated alignment available yet.'}), 400
    galaxy_api_key = (request.form.get('galaxy_api_key') or '').strip()
    if not galaxy_api_key:
        cached = secret_cache.get(job_id)
        galaxy_api_key = (cached or {}).get('galaxy_api_key', '')
    if not galaxy_api_key:
        return jsonify({'error': 'A usegalaxy.eu API key is required.'}), 400
    if job.ml_status == 'running':
        return jsonify({'error': 'An ML tree run is already in progress for this job.'}), 400

    job.ml_status = 'running'
    job.ml_message = 'Submitting alignment to Galaxy RAxML-NG…'
    job.ml_started_at = datetime.now(timezone.utc)
    db.session.commit()

    app_obj = current_app._get_current_object()
    threading.Thread(target=_raxml_thread, args=(app_obj, job_id, galaxy_api_key),
                      daemon=True).start()
    return jsonify({'status': 'running', 'message': job.ml_message})


def _raxml_thread(app, job_id, galaxy_key):
    with app.app_context():
        job = db.session.get(Job, job_id)
        if not job:
            return
        try:
            fragment_files = job.fragment_files or {}
            partition_models = {code: v.get('model', 'GTR+G') for code, v in fragment_files.items()}
            hist, gjob = pipeline.submit_raxml(
                job.concat_path, galaxy_key, partition_spec=job.partition_spec,
                partition_models=partition_models)
            job.galaxy_history_id = hist
            job.galaxy_job_id = gjob
            job.ml_message = 'RAxML-NG running on usegalaxy.eu (ML search + bootstrap)…'
            db.session.commit()

            deadline = time.time() + 3 * 3600
            stage = 'RUNNING'
            while time.time() < deadline:
                stage, msg = pipeline.galaxy_check_status(galaxy_key, gjob)
                job.ml_message = f'{stage}: {(msg or "")[:300]}'
                db.session.commit()
                if stage == 'COMPLETED':
                    break
                if stage in ('FAILED', 'SUSPENDED'):
                    raise RuntimeError(f'Galaxy RAxML-NG job {stage.lower()}: {(msg or "")[:400]}')
                time.sleep(15)
            else:
                raise RuntimeError('Timed out waiting for Galaxy RAxML-NG (3h).')

            dest_dir = os.path.join(job.result_dir, 'raxml')
            pipeline.galaxy_download_results(galaxy_key, gjob, dest_dir)
            tree_path, has_support = pipeline.find_best_tree(dest_dir)
            if not tree_path:
                raise RuntimeError('RAxML-NG finished but produced no readable tree file.')
            with open(tree_path) as fh:
                newick = fh.read().strip()
            job.ml_newick = newick
            job.ml_has_support = has_support

            if job.outgroup_text:
                rooted, rooted_on, not_found = pipeline.reroot_newick(newick, _lines(job.outgroup_text))
                job.ml_rooted_newick = rooted
                job.ml_root_note = f'Rooted on {rooted_on}.' + (
                    f' Not found in tree: {", ".join(not_found)}.' if not_found else '')

            job.ml_status = 'ready'
            job.ml_message = ('ML tree ready.' if has_support else
                               'ML tree ready (no bootstrap-support values were recovered).')
        except Exception as exc:
            job.ml_status = 'error'
            job.ml_message = str(exc)
        db.session.commit()


# ── Divergence-time dating (IQ-TREE + LSD2 on Galaxy) ────────────────────────

@jobs_bp.route('/api/job/<job_id>/build_dating', methods=['POST'])
def build_dating(job_id):
    job = db.session.get(Job, job_id) or abort(404)
    if not job.ml_newick or job.ml_status != 'ready':
        return jsonify({'error': 'Build the ML tree first.'}), 400
    if job.dating_status == 'running':
        return jsonify({'error': 'A dating run is already in progress for this job.'}), 400

    galaxy_api_key = (request.form.get('galaxy_api_key') or '').strip()
    if not galaxy_api_key:
        cached = secret_cache.get(job_id)
        galaxy_api_key = (cached or {}).get('galaxy_api_key', '')
    if not galaxy_api_key:
        return jsonify({'error': 'A usegalaxy.eu API key is required.'}), 400

    try:
        calibrations = json.loads(request.form.get('calibrations', '[]'))
    except (TypeError, ValueError):
        return jsonify({'error': 'Malformed calibration data.'}), 400
    if not isinstance(calibrations, list) or not calibrations:
        return jsonify({'error': 'Add at least one calibration point.'}), 400

    try:
        date_ci = max(0, int(request.form.get('date_ci') or 100))
    except ValueError:
        date_ci = 100
    try:
        clock_sd = float(request.form.get('clock_sd') or 0.2)
    except ValueError:
        clock_sd = 0.2

    base_newick = job.ml_rooted_newick or job.ml_newick
    tip_names = pipeline.newick_tip_names(base_newick)
    date_path = os.path.join(job.result_dir, 'dating_calibrations.txt')
    _, matched_rows, unmatched = pipeline.build_calibration_file(calibrations, tip_names, date_path)
    if not matched_rows:
        detail = f' Not found: {", ".join(unmatched)}.' if unmatched else ''
        return jsonify({'error': f'None of the calibration points matched a tip in the tree.{detail}'}), 400

    job.calibrations = calibrations
    job.dating_status = 'running'
    warn = f' ({len(unmatched)} name(s) not matched: {", ".join(unmatched)})' if unmatched else ''
    job.dating_message = f'Submitting to Galaxy (IQ-TREE / LSD2)…{warn}'
    job.dating_started_at = datetime.now(timezone.utc)
    db.session.commit()

    app_obj = current_app._get_current_object()
    threading.Thread(target=_dating_thread,
                      args=(app_obj, job_id, galaxy_api_key, date_path, date_ci, clock_sd),
                      daemon=True).start()
    return jsonify({'status': 'running', 'message': job.dating_message})


def _dating_thread(app, job_id, galaxy_key, date_path, date_ci, clock_sd):
    with app.app_context():
        job = db.session.get(Job, job_id)
        if not job:
            return
        try:
            base_newick = job.ml_rooted_newick or job.ml_newick
            fragment_files = job.fragment_files or {}
            model = 'GTR+G'
            if job.fragments:
                model = fragment_files.get(job.fragments[0]['code'], {}).get('model') or 'GTR+G'
            outgroup_names = _lines(job.outgroup_text) if job.outgroup_text else None

            hist, gjob = pipeline.submit_dating(
                job.concat_path, base_newick, date_path, galaxy_key, model=model,
                outgroup_names=outgroup_names, date_ci=date_ci, clock_sd=clock_sd)
            job.galaxy_dating_history_id = hist
            job.galaxy_dating_job_id = gjob
            job.dating_message = 'IQ-TREE / LSD2 running on usegalaxy.eu…'
            db.session.commit()

            deadline = time.time() + 3 * 3600
            stage = 'RUNNING'
            while time.time() < deadline:
                stage, msg = pipeline.galaxy_check_status(galaxy_key, gjob)
                job.dating_message = f'{stage}: {(msg or "")[:300]}'
                db.session.commit()
                if stage == 'COMPLETED':
                    break
                if stage in ('FAILED', 'SUSPENDED'):
                    raise RuntimeError(f'Galaxy dating job {stage.lower()}: {(msg or "")[:400]}')
                time.sleep(15)
            else:
                raise RuntimeError('Timed out waiting for Galaxy IQ-TREE/LSD2 (3h).')

            dest_dir = os.path.join(job.result_dir, 'dating')
            pipeline.galaxy_download_results(galaxy_key, gjob, dest_dir)
            timetree_path, report_path = pipeline.find_dating_outputs(dest_dir)
            if not timetree_path:
                raise RuntimeError('IQ-TREE/LSD2 finished but produced no readable dated-tree output.')
            job.dating_newick = pipeline.parse_timetree_newick(timetree_path)
            if report_path:
                with open(report_path, errors='ignore') as fh:
                    job.dating_report = fh.read()[:20000]
            job.dating_status = 'ready'
            job.dating_message = 'Divergence-time estimate ready.'
        except Exception as exc:
            job.dating_status = 'error'
            job.dating_message = str(exc)
        db.session.commit()


# ── Detail / status / reroot / download / delete ────────────────────────────

@jobs_bp.route('/job/<job_id>')
def job_detail(job_id):
    job = db.session.get(Job, job_id) or abort(404)
    return render_template('job.html', job=job, is_owner=_owns(job))


@jobs_bp.route('/api/job/<job_id>/status')
def job_status(job_id):
    job = db.session.get(Job, job_id) or abort(404)
    return jsonify(job.to_status_dict())


@jobs_bp.route('/api/job/<job_id>/reroot', methods=['POST'])
def reroot(job_id):
    job = db.session.get(Job, job_id) or abort(404)
    data = request.get_json(silent=True) or {}
    which = data.get('tree')
    names = [n.strip() for n in (data.get('names') or '').replace(',', '\n').splitlines() if n.strip()]
    base_newick = job.nj_newick if which == 'nj' else job.ml_newick
    if not base_newick:
        return jsonify({'error': 'No tree to reroot.'}), 400
    try:
        rooted, rooted_on, not_found = pipeline.reroot_newick(base_newick, names)
    except Exception as exc:
        return jsonify({'error': f'Rerooting failed: {exc}'}), 400
    note = f'Rooted on {rooted_on}.' + (f' Not found: {", ".join(not_found)}.' if not_found else '')
    if which == 'nj':
        job.nj_rooted_newick = rooted
        job.nj_root_note = note
    else:
        job.ml_rooted_newick = rooted
        job.ml_root_note = note
    db.session.commit()
    return jsonify({'newick': rooted, 'note': note})


@jobs_bp.route('/job/<job_id>/download.zip')
def download_zip(job_id):
    job = db.session.get(Job, job_id) or abort(404)
    fragment_files = job.fragment_files or {}
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, 'w', zipfile.ZIP_DEFLATED) as zf:
        readme = [
            f'phylogen job {job.id}',
            f'created: {job.created_at}',
            f'input mode: {job.input_mode}',
        ]
        if job.input_mode == 'manual':
            readme.append(f'species requested: {job.species_text}')
        else:
            readme.append(f'taxon: {job.taxon_query} (max {job.taxon_max_species} species)')
        if job.outgroup_text:
            readme.append(f'outgroup: {job.outgroup_text}')
        if job.rejected_accessions:
            readme.append(f'user-rejected accessions (excluded from alignment): {", ".join(job.rejected_accessions)}')
        readme.append('')
        readme.append('fragments:')
        for frag in job.fragments:
            code = frag['code']
            info = fragment_files.get(code, {})
            readme.append(f"  {code}: query=\"{frag['query']}\" min_length={frag['min_length']} "
                           f"model={info.get('model', 'n/a')} n_sequences={info.get('n_sequences', 'n/a')}")
            fr = (job.fetch_results or {}).get(code, {})
            if fr.get('missing'):
                readme.append(f"    not found: {', '.join(fr['missing'])}")
        if job.species_coverage:
            readme.append('')
            readme.append('species -> fragments included in the concatenated tree:')
            for sp, cov in sorted(job.species_coverage.items()):
                readme.append(f'  {sp}: {cov}')
        if job.calibrations:
            readme.append('')
            readme.append('divergence-dating calibration points (IQ-TREE/LSD2, millions of years before present):')
            for cal in job.calibrations:
                taxa = ', '.join(cal.get('taxa') or [])
                readme.append(f"  {taxa}: min={cal.get('min_age')} max={cal.get('max_age')}")
            if job.dating_status:
                readme.append(f'dating status: {job.dating_status} — {job.dating_message or ""}')
        zf.writestr('README.txt', '\n'.join(readme) + '\n')

        for code, info in fragment_files.items():
            for label, key in (('raw', 'raw'), ('aligned', 'aligned'), ('trimmed', 'trimmed')):
                path = info.get(key)
                if path and os.path.exists(path):
                    zf.write(path, f'{code}_{label}.fasta')
        if job.concat_path and os.path.exists(job.concat_path):
            zf.write(job.concat_path, 'concat_trimmed.fasta')
        if job.nj_newick:
            zf.writestr('nj_tree.nwk', job.nj_newick)
        if job.nj_rooted_newick:
            zf.writestr('nj_tree_rooted.nwk', job.nj_rooted_newick)
        if job.ml_newick:
            zf.writestr('ml_tree_raxmlng.nwk', job.ml_newick)
        if job.ml_rooted_newick:
            zf.writestr('ml_tree_raxmlng_rooted.nwk', job.ml_rooted_newick)
        if job.dating_newick:
            zf.writestr('dating_timetree.nwk', job.dating_newick)
        if job.dating_report:
            zf.writestr('dating_lsd2_report.txt', job.dating_report)
        cal_path = os.path.join(job.result_dir, 'dating_calibrations.txt')
        if os.path.exists(cal_path):
            zf.write(cal_path, 'dating_calibrations.txt')
    buf.seek(0)
    return send_file(buf, mimetype='application/zip', as_attachment=True,
                      download_name=f'phylogen_{job.id[:8]}.zip')


@jobs_bp.route('/api/job/<job_id>/delete', methods=['POST'])
def delete_job(job_id):
    job = db.session.get(Job, job_id) or abort(404)
    if not _owns(job):
        return jsonify({'error': 'Not your job.'}), 403
    result_dir = job.result_dir
    db.session.delete(job)
    db.session.commit()
    secret_cache.forget(job_id)
    shutil.rmtree(result_dir, ignore_errors=True)
    return jsonify({'status': 'ok'})
