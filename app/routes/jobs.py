import io
import os
import shutil
import threading
import time
import zipfile

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


# ── Create ───────────────────────────────────────────────────────────────────

@jobs_bp.route('/jobs', methods=['POST'])
def create_job():
    f = request.form

    ncbi_email = (f.get('ncbi_email') or '').strip()
    ncbi_api_key = (f.get('ncbi_api_key') or '').strip()
    galaxy_api_key = (f.get('galaxy_api_key') or '').strip()
    input_mode = f.get('input_mode', 'manual')
    marker = f.get('marker', 'COI')
    custom_query = (f.get('custom_query') or '').strip()
    min_length_raw = (f.get('min_length') or '').strip()
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

    if marker == 'custom':
        if not custom_query:
            errors.append('Enter a custom NCBI search query, or pick a preset marker.')
        gene_query = custom_query
        default_min_length = 400
    elif marker in pipeline.MARKER_PRESETS:
        gene_query = pipeline.MARKER_PRESETS[marker]['query']
        default_min_length = pipeline.MARKER_PRESETS[marker]['min_length']
    else:
        errors.append('Unknown marker.')
        gene_query = ''
        default_min_length = 400

    try:
        min_length = int(min_length_raw) if min_length_raw else default_min_length
    except ValueError:
        min_length = default_min_length

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
        marker=marker,
        gene_query=gene_query,
        min_length=min_length,
        outgroup_text='\n'.join(_lines(outgroup_text)) or None,
        status='created',
        status_message='Job created — starting…',
    )
    db.session.add(job)
    db.session.commit()

    secret_cache.put(job.id, ncbi_email=ncbi_email, ncbi_api_key=ncbi_api_key or None,
                      galaxy_api_key=galaxy_api_key)

    app_obj = current_app._get_current_object()
    threading.Thread(target=_pipeline_thread, args=(app_obj, job.id), daemon=True).start()

    return redirect(url_for('jobs.job_detail', job_id=job.id))


# ── Background pipeline ──────────────────────────────────────────────────────

def _pipeline_thread(app, job_id):
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
            _run_pipeline(job, creds)
        except Exception as exc:
            _fail(job_id, str(exc))


def _run_pipeline(job, creds):
    email = creds['ncbi_email']
    ncbi_key = creds.get('ncbi_api_key')
    galaxy_key = creds['galaxy_api_key']

    def progress(msg):
        job.status_message = msg
        db.session.commit()

    job.status = 'fetching'
    db.session.commit()

    if job.input_mode == 'manual':
        species = _lines(job.species_text)
        records, found, missing = pipeline.fetch_species_list(
            species, job.gene_query, email, ncbi_key, job.min_length, progress)
        job.species_found = found
        job.species_missing = missing
    else:
        records, found, n_total = pipeline.fetch_taxon(
            job.taxon_query, job.gene_query, email, ncbi_key, job.min_length,
            job.taxon_max_species, progress)
        job.species_found = found
        job.species_missing = []
        if n_total > job.taxon_max_species:
            progress(f'Found {n_total} species with hits; capped to {job.taxon_max_species}.')
    db.session.commit()

    outgroup_records = []
    if job.outgroup_text:
        names = _lines(job.outgroup_text)
        outgroup_records, _og_found, og_missing = pipeline.fetch_species_list(
            names, job.gene_query, email, ncbi_key, job.min_length, progress)
        if og_missing:
            progress(f'No sequence found for outgroup(s): {", ".join(og_missing)}.')

    all_records = list(records) + list(outgroup_records)
    if len(all_records) < 3:
        raise RuntimeError(
            f'Only {len(all_records)} sequence(s) recovered — at least 3 are needed to '
            f'build a tree. Try a different marker, broaden the species/taxon, or check spelling.')

    raw_path = os.path.join(job.result_dir, 'raw.fasta')
    pipeline.write_fasta(all_records, raw_path)
    job.raw_fasta_path = raw_path
    job.n_sequences = len(all_records)
    db.session.commit()

    job.status = 'aligning'
    progress(f'{len(all_records)} sequences fetched. Aligning on Galaxy (MAFFT)…')

    aligned_path = os.path.join(job.result_dir, 'aligned.fasta')
    trimmed_path = os.path.join(job.result_dir, 'trimmed.fasta')
    _, _, flipped, trim_skipped = pipeline.galaxy_align_trim(
        galaxy_key, raw_path, aligned_path, trimmed_path)
    job.aligned_fasta_path = aligned_path
    job.trimmed_fasta_path = trimmed_path
    db.session.commit()

    note = ''
    if flipped:
        note += f' {len(flipped)} sequence(s) reverse-complemented to match orientation.'
    if trim_skipped:
        note += ' trimAl produced no usable output — using the untrimmed alignment.'

    job.status = 'nj_running'
    progress('Computing neighbor-joining preview tree…' + note)
    try:
        nj = pipeline.local_nj_tree(trimmed_path)
        job.nj_newick = nj
        if job.outgroup_text:
            rooted, rooted_on, not_found = pipeline.reroot_newick(nj, _lines(job.outgroup_text))
            job.nj_rooted_newick = rooted
            job.nj_root_note = f'Rooted on {rooted_on}.' + (
                f' Not found in tree: {", ".join(not_found)}.' if not_found else '')
        job.status = 'nj_ready'
        job.status_message = (f'Tree ready ({job.n_sequences} sequences).{note} '
                               f'Download the results, or build a bootstrapped ML tree below.')
    except Exception as exc:
        job.status = 'trimmed'
        job.status_message = f'Alignment/trim complete, but the NJ preview failed: {exc}.{note}'
    db.session.commit()


# ── ML tree (Galaxy RAxML-NG) ────────────────────────────────────────────────

@jobs_bp.route('/api/job/<job_id>/build_ml', methods=['POST'])
def build_ml(job_id):
    job = db.session.get(Job, job_id) or abort(404)
    if not job.trimmed_fasta_path or not os.path.exists(job.trimmed_fasta_path):
        return jsonify({'error': 'No trimmed alignment available yet.'}), 400
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
            hist, gjob = pipeline.submit_raxml(job.trimmed_fasta_path, galaxy_key)
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
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, 'w', zipfile.ZIP_DEFLATED) as zf:
        readme = [
            f'phylogen job {job.id}',
            f'created: {job.created_at}',
            f'input mode: {job.input_mode}',
            f'marker: {job.marker}',
            f'NCBI query: {job.gene_query}',
            f'min length: {job.min_length}',
        ]
        if job.input_mode == 'manual':
            readme.append(f'species requested: {job.species_text}')
        else:
            readme.append(f'taxon: {job.taxon_query} (max {job.taxon_max_species} species)')
        if job.species_found:
            readme.append(f'species recovered: {", ".join(job.species_found)}')
        if job.species_missing:
            readme.append(f'species NOT found: {", ".join(job.species_missing)}')
        if job.outgroup_text:
            readme.append(f'outgroup: {job.outgroup_text}')
        zf.writestr('README.txt', '\n'.join(readme) + '\n')

        for label, path in (('raw.fasta', job.raw_fasta_path),
                             ('aligned.fasta', job.aligned_fasta_path),
                             ('trimmed.fasta', job.trimmed_fasta_path)):
            if path and os.path.exists(path):
                zf.write(path, label)
        if job.nj_newick:
            zf.writestr('nj_tree.nwk', job.nj_newick)
        if job.nj_rooted_newick:
            zf.writestr('nj_tree_rooted.nwk', job.nj_rooted_newick)
        if job.ml_newick:
            zf.writestr('ml_tree_raxmlng.nwk', job.ml_newick)
        if job.ml_rooted_newick:
            zf.writestr('ml_tree_raxmlng_rooted.nwk', job.ml_rooted_newick)
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
