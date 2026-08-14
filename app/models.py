import uuid
from datetime import datetime, timezone

from app import db


def _uuid():
    return uuid.uuid4().hex


def _iso_utc(dt):
    """isoformat() with an explicit UTC marker. SQLite round-trips datetimes as
    naive (tzinfo dropped even though the value is always UTC — see
    datetime.now(timezone.utc) below), and a bare 'no-offset' ISO string gets
    parsed as *local* time by JS `new Date()` — silently corrupting every
    client-side elapsed-time calculation. Always emit an explicit offset."""
    if dt is None:
        return None
    if dt.tzinfo is not None:
        dt = dt.astimezone(timezone.utc).replace(tzinfo=None)
    return dt.isoformat() + 'Z'


class Job(db.Model):
    """One phylogenetic tree run. Flat — no project/account grouping (see plan).

    Never stores NCBI/Galaxy credentials — those live only in the in-process
    secret cache (app/creds.py) for the lifetime of the pipeline calls that
    need them.

    Pipeline has two user confirmation gates:
      1. after fetch ('fetched' status) — review found/missing sequences per
         fragment before spending Galaxy time aligning them.
      2. after the NJ preview tree ('nj_ready', nj_confirmed flag) — review
         which species/fragments actually made it into the tree before
         committing to a (potentially hours-long) RAxML-NG bootstrap run.
    """
    __tablename__ = 'jobs'

    id = db.Column(db.String(32), primary_key=True, default=_uuid)
    session_id = db.Column(db.String(64), index=True, nullable=False)
    created_at = db.Column(db.DateTime, default=lambda: datetime.now(timezone.utc))
    updated_at = db.Column(db.DateTime, default=lambda: datetime.now(timezone.utc),
                            onupdate=lambda: datetime.now(timezone.utc))

    # Input
    input_mode = db.Column(db.String(10), nullable=False)   # 'manual' | 'taxon'
    species_text = db.Column(db.Text)                        # manual mode: one name per line
    taxon_query = db.Column(db.String(200))                  # taxon mode: higher-taxon name
    taxon_max_species = db.Column(db.Integer, default=40)
    outgroup_text = db.Column(db.Text)                       # optional, one name per line

    # Fragments/markers: list of {code, label, query, min_length}. One or more —
    # multiple fragments are fetched independently per species, then aligned,
    # trimmed, and concatenated (partitioned) for the tree.
    fragments = db.Column(db.JSON, nullable=False)

    # Status
    status = db.Column(db.String(20), default='created')
    status_message = db.Column(db.Text)
    error_message = db.Column(db.Text)
    nj_confirmed = db.Column(db.Boolean, default=False)

    # Per-fragment fetch review: {code: {found:[{species,accession,length}],
    # missing:[...], outgroup_missing:[...], n_sequences}}
    fetch_results = db.Column(db.JSON)
    n_sequences = db.Column(db.Integer)

    # Per-fragment files + estimated model + post-align sequence count:
    # {code: {raw, aligned, trimmed, model, n_sequences}}
    fragment_files = db.Column(db.JSON)

    # Concatenated (partitioned) alignment
    concat_path = db.Column(db.String(500))
    partition_spec = db.Column(db.JSON)     # [{name, start, end}] 1-based inclusive
    species_coverage = db.Column(db.JSON)   # {'genus species': '18S+COI', ...}

    # Trees
    nj_newick = db.Column(db.Text)
    nj_rooted_newick = db.Column(db.Text)
    nj_root_note = db.Column(db.Text)
    ml_newick = db.Column(db.Text)
    ml_rooted_newick = db.Column(db.Text)
    ml_root_note = db.Column(db.Text)
    ml_has_support = db.Column(db.Boolean, default=False)

    # Galaxy RAxML-NG job tracking
    ml_status = db.Column(db.String(20))  # None | running | ready | error
    ml_message = db.Column(db.Text)
    ml_started_at = db.Column(db.DateTime)
    galaxy_history_id = db.Column(db.String(64))
    galaxy_job_id = db.Column(db.String(64))

    @property
    def created_at_iso(self):
        return _iso_utc(self.created_at)

    @property
    def result_dir(self):
        from flask import current_app
        import os
        d = os.path.join(current_app.config['JOBS_DIR'], self.id)
        os.makedirs(d, exist_ok=True)
        return d

    # Four-stage pipeline: 0 fetch, 1 review (user gate), 2 align+trim+concat,
    # 3 build NJ tree. Drives the visual stage bar on the job page — computed
    # from status plus which files exist, so even an error mid-pipeline shows
    # how far it got.
    def stage_progress(self):
        if self.status in ('created', 'fetching'):
            return 0, 'active'
        if self.status == 'fetched':
            return 1, 'waiting'
        if self.status == 'aligning':
            return 2, 'active'
        if self.status == 'nj_running':
            return 3, 'active'
        if self.status == 'nj_ready':
            return 3, 'done'
        if self.status == 'trimmed':
            return 3, 'warn'          # aligned/trimmed/concat fine, NJ preview failed
        if self.status == 'error':
            if self.concat_path:
                return 3, 'error'
            if self.fragment_files and any(v.get('trimmed') for v in self.fragment_files.values()):
                return 2, 'error'
            if self.fetch_results:
                return 1, 'error'
            return 0, 'error'
        return 0, 'active'

    def to_status_dict(self):
        stage_index, stage_state = self.stage_progress()
        fragment_files = self.fragment_files or {}
        fragment_models = {code: v.get('model') for code, v in fragment_files.items() if v.get('model')}
        fragment_counts = {code: v.get('n_sequences') for code, v in fragment_files.items() if v.get('n_sequences')}
        return {
            'id': self.id,
            'status': self.status,
            'status_message': self.status_message,
            'error_message': self.error_message,
            'n_sequences': self.n_sequences,
            'fetch_results': self.fetch_results,
            'fragment_models': fragment_models,
            'fragment_counts': fragment_counts,
            'species_coverage': self.species_coverage,
            'nj_confirmed': self.nj_confirmed,
            'nj_newick': self.nj_rooted_newick or self.nj_newick,
            'nj_root_note': self.nj_root_note,
            'ml_status': self.ml_status,
            'ml_message': self.ml_message,
            'ml_newick': self.ml_rooted_newick or self.ml_newick,
            'ml_has_support': self.ml_has_support,
            'ml_root_note': self.ml_root_note,
            'created_at': self.created_at_iso,
            'ml_started_at': _iso_utc(self.ml_started_at),
            'stage_index': stage_index,
            'stage_state': stage_state,
        }
