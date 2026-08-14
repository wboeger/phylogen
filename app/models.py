import uuid
from datetime import datetime, timezone

from app import db


def _uuid():
    return uuid.uuid4().hex


class Job(db.Model):
    """One phylogenetic tree run. Flat — no project/account grouping (see plan).

    Never stores NCBI/Galaxy credentials — those live only in the in-process
    secret cache (app/secrets.py) for the lifetime of the pipeline calls that
    need them.
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
    marker = db.Column(db.String(10), nullable=False)        # preset code or 'custom'
    gene_query = db.Column(db.Text, nullable=False)          # resolved NCBI query actually used
    min_length = db.Column(db.Integer, default=400)
    outgroup_text = db.Column(db.Text)                       # optional, one name per line

    # Status
    status = db.Column(db.String(20), default='created')
    status_message = db.Column(db.Text)
    error_message = db.Column(db.Text)

    # Discovery / fetch results
    species_found = db.Column(db.JSON)    # list[str] — species actually recovered
    species_missing = db.Column(db.JSON)  # list[str] — requested but not found (manual mode)
    n_sequences = db.Column(db.Integer)

    # Files (paths under JOBS_DIR/<id>/)
    raw_fasta_path = db.Column(db.String(500))
    aligned_fasta_path = db.Column(db.String(500))
    trimmed_fasta_path = db.Column(db.String(500))

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
    galaxy_history_id = db.Column(db.String(64))
    galaxy_job_id = db.Column(db.String(64))

    @property
    def result_dir(self):
        from flask import current_app
        import os
        d = os.path.join(current_app.config['JOBS_DIR'], self.id)
        os.makedirs(d, exist_ok=True)
        return d

    def to_status_dict(self):
        return {
            'id': self.id,
            'status': self.status,
            'status_message': self.status_message,
            'error_message': self.error_message,
            'n_sequences': self.n_sequences,
            'species_found': self.species_found,
            'species_missing': self.species_missing,
            'nj_newick': self.nj_rooted_newick or self.nj_newick,
            'nj_root_note': self.nj_root_note,
            'ml_status': self.ml_status,
            'ml_message': self.ml_message,
            'ml_newick': self.ml_rooted_newick or self.ml_newick,
            'ml_has_support': self.ml_has_support,
            'ml_root_note': self.ml_root_note,
        }
