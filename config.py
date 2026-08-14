import os

BASE_DIR = os.path.abspath(os.path.dirname(__file__))

# Ephemeral storage only (no Railway Volume — see DEPLOY.md). A redeploy wipes
# this; that's fine, since every finished job's real deliverable is the ZIP the
# user downloads, and the homepage job list is a session-cookie convenience,
# not a system of record.
DATA_DIR = os.environ.get('DATA_DIR') or os.path.join(BASE_DIR, 'data')
os.makedirs(DATA_DIR, exist_ok=True)

JOBS_DIR = os.path.join(DATA_DIR, 'jobs')
os.makedirs(JOBS_DIR, exist_ok=True)


class Config:
    SECRET_KEY = os.environ.get('SECRET_KEY', 'dev-key-change-in-production')
    SQLALCHEMY_DATABASE_URI = f"sqlite:///{os.path.join(DATA_DIR, 'db.sqlite')}"
    SQLALCHEMY_TRACK_MODIFICATIONS = False
    SQLALCHEMY_ENGINE_OPTIONS = {'connect_args': {'timeout': 30}}
    JOBS_DIR = JOBS_DIR
    MAX_CONTENT_LENGTH = 2 * 1024 * 1024  # 2 MB max request body (species/taxon text + keys)

    # Galaxy (usegalaxy.eu) tool IDs for MAFFT / trimAl / RAxML-NG. Every job
    # supplies its OWN Galaxy API key at request time — never a server secret —
    # but the tool wiring (which Galaxy tool, which input/param names) is fixed
    # infrastructure, override-able via env if usegalaxy.eu updates versions.
    GALAXY_BASE_URL = os.environ.get('GALAXY_BASE_URL', 'https://usegalaxy.eu')
    GALAXY_MAFFT_TOOL_ID = os.environ.get(
        'GALAXY_MAFFT_TOOL_ID',
        'toolshed.g2.bx.psu.edu/repos/rnateam/mafft/rbc_mafft/7.526+galaxy1')
    GALAXY_MAFFT_INPUT_KEY = os.environ.get('GALAXY_MAFFT_INPUT_KEY', 'inputSequences')
    GALAXY_MAFFT_PARAMS = os.environ.get(
        'GALAXY_MAFFT_PARAMS', '{"cond_flavour|flavourType": "mafft --auto"}')
    GALAXY_TRIMAL_TOOL_ID = os.environ.get(
        'GALAXY_TRIMAL_TOOL_ID',
        'toolshed.g2.bx.psu.edu/repos/iuc/trimal/trimal/1.4.1')
    GALAXY_TRIMAL_INPUT_KEY = os.environ.get('GALAXY_TRIMAL_INPUT_KEY', 'in')
    GALAXY_TRIMAL_PARAMS = os.environ.get(
        'GALAXY_TRIMAL_PARAMS', '{"trimming_mode|mode": "gappyout"}')
    GALAXY_RAXMLNG_TOOL_ID = os.environ.get('GALAXY_RAXMLNG_TOOL_ID', '')
