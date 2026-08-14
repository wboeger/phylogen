"""In-process, memory-only cache for user-supplied NCBI/Galaxy credentials.

Every job brings its own NCBI email/API key and Galaxy API key (see the
homepage notice) — this app has no shared server secret to protect and never
writes a user's key to the database or disk. Keys live here only long enough
for the background pipeline to use them, and expire on their own; a process
restart (e.g. a Railway redeploy) wipes them immediately, same as job data.
"""
import threading
import time

_TTL_SECONDS = 2 * 60 * 60  # 2 hours — enough to fetch/align/trim and then
                             # optionally click "Build ML tree" without retyping
_lock = threading.Lock()
_store = {}  # job_id -> (expires_at, dict)


def put(job_id, **creds):
    with _lock:
        _store[job_id] = (time.time() + _TTL_SECONDS, creds)


def get(job_id):
    with _lock:
        _purge_locked()
        entry = _store.get(job_id)
        return dict(entry[1]) if entry else None


def forget(job_id):
    with _lock:
        _store.pop(job_id, None)


def _purge_locked():
    now = time.time()
    expired = [k for k, (exp, _) in _store.items() if exp < now]
    for k in expired:
        _store.pop(k, None)
