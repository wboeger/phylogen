from flask import Blueprint, g, render_template

from app import db
from app.models import Job
from app.pipeline import MARKER_PRESETS

main_bp = Blueprint('main', __name__)


@main_bp.route('/')
def index():
    jobs = (Job.query.filter_by(session_id=g.session_id)
            .order_by(Job.created_at.desc()).limit(20).all())
    return render_template('index.html', jobs=jobs, markers=MARKER_PRESETS)
