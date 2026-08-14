import secrets

from flask import Flask, g, request
from flask_sqlalchemy import SQLAlchemy

db = SQLAlchemy()

SESSION_COOKIE = 'phylo_sid'
SESSION_MAX_AGE = 60 * 60 * 24 * 30  # 30 days


def create_app():
    app = Flask(__name__)
    app.config.from_object('config.Config')

    db.init_app(app)

    from app.routes.main import main_bp
    from app.routes.jobs import jobs_bp
    app.register_blueprint(main_bp)
    app.register_blueprint(jobs_bp)

    with app.app_context():
        db.create_all()

    @app.before_request
    def _ensure_session():
        sid = request.cookies.get(SESSION_COOKIE)
        g.new_session = sid is None
        g.session_id = sid or secrets.token_urlsafe(24)

    @app.after_request
    def _set_session_cookie(response):
        if getattr(g, 'new_session', False):
            response.set_cookie(
                SESSION_COOKIE, g.session_id,
                max_age=SESSION_MAX_AGE, httponly=True, samesite='Lax',
            )
        return response

    return app
