from flask import Flask


def create_app() -> Flask:
    app = Flask(__name__)
    app.config["SECRET_KEY"] = "some_key"

    from .views import views

    app.register_blueprint(views, url_prefix="/")
    return app
