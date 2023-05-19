import os
import sys
import pytest

from t3.settings import settings as t3_settings
from t3.settings.t3_submit import submit_scripts


@pytest.fixture(autouse=True)
def mock_t3_files(monkeypatch):
    # Mocking the existence of local settings and submit_scripts files
    monkeypatch.setattr(os.path, "isfile", lambda path: True)
    monkeypatch.setattr(os.path, "join", lambda *args: "/mock/path")

    # Mocking the import of local settings
    monkeypatch.setitem(sys.modules, "settings", t3_settings)

    # Mocking the import of local submit_scripts
    monkeypatch.setattr(sys, "path", ["/mock/path"])

    def mock_import_submit_scripts():
        submit_scripts.update({"mock_script": "mock_content"})

    monkeypatch.setattr(submit_scripts, "update", mock_import_submit_scripts)


def test_import_user_settings():
    import user_settings

    # Check that settings dictionary is updated with local settings
    assert "local_settings" in user_settings.settings
    assert user_settings.settings["local_settings"] == {"mock_key": "mock_value"}

    # Check that submit_scripts dictionary is updated with local submit_scripts
    assert "mock_script" in user_settings.submit_scripts
    assert user_settings.submit_scripts["mock_script"] == "mock_content"


def test_import_user_settings_no_local_files(monkeypatch):
    # Mocking the absence of local settings and submit_scripts files
    monkeypatch.setattr(os.path, "isfile", lambda path: False)

    import user_settings

    # Check that settings dictionary remains the same
    assert user_settings.settings == {
        key: val for key, val in vars(t3_settings).items() if "__" not in key
    }

    # Check that submit_scripts dictionary remains the same
    assert user_settings.submit_scripts == submit_scripts
