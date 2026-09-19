from pathlib import Path

import pytest


@pytest.fixture
def validation_data_dir(request):
    """Return the validation-case directory next to the requesting test."""
    test_file_path = Path(request.module.__file__)
    return test_file_path.parent / "validation_cases"
