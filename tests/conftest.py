import pytest
from pathlib import Path

@pytest.fixture
def test_data_dir(request):
    """@brief Fixture that provides the base path to the test files.
    @details
    This function uses the location of the current file (__file__) to 
    dynamically locate the 'data' folder. This ensures that the tests 
    work regardless of where the 'pytest' command is executed 
    in the terminal.

    The following directory structure is assumed:
        -tests/
            - preprocessing_tst/
                - files4test/
            - solver_tst/
                - files4test/
            - postrocessing_tst/
                - files4test/
    
    @param request Internal pytest fixture containing metadata 
            about the requesting test function.
    @return Path pathlib object pointing to the 'tests/XXXX/files4test' directory.
    """
    test_file_path = Path(request.module.__file__)
    return test_file_path.parent / "files4test"
