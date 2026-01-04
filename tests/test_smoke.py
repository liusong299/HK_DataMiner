import subprocess
import os
import pytest

def test_cli_help():
    # Test that the CLI is installed and runnable
    result = subprocess.run(["hkdm", "--help"], capture_output=True)
    assert result.returncode == 0

def test_tutorial_help():
    # Test that tutorial command exists
    result = subprocess.run(["hkdm", "tutorial", "--help"], capture_output=True)
    assert result.returncode == 0
