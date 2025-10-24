#!/usr/bin/env python3
"""
Helper script to run Streamlit app from PyCharm.
This ensures Streamlit runs correctly without configuration issues.
"""
import sys
import subprocess

if __name__ == "__main__":
    # Run streamlit with proper arguments
    cmd = [sys.executable, "-m", "streamlit", "run", "app.py", "--server.port", "8501"]
    subprocess.run(cmd)
