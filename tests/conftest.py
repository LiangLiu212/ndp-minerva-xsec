"""Put the repo root on sys.path so tests import the xsec package in place."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
