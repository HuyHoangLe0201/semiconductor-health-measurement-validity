"""Verify the extracted repository payload without rerunning experiments."""
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parents[1]
manifest=json.loads((ROOT/'bundle_manifest.json').read_text())
for rel,expected in manifest.items():
    path=ROOT/rel
    assert path.is_file(),rel
    assert hashlib.sha256(path.read_bytes()).hexdigest()==expected,rel
print(f'All {len(manifest)} payload hashes match the upload manifest.')
