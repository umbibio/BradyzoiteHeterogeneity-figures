#!/usr/bin/env python
"""Download and verify the input data listed in ``data/MANIFEST.json``.

Every file in the manifest carries a SHA-256 and a list of candidate URLs. The
URLs are tried in order until one yields a file with the expected digest, so a
dataset can be mirrored in several places (Zenodo, git LFS, an institutional
server) without changing any code.

Typical use::

    python scripts/fetch_data.py              # download whatever is missing
    python scripts/fetch_data.py --check      # verify what is already there
    python scripts/fetch_data.py --base-url https://example.org/bz/

Extra mirrors can also be supplied through the ``BZFIG_DATA_URLS`` environment
variable as a whitespace- or comma-separated list of base URLs.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import sys
import tempfile
import urllib.error
import urllib.request
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
MANIFEST = REPO / "data" / "MANIFEST.json"
CHUNK = 1 << 20


def read_manifest_file(datadir: Path = None) -> dict:
    """The manifest for a data directory. Kept out of LFS so it always reads."""
    path = MANIFEST if datadir is None else Path(datadir) / MANIFEST.name
    return json.loads(path.read_text())


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(CHUNK), b""):
            digest.update(block)
    return digest.hexdigest()


def base_urls(extra: list[str], manifest: dict) -> list[str]:
    """Mirror base URLs, highest priority first.

    Command line first, then ``BZFIG_DATA_URLS``, then the manifest's own list —
    so a user can override the published mirrors without editing anything.
    """
    env = os.environ.get("BZFIG_DATA_URLS", "")
    from_env = [u for u in env.replace(",", " ").split() if u]
    from_manifest = manifest.get("mirrors", [])
    urls = [*extra, *from_env, *from_manifest]
    return [u if u.endswith("/") else u + "/" for u in urls]


def candidates(entry: dict, bases: list[str]) -> list[str]:
    """Every URL worth trying for one manifest entry, in priority order."""
    urls = [base + entry["filename"] for base in bases]
    urls += [u for u in entry.get("urls", []) if u]
    seen, ordered = set(), []
    for url in urls:
        if url not in seen:
            seen.add(url)
            ordered.append(url)
    return ordered


def download(url: str, dest: Path) -> None:
    """Fetch ``url`` to ``dest`` via a temporary file, so a failure leaves no partial."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=dest.parent, delete=False) as tmp:
        tmp_path = Path(tmp.name)
    try:
        with urllib.request.urlopen(url, timeout=120) as response, tmp_path.open("wb") as out:
            shutil.copyfileobj(response, out, CHUNK)
        tmp_path.replace(dest)
    except BaseException:
        tmp_path.unlink(missing_ok=True)
        raise


def is_lfs_pointer(path: Path) -> bool:
    """True for a git-lfs pointer left behind by a clone without git-lfs.

    The pointer is a short text stub, so it fails the size check anyway; naming
    it lets the script say something more useful than "wrong size".
    """
    try:
        with path.open("rb") as handle:
            return handle.read(42).startswith(b"version https://git-lfs")
    except OSError:
        return False


def verify(path: Path, entry: dict) -> tuple[bool, str]:
    if not path.exists():
        return False, "missing"
    if is_lfs_pointer(path):
        return False, "git-lfs pointer, not the file itself"
    if path.stat().st_size != entry["bytes"]:
        return False, f"wrong size ({path.stat().st_size} != {entry['bytes']})"
    if sha256(path) != entry["sha256"]:
        return False, "sha256 mismatch"
    return True, "ok"


def fetch(entry: dict, outdir: Path, bases: list[str], force: bool) -> bool:
    path = outdir / entry["filename"]
    if not force:
        ok, why = verify(path, entry)
        if ok:
            print(f"  {entry['filename']}: already present, {why}")
            return True
        if why != "missing":
            print(f"  {entry['filename']}: {why}, re-downloading")

    urls = candidates(entry, bases)
    if not urls:
        print(f"  {entry['filename']}: MISSING and no URL known — see data/README.md")
        return False

    for url in urls:
        print(f"  {entry['filename']}: trying {url}")
        try:
            download(url, path)
        except (urllib.error.URLError, urllib.error.HTTPError, OSError) as exc:
            print(f"    failed: {exc}")
            continue
        ok, why = verify(path, entry)
        if ok:
            print(f"    ok ({entry['bytes']:,} bytes)")
            return True
        print(f"    downloaded but {why}; trying next mirror")
        path.unlink(missing_ok=True)

    print(f"  {entry['filename']}: could not be retrieved from any mirror")
    return False


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=MANIFEST)
    parser.add_argument("--outdir", type=Path, default=None, help="default: the manifest's directory")
    parser.add_argument(
        "--base-url",
        action="append",
        default=[],
        metavar="URL",
        help="mirror base URL; repeatable, tried before the manifest's own URLs",
    )
    parser.add_argument("--check", action="store_true", help="verify only, download nothing")
    parser.add_argument("--force", action="store_true", help="re-download even if the file verifies")
    args = parser.parse_args()

    manifest = json.loads(args.manifest.read_text())
    outdir = args.outdir or args.manifest.parent
    files = manifest["files"]

    if args.check:
        print(f"Verifying {len(files)} file(s) in {outdir}")
        failures = []
        for entry in files:
            ok, why = verify(outdir / entry["filename"], entry)
            print(f"  {entry['filename']}: {why}")
            if not ok:
                failures.append(entry["filename"])
        if failures:
            print(f"\n{len(failures)} file(s) failed verification: {', '.join(failures)}")
            print("Run scripts/fetch_data.py (no arguments) to download them.")
            return 1
        print("\nAll files verified.")
        return 0

    bases = base_urls(args.base_url, manifest)
    print(f"Fetching {len(files)} file(s) into {outdir}")
    if bases:
        print(f"Mirror bases: {', '.join(bases)}")
    missing = [entry["filename"] for entry in files if not fetch(entry, outdir, bases, args.force)]
    if missing:
        print(f"\n{len(missing)} file(s) unavailable: {', '.join(missing)}")
        print("Supply a mirror with --base-url or BZFIG_DATA_URLS, or see data/README.md")
        return 1
    print("\nAll files present and verified.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
