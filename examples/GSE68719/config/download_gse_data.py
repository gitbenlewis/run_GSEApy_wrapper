#!/usr/bin/env python3
"""Download GSE68719 data files into examples/GSE68719/data."""

from __future__ import annotations

import gzip
import logging
from pathlib import Path
import shutil
from urllib.parse import parse_qs, unquote, urlsplit
from urllib.request import Request, urlopen


LOGGER = logging.getLogger(__name__)

GSE68719_URLS = [
    "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE68719&format=file&file=GSE68719%5Fmlpd%5FDESeq2%5Fdiffexp%5Fpmi%5Fage%5Frin%5Fdefault%2Etxt%2Egz",
    "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE68719&format=file&file=GSE68719%5Fmlpd%5FPCG%5FDESeq2%5Fnorm%5Fcounts%2Etxt%2Egz",
    "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE68719&format=file&file=GSE68719%5Fproteomics%5Fdataset%2ERData%2Egz",
]

CONFIG_DIR_PARENT = Path(__file__).resolve().parent.parent
DATA_DIR = CONFIG_DIR_PARENT / "data"


def _filename_from_url(url: str) -> str:
    query = urlsplit(url).query
    params = parse_qs(query)
    if "file" not in params or not params["file"]:
        raise ValueError(f"URL missing 'file' parameter: {url}")
    return unquote(params["file"][0])


def _download_url(url: str, dest_path: Path) -> None:
    if dest_path.exists() and dest_path.stat().st_size > 0:
        LOGGER.info("Skip download (exists): %s", dest_path)
        return
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = dest_path.with_suffix(dest_path.suffix + ".partial")
    req = Request(url, headers={"User-Agent": "run_GSEApy_wrapper/1.0"})
    LOGGER.info("Downloading %s -> %s", url, dest_path)
    with urlopen(req) as response, open(tmp_path, "wb") as out_handle:
        shutil.copyfileobj(response, out_handle)
    tmp_path.replace(dest_path)


def _gunzip_file(gz_path: Path, out_path: Path) -> None:
    if out_path.exists() and out_path.stat().st_size > 0:
        LOGGER.info("Skip unzip (exists): %s", out_path)
        return
    if not gz_path.exists():
        LOGGER.warning("Missing gzip source, cannot unzip: %s", gz_path)
        return
    LOGGER.info("Unzipping %s -> %s", gz_path, out_path)
    with gzip.open(gz_path, "rb") as in_handle, open(out_path, "wb") as out_handle:
        shutil.copyfileobj(in_handle, out_handle)


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    for url in GSE68719_URLS:
        filename = _filename_from_url(url)
        gz_path = DATA_DIR / filename
        _download_url(url, gz_path)
        if gz_path.suffix == ".gz":
            _gunzip_file(gz_path, gz_path.with_suffix(""))
        else:
            LOGGER.info("Skip unzip (not .gz): %s", gz_path)


if __name__ == "__main__":
    main()
