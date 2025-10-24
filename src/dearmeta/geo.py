"""Interfaces for interacting with NCBI GEO to retrieve metadata and IDAT files."""

from __future__ import annotations

import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple
import requests
from pydantic import BaseModel, Field, validator

from .files import download_file, ensure_dir, gunzip_file, compute_md5
from .logging_utils import get_logger

logger = get_logger(__name__)

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
GEO_FTP_BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series"
GEO_GPL_BASE = "https://ftp.ncbi.nlm.nih.gov/geo/platforms"

IDAT_PATTERN = re.compile(r"/(GSM\d+)[^/]*_(Red|Grn)\.idat(\.gz)?$", re.IGNORECASE)

# Known GEO platform accessions for supported methylation arrays.
PLATFORM_ID_ALIASES = {
    "GPL13534": "HM450",  # Illumina HumanMethylation450 BeadChip
    "GPL21145": "EPICv1",  # Illumina HumanMethylationEPIC BeadChip
    "GPL33022": "EPICv2",  # Illumina MethylationEPIC v2.0 BeadChip (new GEO accession)
}


class SeriesMetadata(BaseModel):
    """Minimal GEO series metadata required for validation."""

    gse: str
    title: str
    summary: Optional[str] = None
    organism: str = Field(alias="taxon")
    platform_ids: List[str] = Field(default_factory=list, alias="gpl")
    platform_titles: List[str] = Field(default_factory=list, alias="gpltitle")
    sample_count: int = Field(default=0, alias="n_samples")

    class Config:
        allow_population_by_field_name = True

    @validator("gse", pre=True)
    def _validate_gse(cls, value: str) -> str:
        if not value:
            raise ValueError("Missing GSE accession")
        value = value.strip()
        if not value.upper().startswith("GSE"):
            value = f"GSE{value}"
        if not re.fullmatch(r"GSE\d+", value.upper()):
            raise ValueError(f"Invalid GSE accession: {value}")
        return value.upper()

    @validator("platform_ids", pre=True)
    def _normalise_platform_ids(cls, value):
        if value in (None, ""):
            return []
        if isinstance(value, str):
            items = [item.strip() for item in value.split(";") if item.strip()]
        elif isinstance(value, (list, tuple)):
            items = [str(item).strip() for item in value if str(item).strip()]
        else:
            raise ValueError("Unexpected format for platform IDs")
        normalised = []
        for item in items:
            item = item.upper()
            if not item.startswith("GPL"):
                item = f"GPL{item}" if item else item
            normalised.append(item)
        return normalised

    @validator("platform_titles", pre=True)
    def _ensure_list(cls, value):
        if value in (None, ""):
            return []
        if isinstance(value, str):
            return [item.strip() for item in value.split(";") if item.strip()]
        if isinstance(value, (list, tuple)):
            return [str(item).strip() for item in value if str(item).strip()]
        raise ValueError("Unexpected format for platform titles")

    @validator("sample_count", pre=True, always=True)
    def _clean_sample_count(cls, value):
        if value in (None, "", "NA"):
            return 0
        try:
            return int(value)
        except (TypeError, ValueError):
            return 0

    @property
    def is_human(self) -> bool:
        return self.organism.lower() == "homo sapiens"

    @property
    def platform_version(self) -> str:
        """Return inferred platform label (EPICv1/EPICv2/HM450)."""
        titles = [title.lower() for title in self.platform_titles if title]
        if titles:
            if any("epic v2" in title or "epicv2" in title for title in titles):
                return "EPICv2"
            if any("methylationepic" in title or "epic" in title for title in titles):
                return "EPICv1"
            if any("450k" in title or "methylation450" in title for title in titles):
                return "HM450"
        ids = [pid.upper() for pid in (self.platform_ids or []) if pid]
        for platform_id in ids:
            alias = PLATFORM_ID_ALIASES.get(platform_id)
            if alias:
                return alias
        titles_display = self.platform_titles or ["<missing>"]
        platform_ids_display = self.platform_ids or ["<missing>"]
        raise ValueError(
            "Unable to infer supported platform from metadata (titles: "
            + ", ".join(titles_display)
            + "; platform IDs: "
            + ", ".join(platform_ids_display)
            + ")."
        )

    def ensure_supported_platform(self) -> None:
        try:
            label = self.platform_version
        except ValueError as exc:
            raise ValueError(str(exc)) from exc
        if label not in {"EPICv1", "EPICv2", "HM450"}:
            raise ValueError(
                "Series uses platform(s) "
                + ", ".join(self.platform_titles or self.platform_ids or ["<unknown>"])
                + " which are not supported (expected Illumina EPIC v1/v2 or 450K)."
            )


class SampleMetadata(BaseModel):
    """Representation of a GEO sample (GSM) and associated metadata."""

    gsm: str
    organism: Optional[str] = None
    title: Optional[str] = None
    characteristics: Dict[str, str] = Field(default_factory=dict)
    supplementary: List[str] = Field(default_factory=list)

    @validator("gsm")
    def _validate_gsm(cls, value: str) -> str:
        if not re.fullmatch(r"GSM\d+", value):
            raise ValueError(f"Invalid GSM accession: {value}")
        return value.upper()


@dataclass
class IdatPair:
    gsm: str
    red: str
    green: str
    checksum_red: Optional[str] = None
    checksum_green: Optional[str] = None


class GeoClient:
    """Minimal GEO client built on top of the NCBI E-utilities."""

    def __init__(self, email: Optional[str] = None, api_key: Optional[str] = None):
        self.session = requests.Session()
        self.default_params = {"retmode": "json"}
        if email:
            self.default_params["email"] = email
        if api_key:
            self.default_params["api_key"] = api_key

    # --- Metadata discovery -------------------------------------------------
    def fetch_series(self, gse: str) -> SeriesMetadata:
        gse = gse.upper()
        gds_id = self._esearch_series(gse)
        summary = self._esummary(gds_id)
        metadata = SeriesMetadata.parse_obj(summary)
        if not metadata.platform_titles and metadata.platform_ids:
            metadata.platform_titles = self._fetch_platform_titles(metadata.platform_ids)
        metadata.ensure_supported_platform()
        if not metadata.is_human:
            raise ValueError(f"GEO series {gse} is not annotated as Homo sapiens.")
        return metadata

    def fetch_samples(self, gse: str, cache_dir: Optional[Path] = None, force: bool = False) -> Dict[str, SampleMetadata]:
        """Fetch sample metadata by parsing the series matrix file."""
        cache_dir = cache_dir or Path(".dearmeta_cache")
        matrix_path = self.download_series_matrix(gse, cache_dir=cache_dir, force=force)
        samples = parse_series_matrix(matrix_path)
        return {sample.gsm: sample for sample in samples}

    # --- IDAT discovery -----------------------------------------------------
    def discover_idat_pairs(
        self,
        gse: str,
        samples: Optional[Dict[str, SampleMetadata]] = None,
        cache_dir: Optional[Path] = None,
    ) -> Dict[str, IdatPair]:
        """Discover IDAT pairs for each GSM based on available supplementary files."""
        suppl_urls = self._list_suppl_files(gse)
        pairs = self._collect_idat_pairs(suppl_urls)
        if pairs:
            return pairs

        supplementary_urls: List[str] = []
        sample_map = samples
        if sample_map is None:
            sample_map = self.fetch_samples(gse, cache_dir=cache_dir)
        for sample in sample_map.values():
            supplementary_urls.extend(sample.supplementary)
        fallback_pairs = self._collect_idat_pairs(supplementary_urls)
        return fallback_pairs

    def _collect_idat_pairs(self, urls: Iterable[str]) -> Dict[str, IdatPair]:
        pairs: Dict[str, Dict[str, str]] = defaultdict(dict)
        for url in urls:
            if not url:
                continue
            match = IDAT_PATTERN.search(url)
            if not match:
                continue
            gsm, channel, _ = match.groups()
            channel = channel.lower()
            clean_url = _normalise_geo_url(url)
            if channel == "red":
                pairs[gsm]["red"] = clean_url
            else:
                pairs[gsm]["green"] = clean_url
        return {
            gsm: IdatPair(gsm=gsm, red=channels["red"], green=channels["green"])
            for gsm, channels in pairs.items()
            if {"red", "green"}.issubset(channels.keys())
        }

    # --- Download helpers ---------------------------------------------------
    def download_series_matrix(self, gse: str, cache_dir: Path, force: bool = False) -> Path:
        """Download and return the path to the series matrix file."""
        cache_dir = ensure_dir(cache_dir)
        target_dir = ensure_dir(cache_dir / gse)
        existing = list(target_dir.glob("*_series_matrix.txt"))
        if existing and not force:
            return existing[0]

        matrix_listing = self._list_remote_files(self._matrix_url(gse))
        matrix_candidates = [f for f in matrix_listing if f.endswith("_series_matrix.txt.gz")]
        if not matrix_candidates:
            raise RuntimeError(f"No series matrix found for {gse}")
        matrix_file = matrix_candidates[0]
        dest = target_dir / matrix_file
        checksum = download_file(f"{self._matrix_url(gse)}/{matrix_file}", dest)
        logger.info("Downloaded %s (md5=%s)", dest.name, checksum)
        extracted = gunzip_file(dest, remove_original=False)
        return extracted

    def download_platform_table(self, gpl: str, destination_dir: Path) -> Path:
        """Download GPL table to destination directory."""
        destination_dir = ensure_dir(destination_dir)
        group_dir = _platform_group_dir(gpl)
        url = f"{GEO_GPL_BASE}/{group_dir}/{gpl}/annot/"
        files = self._list_remote_files(url)
        candidates = [name for name in files if name.endswith("_annot.txt.gz")]
        if not candidates:
            raise RuntimeError(f"Could not locate annotation file for {gpl}")
        annot = candidates[0]
        dest = destination_dir / annot
        download_file(f"{url}{annot}", dest)
        return gunzip_file(dest, remove_original=False)

    def download_idat_pair(self, pair: IdatPair, dest_root: Path) -> IdatPair:
        """Download an IDAT pair for a GSM into ``dest_root``."""
        target_dir = ensure_dir(dest_root / pair.gsm)
        red_path = target_dir / Path(pair.red).name
        green_path = target_dir / Path(pair.green).name
        reused_red, md5_red = _reuse_existing_idat(red_path)
        reused_green, md5_green = _reuse_existing_idat(green_path)
        if reused_red:
            logger.debug("Reusing existing IDAT archive %s", red_path)
        else:
            md5_red = download_file(pair.red, red_path)
        if reused_green:
            logger.debug("Reusing existing IDAT archive %s", green_path)
        else:
            md5_green = download_file(pair.green, green_path)
        extracted_red = gunzip_file(red_path)
        extracted_green = gunzip_file(green_path)
        return IdatPair(
            gsm=pair.gsm,
            red=str(extracted_red),
            green=str(extracted_green),
            checksum_red=md5_red,
            checksum_green=md5_green,
        )

    # --- Internal helpers ---------------------------------------------------
    def _esearch_series(self, gse: str) -> str:
        params = {**self.default_params, "db": "gds", "term": f"{gse}[Accession]"}
        response = self.session.get(f"{EUTILS_BASE}/esearch.fcgi", params=params, timeout=30)
        response.raise_for_status()
        payload = response.json()
        ids = payload.get("esearchresult", {}).get("idlist", [])
        if not ids:
            raise ValueError(f"GEO series not found: {gse}")
        return ids[0]

    def _esummary(self, gds_id: str) -> dict:
        params = {**self.default_params, "db": "gds", "id": gds_id}
        response = self.session.get(f"{EUTILS_BASE}/esummary.fcgi", params=params, timeout=30)
        response.raise_for_status()
        payload = response.json()
        return payload["result"][gds_id]

    def _fetch_platform_titles(self, platform_ids: List[str]) -> List[str]:
        titles: List[str] = []
        for platform_id in platform_ids:
            try:
                params = {**self.default_params, "db": "gds", "id": platform_id}
                response = self.session.get(f"{EUTILS_BASE}/esummary.fcgi", params=params, timeout=30)
                response.raise_for_status()
                payload = response.json()
                summary = payload["result"].get(platform_id)
                if summary and summary.get("title"):
                    titles.append(summary["title"])
                else:  # pragma: no cover - fallback path
                    titles.append(platform_id)
            except requests.RequestException:  # pragma: no cover - network failure
                titles.append(platform_id)
        return titles

    def _list_suppl_files(self, gse: str) -> List[str]:
        url = f"{self._series_base(gse)}/suppl/"
        files = self._list_remote_files(url)
        return [f"{url}{name}" for name in files]

    def _list_remote_files(self, url: str) -> List[str]:
        response = self.session.get(url, timeout=30)
        response.raise_for_status()
        hrefs = re.findall(r'href="([^"]+)"', response.text)
        # remove parent directory references
        return [href for href in hrefs if not href.startswith("?") and href not in ("../", "./")]

    def _series_base(self, gse: str) -> str:
        group_dir = _series_group_dir(gse)
        return f"{GEO_FTP_BASE}/{group_dir}/{gse}"

    def _matrix_url(self, gse: str) -> str:
        return f"{self._series_base(gse)}/matrix"


def _normalise_geo_url(url: str) -> str:
    if url.startswith("ftp://"):
        return "https://" + url[len("ftp://") :]
    return url


def _reuse_existing_idat(path: Path) -> Tuple[bool, Optional[str]]:
    """Return whether a previously downloaded IDAT archive can be reused."""
    extracted = path.with_suffix("") if path.suffix == ".gz" else path
    if path.exists() and extracted.exists():
        return True, compute_md5(path)
    return False, None


def _series_group_dir(accession: str) -> str:
    """Return the GEO FTP grouping directory for an accession."""
    number = int(re.sub(r"\D", "", accession))
    group = number // 1000
    if group == 0:
        return "GSEnnn"
    return f"GSE{group}nnn"


def _platform_group_dir(accession: str) -> str:
    number = int(re.sub(r"\D", "", accession))
    group = number // 1000
    if group == 0:
        return "GPLnnn"
    return f"GPL{group}nnn"


def parse_series_matrix(path: Path) -> List[SampleMetadata]:
    """Parse a GEO series matrix file to extract sample metadata."""
    characteristics: Dict[str, Dict[str, str]] = defaultdict(dict)
    titles: Dict[str, str] = {}
    organisms: Dict[str, str] = {}
    supplementary: Dict[str, List[str]] = defaultdict(list)

    gsm_order: List[str] = []

    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith("!Sample_geo_accession"):
                gsm_order = _split_series_line(line)
            elif line.startswith("!Sample_title"):
                titles = dict(zip(gsm_order, _split_series_line(line)))
            elif line.startswith("!Sample_organism_ch1"):
                organisms = dict(zip(gsm_order, _split_series_line(line)))
            elif line.startswith("!Sample_characteristics_ch1"):
                values = _split_series_line(line)
                for gsm, value in zip(gsm_order, values):
                    key, val = _split_characteristic(value)
                    if key:
                        characteristics[gsm][key] = val
            elif line.startswith("!Sample_supplementary_file"):
                values = _split_series_line(line)
                for gsm, value in zip(gsm_order, values):
                    if value and value != "NONE":
                        supplementary[gsm].append(value)

    samples: List[SampleMetadata] = []
    for gsm in gsm_order:
        samples.append(
            SampleMetadata(
                gsm=gsm,
                title=titles.get(gsm),
                organism=organisms.get(gsm),
                characteristics=characteristics.get(gsm, {}),
                supplementary=supplementary.get(gsm, []),
            )
        )
    return samples


def _split_series_line(line: str) -> List[str]:
    parts = line.split("\t")
    # drop header key
    values: List[str] = []
    for raw in parts[1:]:
        value = raw.strip()
        if value.startswith('"') and value.endswith('"') and len(value) >= 2:
            value = value[1:-1]
        elif value.startswith("'") and value.endswith("'") and len(value) >= 2:
            value = value[1:-1]
        values.append(value.strip())
    return values


def _split_characteristic(value: str) -> Tuple[str, str]:
    if ":" in value:
        key, val = value.split(":", 1)
        key = key.strip().lower().replace(" ", "_")
        return key, val.strip()
    return value.strip().lower().replace(" ", "_"), ""
