#!/usr/bin/env python3
"""Small, dependency-free clients for external variant annotation services.

The public functions return dictionaries and retain provider responses under
``data``.  They deliberately do not translate provider fields into the local
database schema; that belongs in a separate annotation/import layer.

Ensembl and ClinVar provide documented JSON APIs.  MITOMAP currently exposes a
human-facing allele search rather than a documented JSON API, so its adapter
returns the response text and content type without attempting fragile parsing.
"""

import argparse
import json
import os
import time
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.parse import quote, urlencode
from urllib.request import Request, urlopen


ENSEMBL_REST_URL = "https://rest.ensembl.org"
NCBI_EUTILS_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
MITOMAP_ALLELE_SEARCH_URL = "https://www.mitomap.org/allelesearch.html"
DEFAULT_TIMEOUT = 30.0
USER_AGENT = "Rohrer-Barb-Mito-Annotation/1.0"


class AnnotationAPIError(RuntimeError):
    """Raised when an annotation provider cannot return a usable response."""


def _request(
    url: str,
    *,
    method: str = "GET",
    payload: dict[str, Any] | None = None,
    headers: dict[str, str] | None = None,
    timeout: float = DEFAULT_TIMEOUT,
) -> tuple[Any, str]:
    """Request a URL and return its decoded body and content type."""
    request_headers = {
        "Accept": "application/json",
        "User-Agent": USER_AGENT,
        **(headers or {}),
    }
    body = None
    if payload is not None:
        body = json.dumps(payload).encode("utf-8")
        request_headers["Content-Type"] = "application/json"

    request = Request(url, data=body, headers=request_headers, method=method)
    try:
        with urlopen(request, timeout=timeout) as response:
            raw = response.read().decode("utf-8", errors="replace")
            content_type = response.headers.get_content_type()
    except HTTPError as error:
        detail = error.read().decode("utf-8", errors="replace")[:1000]
        raise AnnotationAPIError(
            f"HTTP {error.code} from {url}: {detail or error.reason}"
        ) from error
    except (URLError, TimeoutError) as error:
        raise AnnotationAPIError(f"Request failed for {url}: {error}") from error

    if content_type == "application/json" or raw.lstrip().startswith(("{", "[")):
        try:
            return json.loads(raw), content_type
        except json.JSONDecodeError as error:
            raise AnnotationAPIError(f"Invalid JSON returned by {url}") from error
    return raw, content_type


def _envelope(source: str, url: str, data: Any, **details: Any) -> dict[str, Any]:
    """Wrap a provider response with enough provenance to cache it safely."""
    return {
        "source": source,
        "retrieved_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "request_url": url,
        **details,
        "data": data,
    }


def fetch_ensembl_vep(
    chrom: str,
    position: int,
    ref: str,
    alt: str,
    *,
    species: str = "homo_sapiens",
    timeout: float = DEFAULT_TIMEOUT,
) -> dict[str, Any]:
    """Return raw Ensembl VEP annotations for one VCF-style allele.

    The main Ensembl REST service annotates the current assembly (GRCh38 for
    human). Use a release-specific/offline VEP cache when reproducibility across
    Ensembl releases is required.
    """
    chrom = str(chrom).removeprefix("chr")
    variant = f"{chrom} {int(position)} . {ref.upper()} {alt.upper()} . . ."
    url = f"{ENSEMBL_REST_URL}/vep/{quote(species, safe='')}/region"
    data, _ = _request(
        url,
        method="POST",
        payload={"variants": [variant]},
        timeout=timeout,
    )
    return _envelope(
        "ensembl_vep",
        url,
        data,
        query={"chrom": chrom, "position": int(position), "ref": ref, "alt": alt},
    )


def mitochondrial_hgvs(position: int, ref: str, alt: str) -> str:
    """Build an rCRS mitochondrial HGVS expression for a substitution."""
    if len(ref) != 1 or len(alt) != 1:
        raise ValueError("Automatic mitochondrial HGVS currently supports SNVs only")
    return f"NC_012920.1:m.{int(position)}{ref.upper()}>{alt.upper()}"


def fetch_clinvar(
    query: str,
    *,
    api_key: str | None = None,
    email: str | None = None,
    tool: str = "rohrer_barb_mito_annotation",
    timeout: float = DEFAULT_TIMEOUT,
) -> dict[str, Any]:
    """Search ClinVar and return JSON document summaries for matching records."""
    common = {"db": "clinvar", "retmode": "json", "tool": tool}
    api_key = api_key or os.environ.get("NCBI_API_KEY")
    email = email or os.environ.get("NCBI_EMAIL")
    if api_key:
        common["api_key"] = api_key
    if email:
        common["email"] = email

    search_params = {**common, "term": query, "retmax": "100"}
    search_url = f"{NCBI_EUTILS_URL}/esearch.fcgi?{urlencode(search_params)}"
    search_data, _ = _request(search_url, timeout=timeout)
    ids = search_data.get("esearchresult", {}).get("idlist", [])

    summaries: dict[str, Any] = {"result": {"uids": []}}
    summary_url = None
    if ids:
        summary_params = {**common, "id": ",".join(ids)}
        summary_url = f"{NCBI_EUTILS_URL}/esummary.fcgi?{urlencode(summary_params)}"
        summaries, _ = _request(summary_url, timeout=timeout)

    return _envelope(
        "clinvar",
        search_url,
        {"search": search_data, "summaries": summaries},
        query=query,
        summary_url=summary_url,
    )


def fetch_clinvar_mito_variant(
    position: int,
    ref: str,
    alt: str,
    **kwargs: Any,
) -> dict[str, Any]:
    """Search ClinVar for an rCRS mitochondrial SNV by its HGVS name."""
    hgvs = mitochondrial_hgvs(position, ref, alt)
    return fetch_clinvar(f'"{hgvs}"[Variant Name]', **kwargs)


def fetch_mitomap(
    position: int,
    ref: str | None = None,
    alt: str | None = None,
    *,
    timeout: float = DEFAULT_TIMEOUT,
) -> dict[str, Any]:
    """Retrieve MITOMAP's public allele-search page as unparsed text.

    MITOMAP does not document a JSON annotation API. This intentionally avoids
    presenting scraped HTML as stable structured annotation data. ``data`` is
    therefore a dictionary containing the returned HTML and content type.
    """
    notation = str(int(position))
    if ref and alt:
        notation = f"m.{int(position)}{ref.upper()}>{alt.upper()}"
    params = urlencode({"search": notation})
    url = f"{MITOMAP_ALLELE_SEARCH_URL}?{params}"
    raw, content_type = _request(
        url,
        headers={"Accept": "text/html,application/xhtml+xml"},
        timeout=timeout,
    )
    return _envelope(
        "mitomap",
        url,
        {
            "content_type": content_type,
            "html": raw,
            "structured_api": False,
            "warning": "MITOMAP has no documented public JSON annotation API.",
        },
        query=notation,
    )


def fetch_all_mito_annotations(
    position: int,
    ref: str,
    alt: str,
    *,
    continue_on_error: bool = True,
    timeout: float = DEFAULT_TIMEOUT,
) -> dict[str, Any]:
    """Fetch all three providers, isolating provider failures by default."""
    providers = {
        "ensembl": lambda: fetch_ensembl_vep(
            "MT", position, ref, alt, timeout=timeout
        ),
        "clinvar": lambda: fetch_clinvar_mito_variant(
            position, ref, alt, timeout=timeout
        ),
        "mitomap": lambda: fetch_mitomap(position, ref, alt, timeout=timeout),
    }
    results: dict[str, Any] = {}
    for name, fetch in providers.items():
        try:
            results[name] = fetch()
        except (AnnotationAPIError, ValueError) as error:
            if not continue_on_error:
                raise
            results[name] = {"source": name, "error": str(error)}
    return results


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Fetch raw mitochondrial variant annotations as JSON."
    )
    parser.add_argument("position", type=int)
    parser.add_argument("ref")
    parser.add_argument("alt")
    parser.add_argument("--source", choices=("all", "ensembl", "clinvar", "mitomap"), default="all")
    parser.add_argument("--timeout", type=float, default=DEFAULT_TIMEOUT)
    return parser


def main() -> None:
    args = build_argument_parser().parse_args()
    if args.source == "all":
        result = fetch_all_mito_annotations(
            args.position, args.ref, args.alt, timeout=args.timeout
        )
    elif args.source == "ensembl":
        result = fetch_ensembl_vep(
            "MT", args.position, args.ref, args.alt, timeout=args.timeout
        )
    elif args.source == "clinvar":
        result = fetch_clinvar_mito_variant(
            args.position, args.ref, args.alt, timeout=args.timeout
        )
    else:
        result = fetch_mitomap(
            args.position, args.ref, args.alt, timeout=args.timeout
        )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
