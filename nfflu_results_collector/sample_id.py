"""Parses a FastQID / sample name into its component fields: CID, Plate,
Index, Well.

Sample names are matched against three formats, in order, falling back to
the next if the previous doesn't match. This ordered fallback is the
single biggest source of historical bugs in this tool (see git log), so
each pattern is named and documented here rather than left as an inline
regex chain.

  1. standard: "{CID}-{Plate}-{Index}-{Well}"
     e.g. "C0000000001-2024-A1-A01"
  2. control:  "{CID}-{label}-{Plate}-{Index}-{Well}", where CID has no
     leading digit-only run (control/NTC samples are often labelled,
     e.g. "NTC20240101-Control-2024-B2-B02")
  3. salvage: neither of the above matched exactly (extra/missing
     characters, unexpected casing, etc). Recovers Plate/Index/Well from
     the trailing "-{plate}-{index}-{well}" suffix if present, and CID
     from a leading standard- or control-style prefix if present.
     Fields that can't be recovered are pd.NA; if even the CID can't be
     recovered, the full sample name is used as CID.
"""
import re

import pandas as pd

_STANDARD_RE = re.compile(r"^[A-Z][0-9]{10}-[0-9]{4,}-[A-Z0-9]{1,}-[A-Z][0-9]{2}$")
_CONTROL_RE = re.compile(r"^([A-Za-z]{3,}[0-9]{8})-[A-Za-z]+-([0-9]{4,})-([A-Z0-9]{1,})-([A-Z][0-9]{2})$")
_SALVAGE_SUFFIX_RE = re.compile(r"-[0-9]{4,9}-[A-Z0-9]{1,}-[A-Z][0-9]{2}$")
_SALVAGE_CID_SAMPLE_RE = re.compile(r"^([A-Z][0-9]{10}).*")
_SALVAGE_CID_CONTROL_RE = re.compile(r"^([A-Za-z]{3,}[0-9]{8}).*")

FIELDS = ["CID", "Plate", "Index", "Well"]


def parse_sample_id(sample_name):
    """Parse a sample name into {CID, Plate, Index, Well}. Plate/Index/Well
    are pd.NA if they can't be recovered; CID falls back to the full
    sample name if no recognizable prefix is found."""

    if _STANDARD_RE.match(sample_name):
        cid, plate, index, well = sample_name.split("-")
        return {"CID": cid, "Plate": plate, "Index": index, "Well": well}

    control_match = _CONTROL_RE.match(sample_name)
    if control_match:
        cid, plate, index, well = control_match.groups()
        return {"CID": cid, "Plate": plate, "Index": index, "Well": well}

    salvage_suffix_match = _SALVAGE_SUFFIX_RE.search(sample_name)
    if salvage_suffix_match:
        plate, index, well = salvage_suffix_match.group(0).lstrip("-").split("-")
    else:
        plate, index, well = pd.NA, pd.NA, pd.NA

    salvage_sample_cid_match = _SALVAGE_CID_SAMPLE_RE.search(sample_name)
    salvage_control_cid_match = _SALVAGE_CID_CONTROL_RE.search(sample_name)
    if salvage_sample_cid_match:
        cid = salvage_sample_cid_match.group(1)
    elif salvage_control_cid_match:
        cid = salvage_control_cid_match.group(1)
    else:
        cid = sample_name

    return {"CID": cid, "Plate": plate, "Index": index, "Well": well}
