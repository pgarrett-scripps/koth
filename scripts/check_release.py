#!/usr/bin/env python3
"""Check static citation policy, documentation assets, and release identity."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import tomllib
from urllib.parse import unquote, urlsplit

import yaml

ROOT = Path(__file__).resolve().parents[1]
REPOSITORY = "https://github.com/pgarrett-scripps/koth"
REPOSITORY_NAME = "pgarrett-scripps/koth"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def check_metadata(root):
    require((root / "LICENSE").read_bytes() == (root / "koth_ff/LICENSE").read_bytes(),
            "Root and crate licenses differ")
    zen = json.loads((root / ".zenodo.json").read_text())
    cff = yaml.safe_load((root / "CITATION.cff").read_text())
    # Zenodo takes the version and date from the release. CITATION.cff carries
    # them, written by `just cite-sync` and checked against Cargo.toml below.
    funding = {"grants", "funding", "funders", "funding-references", "funding_references"}
    forbidden = funding | {"version", "date-released", "publication_date"}

    def check_keys(value, banned):
        if isinstance(value, dict):
            require(not banned.intersection(value),
                    "Citation metadata must not hardcode " + ", ".join(sorted(banned.intersection(value))))
            for item in value.values():
                check_keys(item, banned)
        elif isinstance(value, list):
            for item in value:
                check_keys(item, banned)

    check_keys(zen, forbidden)
    check_keys(cff, funding)
    version = tomllib.loads((root / "Cargo.toml").read_text())["workspace"]["package"]["version"]
    require(str(cff.get("version")) == version,
            "CITATION.cff version differs from Cargo.toml; run `just cite-sync`")
    require(re.fullmatch(r"\d{4}-\d{2}-\d{2}", str(cff.get("date-released", ""))),
            "CITATION.cff date-released is missing; run `just cite-sync`")
    require(cff["cff-version"] == "1.2.0", "Unsupported CFF format")
    require(zen["upload_type"] == cff["type"] == "software", "Metadata must describe software")
    require(zen["title"] == cff["title"], "Citation titles differ")
    require(zen["description"] == cff["abstract"], "Citation descriptions differ")
    require(zen["license"] == cff["license"] == "MIT", "License mismatch")
    require(zen["keywords"] == cff["keywords"], "Citation keywords differ")
    require(cff["repository-code"] == cff["url"] == REPOSITORY, "Incorrect repository URL")
    require(any(x["identifier"] == REPOSITORY for x in zen["related_identifiers"]),
            "Zenodo must reference this repository")
    require(len(cff["authors"]) == len(zen["creators"]) > 0, "Creator lists differ")
    for author, creator in zip(cff["authors"], zen["creators"]):
        name = author["family-names"] + ", " + author["given-names"]
        if "name-suffix" in author:
            name += ", " + author["name-suffix"]
        require(name == creator["name"], "Creator names/order differ")
        orcid = creator.get("orcid")
        require(author.get("orcid") == (orcid and "https://orcid.org/" + orcid), "ORCIDs differ")
        require(author["affiliation"] == creator["affiliation"], "Affiliations differ")


def check_docs(root):
    # Check file destinations, including this repository's absolute GitHub links.
    paths = [root / "README.md", root / "RELEASE.md", root / "CONTRIBUTING.md"]
    paths += sorted((root / "docs").glob("*.md"))
    for path in paths:
        for destination in re.findall(r"\]\(([^)]+)\)", path.read_text()):
            destination = destination.split(' "', 1)[0]
            prefix = REPOSITORY + "/blob/master/"
            image_prefix = "https://raw.githubusercontent.com/" + REPOSITORY_NAME + "/master/"
            if destination.startswith(prefix):
                base, destination = root, destination[len(prefix):]
            elif destination.startswith(image_prefix):
                base, destination = root, destination[len(image_prefix):]
            elif urlsplit(destination).scheme or destination.startswith("#"):
                continue
            else:
                base = path.parent
            target = base / unquote(destination.split("#", 1)[0])
            require(target.exists(), f"Broken documentation link in {path.name}: {destination}")
    provenance = json.loads((root / "docs/assets/provenance.json").read_text())
    for item in provenance["figures"]:
        digest = hashlib.sha256((root / item["path"]).read_bytes()).hexdigest()
        require(digest == item["sha256"], f"Figure differs from provenance: {item['path']}")


def check_event(root, event):
    require(event.get("action") == "published", "Only published GitHub releases may publish packages")
    require(event["repository"]["full_name"] == REPOSITORY_NAME, "Unexpected release repository")
    require(event["repository"]["private"] is False, "Make the repository public before a public package release")
    release = event["release"]
    require(release["draft"] is False, "Draft releases must not publish")
    version = tomllib.loads((root / "Cargo.toml").read_text())["workspace"]["package"]["version"]
    require(release["tag_name"] == "v" + version, "Release tag must equal v<Cargo.toml version>")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--event", type=Path, help="GitHub published-release payload to validate")
    args = parser.parse_args()
    check_metadata(ROOT)
    check_docs(ROOT)
    if args.event:
        check_event(ROOT, json.loads(args.event.read_text()))
    print("Citation policy, documentation links, figure hashes, and requested release checks passed.")


if __name__ == "__main__":
    main()
