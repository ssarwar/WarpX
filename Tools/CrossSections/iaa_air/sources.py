# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Read physical source records without importing or executing the source package."""

import hashlib
import re
import tarfile

import numpy as np

ARCHIVE_SHA256 = "b864381086120fff37eb8c658dfcb0f150b80969cf9add626a38d89bab0a3c38"
N2_STATES = {
    "A3Σu+": "A",
    "B3Πg": "B",
    "W3Δu": "W",
    "B'3Σu-": "B_prime",
    "a'1Σu-": "a_prime",
    "a1Πg": "a",
    "w1Δu": "w",
    "C3Πu": "C",
    "E3Σg+": "E",
    'a"1Σg+': "a_double_prime",
    "b1Πu": "b",
    "c31Πu": "c3",
    "o31Πu": "o3",
    "b'1Σu+": "b_prime",
    "c'41Σu+": "c4_prime",
    "G3Πu": "G",
    "F3Πu": "F",
}
O2_STATES = {
    "a1Δg": "a",
    "b1Σg+": "b",
    "HC": "herzberg",
    "SR": "schumann_runge",
    "LB": "longest_band",
    "2B": "second_band",
    "?": "complement_11",
}


def digest(data):
    return hashlib.sha256(data).hexdigest()


def file_stem(channel):
    """Distinguish A/a, B/b and W/w on case-insensitive filesystems."""
    return "triplet_" + channel.lower() if channel[0].isupper() else channel


class Archive:
    def __init__(self, path):
        if digest(path.read_bytes()) != ARCHIVE_SHA256:
            raise ValueError(
                "Unrecognized elmolcs archive: review source mappings before export"
            )
        self.archive = tarfile.open(path)

    def read(self, path):
        return self.archive.extractfile("elmolcs/" + path).read().decode()


def blocks(text):
    """Read the numerical IAA/LXCat blocks, preserving original metadata."""
    lines = text.splitlines()
    for i, line in enumerate(lines):
        match = re.fullmatch(r"SPECIES:\s*e\s*/\s*(N2|O2)\s*", line)
        if not match:
            continue
        j = i
        while j < len(lines) and not re.match(r"-{5,}", lines[j]):
            j += 1
        k = j + 1
        while k < len(lines) and not re.match(r"-{5,}", lines[k]):
            k += 1
        header = "\n".join(lines[i:j])
        process = re.search(r"^PROCESS:\s*(.*)$", header, re.M)
        if not process:
            continue
        table = []
        for row in lines[j + 1 : k]:
            row = row.partition("#")[0].strip()
            if not row:
                continue
            try:
                values = [float(x) for x in row.split()]
            except ValueError:
                continue
            if len(values) == 2:
                table.append(values)
        if len(table) < 2:
            continue
        previous = lines[i - 1].strip()
        try:
            loss = float(previous)
        except ValueError:
            loss = 0.0
        yield {
            "target": match[1],
            "original_process": process[1],
            "header": header,
            "source_line": i + 1,
            "loss_eV": loss,
            "table": table,
        }


def read_sources(archive, attachment):
    """Return physical records. Aggregate and rotational rates are excluded."""
    text = archive.read("Data/cs/Cross section.txt")
    result = {"N2": [], "O2": []}
    for record in blocks(text):
        target = record["target"]
        process = record["original_process"]
        kind = process.rsplit(",", 1)[-1].strip()
        header = record["header"]
        record["source_file"] = "Data/cs/Cross section.txt"
        record["source_sha256"] = digest(text.encode())
        record["units"] = "m2"
        record["type"] = "excitation"
        record["model"] = "table"
        record["alternative_group"] = None
        if kind == "Vibrational":
            level = re.search(r"v=0-(\d+|None)", process)[1]
            record["id"] = (
                "vibrational_dissociation" if level == "None" else f"v{int(level):02d}"
            )
            record["family"] = "dissociation" if level == "None" else "vibrational"
        elif kind == "Dissociation":
            record["id"] = "residual_dissociation"
            record["family"] = "dissociation"
        elif kind == "Dissociative":
            record.update(
                id="attachment_iaa",
                family="attachment",
                type="attachment",
                alternative_group="dissociative_attachment",
                loss_eV=0.0,
            )
        elif kind == "Complementary":
            record["id"] = f"complement_{int(record['loss_eV'])}"
            record["family"] = "complementary"
        elif kind == "Electronic":
            state = re.search(r"\*\(([^)]+)\)", process)[1]
            record["id"] = (N2_STATES if target == "N2" else O2_STATES)[state]
            if state == "SR":
                record["id"] += "_onset" if "onset" in header else "_bulk"
            record["family"] = "complementary" if state == "?" else "electronic"
            coefficients = re.search(r"c\s*=\s*[\[(]([^\])]+)[\])]", header, re.S)[1]
            record["fit"] = {
                "coefficients": [
                    float(x) for x in coefficients.replace(",", " ").split()
                ],
                "threshold_eV": float(re.search(r"Eth\s*=\s*([\d.eE+-]+)", header)[1]),
                "asymptotic_order": int(re.search(r"asymp\s*=\s*(\d+)", header)[1]),
                "equation": re.findall(
                    r"[Ee]q\.?(?:\s*)=(?:\s*)['\"]?([\w+]+)", header
                )[-1],
            }
            record["model"] = "analytic"
        else:
            continue
        energy, sigma = np.asarray(record["table"]).T
        if np.any(np.diff(energy) <= 0) or np.any(sigma < 0):
            raise ValueError(f"Invalid source rows for {target}/{record['id']}")
        record["source_range_eV"] = [float(energy[0]), float(energy[-1])]
        record["tabulated_onset_eV"] = float(energy[0])
        nominal = re.search(r"Eth\s*=\s*([\d.eE+-]+)", header)
        record["nominal_threshold_eV"] = (
            float(nominal[1]) if nominal else record["loss_eV"]
        )
        if "fit" in record:
            fit = record["fit"]
            coefficients = fit["coefficients"]
            record["fitted_threshold_eV"] = (
                coefficients[2]
                if fit["equation"] == "f0m"
                else coefficients[4]
                if "a" in fit["equation"] and len(coefficients) == 5
                else fit["threshold_eV"]
            )
        record["repairs"] = []
        if target == "N2" and record["id"] == "v01":
            # Source file Data/cs/N2/Cross section.e-N2, Linert/Tanaka tail.
            original = archive.read("Data/cs/N2/Cross section.e-N2")
            assert re.search(r"^\s*30\s+2\.66e-22\s*$", original, re.M)
            assert re.search(r"^\s*20\s+122\.5e-23", original, re.M)
            record["table"].extend([[20.0, 122.5e-23], [30.0, 2.66e-22]])
            record["table"].sort()
            record["source_range_eV"][1] = 30.0
            record["repairs"].append(
                "Restore 20 eV, 1.225e-21 m2 and 30 eV, 2.66e-22 m2 from analytic-source file"
            )
        if record["id"].startswith("schumann_runge"):
            record["partition"] = "onset" if record["id"].endswith("onset") else "bulk"
            record["repairs"].append(
                "Partition one analytic total over 9.955581--10.213826 eV"
            )
        result[target].append(record)

    for source, text in [("muroran", attachment.read_text())] + [
        (name, archive.read(f"Data/cs/O2/others/Cross section ({name}).e-O2"))
        for name in ["phelps", "biagi"]
    ]:
        for record in blocks(text):
            if not record["original_process"].endswith("Attachment"):
                continue
            three_body = "O2-" in record["original_process"]
            if not three_body and source != "muroran":
                continue
            factor = 1e-6 if three_body and source != "muroran" else 1.0
            record["table"] = [[e, v * factor] for e, v in record["table"]]
            record.update(
                id=("three_body_" if three_body else "attachment_") + source,
                family="attachment",
                type="attachment",
                loss_eV=0.0,
                model="table",
                units="m5" if three_body else "m2",
                third_body="O2" if three_body else None,
                alternative_group="three_body_attachment"
                if three_body
                else "dissociative_attachment",
                source_file=attachment.name
                if source == "muroran"
                else f"Data/cs/O2/others/Cross section ({source}).e-O2",
                source_sha256=digest(text.encode()),
                conversion_factor=factor,
                repairs=[],
                tabulated_onset_eV=record["table"][0][0],
                nominal_threshold_eV=0.0 if three_body else 3.64,
                source_range_eV=[record["table"][0][0], record["table"][-1][0]],
            )
            result["O2"].append(record)
    for gas, count in [("N2", 35), ("O2", 55)]:
        records = result[gas]
        if len(records) != count or len({r["id"] for r in records}) != count:
            raise ValueError(
                f"Unexpected {gas} inventory: {[r['id'] for r in records]}"
            )
    return result
