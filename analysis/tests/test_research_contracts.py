"""Checks that guard against evidence drift and overwritten run history."""
import json
import sys
import tempfile
import unittest
from pathlib import Path

ANALYSIS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ANALYSIS))
sys.path.insert(0, str(ANALYSIS / "scripts"))
from lib.provenance import archive_existing_record, sha256_file, write_json_atomic
from claim_contract import claim_family, parse_register, verify_bindings, generated


class ResearchContractTests(unittest.TestCase):
    def test_previous_record_bytes_survive_rerun(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "run.json"
            original = b'{"results": {"effect": 12}}\n'
            path.write_bytes(original)
            archive = archive_existing_record(path)
            write_json_atomic(path, {"results": {"effect": 4}})
            self.assertEqual(archive.read_bytes(), original)
            self.assertEqual(json.loads(path.read_text())["results"]["effect"], 4)
            self.assertNotEqual(sha256_file(path), sha256_file(archive))

    def test_changed_table_is_caught_without_rerunning_analysis(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "result.csv").write_text("arm,value\na,2\na,4\nb,99\n")
            binding = {"id": "test", "kind": "csv", "path": "result.csv", "where": {"arm": "a"},
                       "field": "value", "operation": "median", "expected": 3}
            self.assertFalse(verify_bindings([binding], root))
            (root / "result.csv").write_text("arm,value\na,2\na,8\nb,99\n")
            self.assertTrue(verify_bindings([binding], root))

    def test_missing_selected_rows_fail_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "result.csv").write_text("arm,value\nb,2\n")
            binding = {"id": "test", "kind": "csv", "path": "result.csv", "where": {"arm": "a"},
                       "field": "value", "operation": "median", "expected": 0}
            self.assertTrue(verify_bindings([binding], root))

    def test_historical_heading_does_not_reassign_cardoso_claims(self):
        for number in (65, 84, 105, 115):
            self.assertEqual(claim_family(number), "Cardoso and ligand extensions")

    def test_duplicate_claim_ids_fail(self):
        row = "| C1 | claim | design | `file.csv` | Descriptive only | caveat |\n"
        with self.assertRaises(ValueError):
            parse_register(row + row, Path.cwd())

    def test_displaced_section_is_not_an_active_validated_claim(self):
        text = "## Displaced\n| C2 | old claim | design | `file.csv` | Validated | outside work |\n"
        rows = parse_register(text, Path.cwd())
        self.assertEqual(rows[0]["status_group"], "Displaced")
        self.assertEqual(rows[0]["status"], "Validated")

    def test_narrative_number_cannot_drift_from_bound_artifact(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "analysis/claims").mkdir(parents=True)
            (root / "result.csv").write_text("value\n3\n")
            binding = {"id": "test", "claims": ["C1"], "kind": "csv", "path": "result.csv",
                       "field": "value", "operation": "single", "expected": 3,
                       "claim_text_template": "effect {value:.1f}"}
            (root / "analysis/claims/numeric_bindings.json").write_text(json.dumps({"bindings": [binding]}))
            (root / "CLAIMS.md").write_text("| C1 | effect 3.0 | design | `result.csv` | Descriptive only | caveat |\n")
            generated(root)
            (root / "CLAIMS.md").write_text("| C1 | effect 4.0 | design | `result.csv` | Descriptive only | caveat |\n")
            with self.assertRaises(ValueError):
                generated(root)


if __name__ == "__main__":
    unittest.main()
