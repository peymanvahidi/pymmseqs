import os
import unittest
import subprocess
from pathlib import Path
import tempfile

from pymmseqs.config import CreateDBConfig, EasyLinClustConfig
from pymmseqs.parsers.base_cluster_parser import BaseClusterParser


class TestCreateDB(unittest.TestCase):
    def test_createdb_output_matches_cli(self):
        """
        Test that our createdb function produces identical output to mmseqs CLI
        """
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)

            # Create test FASTA file
            fasta_content = ">seq1\nAAAA\n>seq2\nCCCC\n"
            fasta_file = tmp_path / "input_test.fasta"
            fasta_file.write_text(fasta_content)

            # Define output paths
            func_output = tmp_path / "func_output" / "mydb"
            cli_output = tmp_path / "cli_output" / "mydb"

            # 1. Run our Python function implementation
            config = CreateDBConfig(
                fasta_file=fasta_file,
                sequence_db=func_output,
                write_lookup=1  # Match CLI default behavior
            )
            config.run()

            # 2. Run mmseqs CLI command
            # Create parent directory for CLI output
            cli_output.parent.mkdir(parents=True, exist_ok=True)

            subprocess.run(
                [
                    "mmseqs",
                    "createdb",
                    str(fasta_file),
                    str(cli_output),
                    "--write-lookup", "1"  # Match our function's default
                ],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE
            )

            # 3. Compare outputs
            # Get all generated files from both implementations
            func_files = set(func_output.parent.glob("mydb*"))
            cli_files = set(cli_output.parent.glob("mydb*"))

            # Check same number of files created
            self.assertEqual(
                len(func_files),
                len(cli_files),
                "Different number of output files generated"
            )

            # Compare each pair of files
            for func_file in func_files:
                filename = func_file.name
                cli_file = cli_output.parent / filename

                # Verify file exists in CLI output
                self.assertTrue(
                    cli_file.exists(),
                    f"File {filename} missing in CLI output"
                )

                # Compare file contents
                with open(func_file, "rb") as f1, open(cli_file, "rb") as f2:
                    self.assertEqual(
                        f1.read(),
                        f2.read(),
                        f"Content mismatch in {filename}"
                    )


class TestEasyLinClust(unittest.TestCase):
    def test_easy_linclust_output_matches_cli(self):
        """
        Test that our easy_linclust config produces identical output to mmseqs CLI
        """
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)

            # Create test FASTA file with enough sequences to form clusters
            fasta_content = (
                ">seq1\nACDEFGHIKLMNPQRSTVWY\n"
                ">seq2\nACDEFGHIKLMNPQRSTVWY\n"
                ">seq3\nWWWWWWWWWWWWWWWWWWWW\n"
            )
            fasta_file = tmp_path / "input_test.fasta"
            fasta_file.write_text(fasta_content)

            # Define output paths
            func_output = tmp_path / "func_output" / "result"
            cli_output = tmp_path / "cli_output" / "result"

            func_tmp = tmp_path / "func_tmp"
            cli_tmp = tmp_path / "cli_tmp"

            # 1. Run our Python config implementation
            config = EasyLinClustConfig(
                fasta_files=fasta_file,
                cluster_prefix=func_output,
                tmp_dir=func_tmp,
                min_seq_id=0.3,
            )
            config.run()

            # 2. Run mmseqs CLI command
            cli_output.parent.mkdir(parents=True, exist_ok=True)
            cli_tmp.mkdir(parents=True, exist_ok=True)

            subprocess.run(
                [
                    "mmseqs",
                    "easy-linclust",
                    str(fasta_file),
                    str(cli_output),
                    str(cli_tmp),
                    "--min-seq-id", "0.3",
                ],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE
            )

            # 3. Compare outputs
            func_files = set(func_output.parent.glob("result*"))
            cli_files = set(cli_output.parent.glob("result*"))

            self.assertEqual(
                len(func_files),
                len(cli_files),
                "Different number of output files generated"
            )

            for func_file in func_files:
                filename = func_file.name
                cli_file = cli_output.parent / filename

                self.assertTrue(
                    cli_file.exists(),
                    f"File {filename} missing in CLI output"
                )

                with open(func_file, "rb") as f1, open(cli_file, "rb") as f2:
                    self.assertEqual(
                        f1.read(),
                        f2.read(),
                        f"Content mismatch in {filename}"
                    )


class TestBaseClusterParser(unittest.TestCase):
    """Unit tests for BaseClusterParser bug fixes (no mmseqs2 binary needed)."""

    def _make_parser(self, cluster_prefix=""):
        """Create a parser with mocked to_rep_list for testing split logic."""

        class MockConfig:
            pass

        MockConfig.cluster_prefix = cluster_prefix
        parser = BaseClusterParser(MockConfig())

        # Provide 100 items for splitting
        items = [(f"seq{i}", "ACGT") for i in range(100)]

        def mock_rep_list(with_seq=True):
            if with_seq:
                return list(items)
            return [s[0] for s in items]

        parser.to_rep_list = mock_rep_list
        return parser

    def test_split_normalization(self):
        """Proportions that don't sum to 1.0 are normalized correctly."""
        parser = self._make_parser()
        train, val, test = parser.split_rep_as_list(
            train=8, val=1, test=1, seed=42
        )
        total = len(train) + len(val) + len(test)
        self.assertEqual(total, 100, "All items must be accounted for")
        # 8:1:1 ratio → ~80/10/10
        self.assertTrue(
            75 <= len(train) <= 85,
            f"Train size {len(train)} not in expected range [75,85]"
        )

    def test_split_rep_as_fasta_no_unbound_error(self):
        """split_rep_as_fasta returns None for empty splits instead of raising."""
        with tempfile.TemporaryDirectory() as td:
            parser = self._make_parser(cluster_prefix=os.path.join(td, "result"))
            train_path, val_path, test_path = parser.split_rep_as_fasta(
                train=1.0, val=0, test=0
            )
            self.assertIsNotNone(train_path)
            self.assertIsNone(val_path)
            self.assertIsNone(test_path)

    def test_to_path_uses_correct_filename(self):
        """to_path returns _rep_seq.fasta (not _rep_seqs.fasta)."""
        parser = self._make_parser(cluster_prefix="/tmp/result")
        paths = parser.to_path()
        rep_seq_path = paths[2]
        self.assertTrue(
            rep_seq_path.endswith("_rep_seq.fasta"),
            f"Expected _rep_seq.fasta, got {rep_seq_path}"
        )


if __name__ == "__main__":
    unittest.main()
