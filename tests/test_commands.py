import os
import unittest
import subprocess
from pathlib import Path
import tempfile

from pymmseqs.config import CreateDBConfig, EasyLinClustConfig, EasyRbhConfig, EasyTaxonomyConfig
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


class TestEasyRbh(unittest.TestCase):
    """Test easy_rbh CLI argument generation."""

    def test_easy_rbh_args_match_cli(self):
        """Verify generated CLI args match expected mmseqs command."""
        config = EasyRbhConfig(
            query_fasta="/tmp/test_query.fasta",
            target_fasta_or_db="/tmp/test_target.fasta",
            alignment_file="/tmp/test_result",
            tmp_dir="/tmp/test_tmp",
            # Non-default params to verify they appear
            s=7.5,
            e=0.01,
            realign=True,
        )

        # Skip file existence checks for this test
        config._check_required_files = lambda: None
        config._caller_dir = Path("/tmp")

        args = config._get_command_args("easy_rbh")

        # Verify command name (underscores -> hyphens)
        self.assertEqual(args[0], "easy-rbh")

        # Verify required positional args in order
        self.assertEqual(args[1], "/tmp/test_query.fasta")
        self.assertEqual(args[2], "/tmp/test_target.fasta")
        self.assertEqual(args[3], "/tmp/test_result")
        self.assertEqual(args[4], "/tmp/test_tmp")

        # Verify non-default sensitivity appears
        self.assertIn("-s", args)
        idx = args.index("-s")
        self.assertEqual(args[idx + 1], "7.5")

        # Verify non-default e-value appears
        self.assertIn("-e", args)
        idx = args.index("-e")
        self.assertEqual(args[idx + 1], "0.01")

        # Verify non-default bool (realign) appears
        self.assertIn("--realign", args)
        idx = args.index("--realign")
        self.assertEqual(args[idx + 1], "1")

        # Verify default params are NOT in args
        self.assertNotIn("--min-seq-id", args)
        self.assertNotIn("--alignment-mode", args)
        self.assertNotIn("--cov-mode", args)

    def test_easy_rbh_default_createdb_mode(self):
        """Verify createdb_mode default is 1 (different from easy_search's 0)."""
        config = EasyRbhConfig(
            query_fasta="/tmp/test_query.fasta",
            target_fasta_or_db="/tmp/test_target.fasta",
            alignment_file="/tmp/test_result",
            tmp_dir="/tmp/test_tmp",
        )

        config._check_required_files = lambda: None
        config._caller_dir = Path("/tmp")

        args = config._get_command_args("easy_rbh")

        # createdb_mode=1 is the default for easy_rbh, so it should NOT appear
        self.assertNotIn("--createdb-mode", args)


class TestEasyTaxonomy(unittest.TestCase):
    """Test easy_taxonomy CLI argument generation."""

    def test_easy_taxonomy_args_match_cli(self):
        """Verify generated CLI args match expected mmseqs command."""
        config = EasyTaxonomyConfig(
            fasta_file="/tmp/test_query.fasta",
            target_db="/tmp/test_targetDB",
            tax_reports="/tmp/test_result",
            tmp_dir="/tmp/test_tmp",
            # Non-default params to verify they appear
            lca_mode=4,
            s=7.5,
            e=0.01,
        )

        # Skip file existence checks for this test
        config._check_required_files = lambda: None
        config._caller_dir = Path("/tmp")

        args = config._get_command_args("easy_taxonomy")

        # Verify command name (underscores -> hyphens)
        self.assertEqual(args[0], "easy-taxonomy")

        # Verify required positional args in order
        self.assertEqual(args[1], "/tmp/test_query.fasta")
        self.assertEqual(args[2], "/tmp/test_targetDB")
        self.assertEqual(args[3], "/tmp/test_result")
        self.assertEqual(args[4], "/tmp/test_tmp")

        # Verify non-default lca_mode appears
        self.assertIn("--lca-mode", args)
        idx = args.index("--lca-mode")
        self.assertEqual(args[idx + 1], "4")

        # Verify non-default sensitivity appears
        self.assertIn("-s", args)
        idx = args.index("-s")
        self.assertEqual(args[idx + 1], "7.5")

        # Verify non-default e-value appears
        self.assertIn("-e", args)
        idx = args.index("-e")
        self.assertEqual(args[idx + 1], "0.01")

        # Verify default params are NOT in args
        self.assertNotIn("--report-mode", args)
        self.assertNotIn("--vote-mode", args)
        self.assertNotIn("--min-seq-id", args)

    def test_easy_taxonomy_variadic_fasta(self):
        """Verify multiple FASTA files are handled correctly."""
        config = EasyTaxonomyConfig(
            fasta_file=["/tmp/file1.fasta", "/tmp/file2.fasta"],
            target_db="/tmp/test_targetDB",
            tax_reports="/tmp/test_result",
            tmp_dir="/tmp/test_tmp",
        )

        config._check_required_files = lambda: None
        config._caller_dir = Path("/tmp")

        args = config._get_command_args("easy_taxonomy")

        self.assertEqual(args[0], "easy-taxonomy")
        self.assertEqual(args[1], "/tmp/file1.fasta")
        self.assertEqual(args[2], "/tmp/file2.fasta")
        self.assertEqual(args[3], "/tmp/test_targetDB")
        self.assertEqual(args[4], "/tmp/test_result")
        self.assertEqual(args[5], "/tmp/test_tmp")


class TestEasyTaxonomyParser(unittest.TestCase):
    """Test EasyTaxonomyParser methods."""

    def test_parser_import(self):
        """Verify parser can be imported."""
        from pymmseqs.parsers.easy_taxonomy_parser import EasyTaxonomyParser
        self.assertTrue(EasyTaxonomyParser is not None)

    def test_parser_has_all_methods(self):
        """Verify parser class has all expected methods."""
        from pymmseqs.parsers.easy_taxonomy_parser import EasyTaxonomyParser
        expected_methods = [
            'to_pandas', 'to_list', 'to_gen', 'to_path',
            'to_json', 'to_csv',
            'report', 'lca_assignments', 'top_hits',
            'summary', 'composition', 'diversity',
            'rank_summary', 'unclassified_report', 'filter_by_taxon',
            'plot_composition', 'plot_rank_resolution',
            'plot_composition_pie', 'plot_classified_vs_unclassified',
            'plot_diversity_comparison',
        ]
        for method in expected_methods:
            self.assertTrue(
                hasattr(EasyTaxonomyParser, method),
                f"Missing method: {method}"
            )

    def test_parser_init_extracts_config_attrs(self):
        """Verify parser extracts correct attributes from config."""
        from pymmseqs.parsers.easy_taxonomy_parser import EasyTaxonomyParser

        config = EasyTaxonomyConfig(
            fasta_file="/tmp/test.fasta",
            target_db="/tmp/test_targetDB",
            tax_reports="/tmp/test_result",
            tmp_dir="/tmp/test_tmp",
            lca_ranks="superkingdom,phylum",
            tax_lineage=1,
        )
        config._check_required_files = lambda: None

        parser = EasyTaxonomyParser(config)
        self.assertEqual(parser._prefix, "/tmp/test_result")
        self.assertEqual(parser._lca_ranks, "superkingdom,phylum")
        self.assertEqual(parser._tax_lineage, 1)

    def test_to_path_returns_all_outputs(self):
        """Verify to_path returns dict with all 4 output files."""
        from pymmseqs.parsers.easy_taxonomy_parser import EasyTaxonomyParser

        config = EasyTaxonomyConfig(
            fasta_file="/tmp/test.fasta",
            target_db="/tmp/test_targetDB",
            tax_reports="/tmp/result",
            tmp_dir="/tmp/test_tmp",
        )
        config._check_required_files = lambda: None

        parser = EasyTaxonomyParser(config)
        paths = parser.to_path()

        self.assertIsInstance(paths, dict)
        self.assertIn('lca', paths)
        self.assertIn('report', paths)
        self.assertIn('tophit_aln', paths)
        self.assertIn('tophit_report', paths)
        self.assertEqual(paths['lca'], '/tmp/result_lca.tsv')
        self.assertEqual(paths['report'], '/tmp/result_report')

    def test_rank_constants(self):
        """Verify rank order and aliases are consistent."""
        from pymmseqs.parsers.easy_taxonomy_parser import (
            RANK_ORDER, RANK_ALIASES
        )
        # Core ranks must be in RANK_ORDER
        for rank in ['phylum', 'class', 'order', 'family', 'genus', 'species']:
            self.assertIn(rank, RANK_ORDER)
        # Aliases must map to ranks in RANK_ORDER
        for alias, canonical in RANK_ALIASES.items():
            self.assertTrue(
                alias in RANK_ORDER or canonical in RANK_ORDER,
                f"Alias {alias}->{canonical} not in RANK_ORDER"
            )


if __name__ == "__main__":
    unittest.main()
