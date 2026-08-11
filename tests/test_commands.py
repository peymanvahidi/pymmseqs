import os
import unittest
import subprocess
from pathlib import Path
import tempfile

from pymmseqs.config import (
    ConvertAlisConfig,
    Convert2FastaConfig,
    CreateDBConfig,
    EasyClusterConfig,
    EasyLinClustConfig,
    EasyLinSearchConfig,
    EasyRbhConfig,
    EasyTaxonomyConfig,
    ExtractOrfsConfig,
)
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

    def test_createdb_v18_args(self):
        """Verify v18 params (gpu, createdb_mode=2, mask family) generate correct CLI args."""
        config = CreateDBConfig(
            fasta_file="/tmp/test_input.fasta",
            sequence_db="/tmp/test_db",
            createdb_mode=2,
            gpu=1,
            mask=0,
        )
        # Skip file-existence checks for arg-generation test
        config._check_required_files = lambda: None

        args = config._get_command_args("createdb")

        self.assertEqual(args[0], "createdb")

        # Non-default v18 params appear as flags
        self.assertIn("--gpu", args)
        self.assertEqual(args[args.index("--gpu") + 1], "1")

        self.assertIn("--createdb-mode", args)
        self.assertEqual(args[args.index("--createdb-mode") + 1], "2")

        self.assertIn("--mask", args)
        self.assertEqual(args[args.index("--mask") + 1], "0")

        # Params left at default are omitted
        self.assertNotIn("--mask-prob", args)
        self.assertNotIn("--mask-n-repeat", args)


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

    def test_easy_linclust_v18_args(self):
        """Verify v18 resync params (gpu + previously-dropped sub_mat/max_seq_len/db_load_mode) emit."""
        config = EasyLinClustConfig(
            fasta_files="/tmp/in.fasta",
            cluster_prefix="/tmp/result",
            tmp_dir="/tmp/tmp",
            gpu=1,
            max_seq_len=1000,
            db_load_mode=2,
        )
        config._check_required_files = lambda: None
        config._caller_dir = Path("/tmp")
        args = config._get_command_args("easy-linclust")

        self.assertEqual(args[0], "easy-linclust")
        self.assertIn("--gpu", args)
        self.assertEqual(args[args.index("--gpu") + 1], "1")
        self.assertIn("--max-seq-len", args)
        self.assertEqual(args[args.index("--max-seq-len") + 1], "1000")
        self.assertIn("--db-load-mode", args)
        self.assertEqual(args[args.index("--db-load-mode") + 1], "2")

    def test_easy_linclust_defaults_omit_v18_args(self):
        """At defaults, the resynced params must NOT be emitted (byte-for-byte safety)."""
        config = EasyLinClustConfig(
            fasta_files="/tmp/in.fasta",
            cluster_prefix="/tmp/result",
            tmp_dir="/tmp/tmp",
        )
        config._check_required_files = lambda: None
        config._caller_dir = Path("/tmp")
        args = config._get_command_args("easy-linclust")
        for flag in ("--gpu", "--sub-mat", "--max-seq-len", "--db-load-mode"):
            self.assertNotIn(flag, args)


class TestEasyCluster(unittest.TestCase):
    """Test easy_cluster CLI argument generation (v18 --gpu resync)."""

    def test_easy_cluster_v18_args(self):
        """Verify --gpu and --createdb-mode emit when set non-default."""
        config = EasyClusterConfig(
            fasta_files="/tmp/in.fasta",
            cluster_prefix="/tmp/result",
            tmp_dir="/tmp/tmp",
            gpu=1,
            createdb_mode=2,
        )
        config._check_required_files = lambda: None
        config._caller_dir = Path("/tmp")
        args = config._get_command_args("easy-cluster")

        self.assertEqual(args[0], "easy-cluster")
        self.assertIn("--gpu", args)
        self.assertEqual(args[args.index("--gpu") + 1], "1")
        self.assertIn("--createdb-mode", args)
        self.assertEqual(args[args.index("--createdb-mode") + 1], "2")

    def test_easy_cluster_default_gpu_omitted(self):
        """gpu defaults to 0 and must be omitted."""
        config = EasyClusterConfig(
            fasta_files="/tmp/in.fasta",
            cluster_prefix="/tmp/result",
            tmp_dir="/tmp/tmp",
        )
        config._check_required_files = lambda: None
        config._caller_dir = Path("/tmp")
        args = config._get_command_args("easy-cluster")
        self.assertNotIn("--gpu", args)


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

    def test_split_matches_sklearn_exactly(self):
        """
        The split is pinned to what sklearn's train_test_split produced.

        scikit-learn was dropped as a dependency (it ships no musllinux wheel,
        which broke `pip install` on Alpine) and replaced with a numpy
        reimplementation. These values were captured from the sklearn version,
        so users with a fixed seed keep getting the same split.
        """
        from pymmseqs.parsers.base_cluster_parser import _train_test_split

        items = [f"seq{i}" for i in range(100)]
        train, test = _train_test_split(items, test_size=0.2, shuffle=True, random_state=42)
        self.assertEqual((len(train), len(test)), (80, 20))
        self.assertEqual(train[:5], ["seq55", "seq88", "seq26", "seq42", "seq69"])
        self.assertEqual(test[:5], ["seq83", "seq53", "seq70", "seq45", "seq44"])

    def test_split_test_size_float_edge_case(self):
        """
        n_train must be n - n_test, not floor(n * (1 - test_size)).

        At test_size=0.9 the latter yields 0 instead of 1, because
        1 - 0.9 == 0.09999999999999998. This silently disagreed with sklearn
        on 16 of 400 combinations before it was caught.
        """
        from pymmseqs.parsers.base_cluster_parser import _train_test_split

        train, test = _train_test_split(
            [f"s{i}" for i in range(10)], test_size=0.9, shuffle=True, random_state=0
        )
        self.assertEqual(len(train), 1, "one item must remain in train, not zero")
        self.assertEqual(len(test), 9)
        self.assertEqual(len(train) + len(test), 10, "no item may be dropped")

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


class TestEasyLinSearch(unittest.TestCase):
    """Test easy_linsearch CLI argument generation and parser integration."""

    def test_easy_linsearch_args_match_cli(self):
        """Verify generated CLI args match expected mmseqs command."""
        config = EasyLinSearchConfig(
            query_fasta="/tmp/test_query.fasta",
            target_fasta_or_db="/tmp/test_target.fasta",
            alignment_file="/tmp/test_result",
            tmp_dir="/tmp/test_tmp",
            # Non-default params to verify they appear
            format_mode=4,
            e=0.01,
            min_seq_id=0.5,
        )

        config._check_required_files = lambda: None
        config._caller_dir = Path("/tmp")

        args = config._get_command_args("easy_linsearch")

        # Verify command name (underscores -> hyphens)
        self.assertEqual(args[0], "easy-linsearch")

        # Verify required positional args in order
        self.assertEqual(args[1], "/tmp/test_query.fasta")
        self.assertEqual(args[2], "/tmp/test_target.fasta")
        self.assertEqual(args[3], "/tmp/test_result")
        self.assertEqual(args[4], "/tmp/test_tmp")

        # format_mode=4 is required by EasySearchParser
        self.assertIn("--format-mode", args)
        self.assertEqual(args[args.index("--format-mode") + 1], "4")

        self.assertIn("-e", args)
        self.assertEqual(args[args.index("-e") + 1], "0.01")

        self.assertIn("--min-seq-id", args)
        self.assertEqual(args[args.index("--min-seq-id") + 1], "0.5")

    def test_easy_linsearch_end_to_end(self):
        """Run easy_linsearch and verify it returns a working EasySearchParser.

        Note: easy-linsearch is the fast/less-sensitive linear search and may
        legitimately return zero hits for a small protein input, so this asserts
        the command + parser wiring (format_mode=4 header parsed, output produced)
        rather than the presence of a hit.
        """
        from pymmseqs.commands.easy_linsearch import easy_linsearch

        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            query = tmp_path / "query.fasta"
            target = tmp_path / "target.fasta"
            # Diverse ~127-residue protein fragment (avoids low-complexity masking)
            seq = (
                "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKAL"
                "PDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGE"
            )
            query.write_text(f">q1\n{seq}\n")
            target.write_text(f">t1\n{seq}\n")
            out = tmp_path / "aln.m8"

            parser = easy_linsearch(query, target, out, tmp_dir=tmp_path / "tmp")

            # Parser is wired to the alignment file and parses the format_mode=4 header
            self.assertTrue(Path(parser.to_path()).exists())
            df = parser.to_pandas()
            for col in ("query", "target", "fident", "evalue", "bits"):
                self.assertIn(col, df.columns)


class TestConvertAlis(unittest.TestCase):
    """Test convertalis CLI argument generation and parser integration."""

    def test_convertalis_args_match_cli(self):
        """Verify generated CLI args match expected mmseqs command."""
        config = ConvertAlisConfig(
            query_db="/tmp/qdb",
            target_db="/tmp/tdb",
            alignment_db="/tmp/alndb",
            alignment_file="/tmp/out.m8",
            # Non-default params to verify they appear
            format_mode=4,
            search_type=3,
        )

        config._check_required_files = lambda: None
        config._caller_dir = Path("/tmp")

        args = config._get_command_args("convertalis")

        self.assertEqual(args[0], "convertalis")

        # Verify required positional args in order
        self.assertEqual(args[1], "/tmp/qdb")
        self.assertEqual(args[2], "/tmp/tdb")
        self.assertEqual(args[3], "/tmp/alndb")
        self.assertEqual(args[4], "/tmp/out.m8")

        # format_mode=4 is required by EasySearchParser
        self.assertIn("--format-mode", args)
        self.assertEqual(args[args.index("--format-mode") + 1], "4")

        self.assertIn("--search-type", args)
        self.assertEqual(args[args.index("--search-type") + 1], "3")

    def test_convertalis_end_to_end(self):
        """createdb -> search -> convertalis should yield a parsable alignment table."""
        from pymmseqs.commands.convertalis import convertalis

        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            seq = (
                "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKAL"
                "PDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGE"
            )
            query = tmp_path / "q.fasta"
            target = tmp_path / "t.fasta"
            query.write_text(f">q1\n{seq}\n")
            target.write_text(f">t1\n{seq}\n")

            qdb = tmp_path / "qdb"
            tdb = tmp_path / "tdb"
            alndb = tmp_path / "alndb"
            mmtmp = tmp_path / "mmtmp"
            mmtmp.mkdir()

            # Build inputs with the mmseqs CLI (isolates the convertalis wrapper under test)
            for src, db in ((query, qdb), (target, tdb)):
                subprocess.run(
                    ["mmseqs", "createdb", str(src), str(db)],
                    check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                )
            subprocess.run(
                ["mmseqs", "search", str(qdb), str(tdb), str(alndb), str(mmtmp)],
                check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            )

            out = tmp_path / "out.m8"
            parser = convertalis(qdb, tdb, alndb, out)
            df = parser.to_pandas()

            self.assertIn("query", df.columns)
            self.assertIn("target", df.columns)
            self.assertGreaterEqual(len(df), 1)


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


class TestEasyRbhParser(unittest.TestCase):
    """Test EasyRbhParser methods."""

    def test_parser_import(self):
        """Verify parser can be imported."""
        from pymmseqs.parsers.easy_rbh_parser import EasyRbhParser
        self.assertTrue(EasyRbhParser is not None)

    def test_parser_has_all_methods(self):
        """Verify parser class has all expected methods."""
        from pymmseqs.parsers.easy_rbh_parser import EasyRbhParser
        expected_methods = [
            'to_pandas', 'to_list', 'to_gen', 'to_path',
            'to_json', 'to_csv', 'summary',
            'plot_identity_distribution', 'plot_evalue_distribution',
            'plot_alignment_length', 'plot_score_distribution',
        ]
        for method in expected_methods:
            self.assertTrue(
                hasattr(EasyRbhParser, method),
                f"Missing method: {method}"
            )

    def test_parser_init_extracts_config_attrs(self):
        """Verify parser extracts correct attributes from config."""
        from pymmseqs.parsers.easy_rbh_parser import EasyRbhParser

        config = EasyRbhConfig(
            query_fasta="/tmp/test_query.fasta",
            target_fasta_or_db="/tmp/test_target.fasta",
            alignment_file="/tmp/test_rbh_result.tsv",
            tmp_dir="/tmp/test_tmp",
        )
        config._check_required_files = lambda: None

        parser = EasyRbhParser(config)
        self.assertEqual(str(parser._alignment_file), "/tmp/test_rbh_result.tsv")
        self.assertEqual(parser._format_mode, 0)
        self.assertEqual(
            parser._format_output,
            "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits"
        )

    def test_to_path_returns_string(self):
        """Verify to_path returns alignment file path as string."""
        from pymmseqs.parsers.easy_rbh_parser import EasyRbhParser

        config = EasyRbhConfig(
            query_fasta="/tmp/test_query.fasta",
            target_fasta_or_db="/tmp/test_target.fasta",
            alignment_file="/tmp/test_rbh_result.tsv",
            tmp_dir="/tmp/test_tmp",
        )
        config._check_required_files = lambda: None

        parser = EasyRbhParser(config)
        path = parser.to_path()
        self.assertIsInstance(path, str)
        self.assertEqual(path, "/tmp/test_rbh_result.tsv")

    def test_get_columns_parses_format_output(self):
        """Verify _get_columns correctly parses format_output string."""
        from pymmseqs.parsers.easy_rbh_parser import EasyRbhParser

        config = EasyRbhConfig(
            query_fasta="/tmp/test_query.fasta",
            target_fasta_or_db="/tmp/test_target.fasta",
            alignment_file="/tmp/test_rbh_result.tsv",
            tmp_dir="/tmp/test_tmp",
            format_output="query,target,fident,evalue",
        )
        config._check_required_files = lambda: None

        parser = EasyRbhParser(config)
        columns = parser._get_columns()
        self.assertEqual(columns, ["query", "target", "fident", "evalue"])

    def test_parser_format_mode_4(self):
        """Verify parser stores format_mode=4 from config."""
        from pymmseqs.parsers.easy_rbh_parser import EasyRbhParser

        config = EasyRbhConfig(
            query_fasta="/tmp/test_query.fasta",
            target_fasta_or_db="/tmp/test_target.fasta",
            alignment_file="/tmp/test_rbh_result.tsv",
            tmp_dir="/tmp/test_tmp",
            format_mode=4,
        )
        config._check_required_files = lambda: None

        parser = EasyRbhParser(config)
        self.assertEqual(parser._format_mode, 4)

    def test_to_pandas_format_mode_0(self):
        """Verify to_pandas uses format_output columns for format_mode=0."""
        from pymmseqs.parsers.easy_rbh_parser import EasyRbhParser

        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            f.write("q1\tt1\t0.95\t100\t5\t0\t1\t100\t1\t100\t1e-50\t200\n")
            f.write("q2\tt2\t0.80\t80\t16\t0\t1\t80\t1\t80\t1e-30\t150\n")
            tmp_file = f.name

        try:
            config = EasyRbhConfig(
                query_fasta="/tmp/test_query.fasta",
                target_fasta_or_db="/tmp/test_target.fasta",
                alignment_file=tmp_file,
                tmp_dir="/tmp/test_tmp",
            )
            config._check_required_files = lambda: None
            parser = EasyRbhParser(config)
            df = parser.to_pandas()
            self.assertEqual(len(df), 2)
            self.assertIn("query", df.columns)
            self.assertIn("fident", df.columns)
            self.assertIn("evalue", df.columns)
        finally:
            os.unlink(tmp_file)

    def test_to_pandas_format_mode_4(self):
        """Verify to_pandas reads headers from file for format_mode=4."""
        from pymmseqs.parsers.easy_rbh_parser import EasyRbhParser

        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            f.write("query\ttarget\tfident\tevalue\n")
            f.write("q1\tt1\t0.95\t1e-50\n")
            tmp_file = f.name

        try:
            config = EasyRbhConfig(
                query_fasta="/tmp/test_query.fasta",
                target_fasta_or_db="/tmp/test_target.fasta",
                alignment_file=tmp_file,
                tmp_dir="/tmp/test_tmp",
                format_mode=4,
                format_output="query,target,fident,evalue",
            )
            config._check_required_files = lambda: None
            parser = EasyRbhParser(config)
            df = parser.to_pandas()
            self.assertEqual(len(df), 1)
            self.assertEqual(list(df.columns), ["query", "target", "fident", "evalue"])
        finally:
            os.unlink(tmp_file)

    def test_summary_with_default_columns(self):
        """Verify summary returns expected keys with default format_output."""
        from pymmseqs.parsers.easy_rbh_parser import EasyRbhParser

        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            f.write("q1\tt1\t0.95\t100\t5\t0\t1\t100\t1\t100\t1e-50\t200\n")
            f.write("q2\tt2\t0.80\t80\t16\t0\t1\t80\t1\t80\t1e-30\t150\n")
            tmp_file = f.name

        try:
            config = EasyRbhConfig(
                query_fasta="/tmp/test_query.fasta",
                target_fasta_or_db="/tmp/test_target.fasta",
                alignment_file=tmp_file,
                tmp_dir="/tmp/test_tmp",
            )
            config._check_required_files = lambda: None
            parser = EasyRbhParser(config)
            stats = parser.summary()
            self.assertEqual(stats["total_pairs"], 2)
            self.assertIn("mean_identity", stats)
            self.assertIn("median_evalue", stats)
            self.assertIn("mean_alignment_length", stats)
            self.assertIn("median_bit_score", stats)
        finally:
            os.unlink(tmp_file)

    def test_summary_with_minimal_columns(self):
        """Verify summary gracefully handles missing columns."""
        from pymmseqs.parsers.easy_rbh_parser import EasyRbhParser

        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            f.write("q1\tt1\n")
            f.write("q2\tt2\n")
            tmp_file = f.name

        try:
            config = EasyRbhConfig(
                query_fasta="/tmp/test_query.fasta",
                target_fasta_or_db="/tmp/test_target.fasta",
                alignment_file=tmp_file,
                tmp_dir="/tmp/test_tmp",
                format_output="query,target",
            )
            config._check_required_files = lambda: None
            parser = EasyRbhParser(config)
            stats = parser.summary()
            self.assertEqual(stats["total_pairs"], 2)
            self.assertNotIn("mean_identity", stats)
            self.assertNotIn("median_evalue", stats)
        finally:
            os.unlink(tmp_file)

    def test_to_list_returns_dicts(self):
        """Verify to_list returns list of dicts with correct keys."""
        from pymmseqs.parsers.easy_rbh_parser import EasyRbhParser

        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            f.write("q1\tt1\t0.95\t100\t5\t0\t1\t100\t1\t100\t1e-50\t200\n")
            tmp_file = f.name

        try:
            config = EasyRbhConfig(
                query_fasta="/tmp/test_query.fasta",
                target_fasta_or_db="/tmp/test_target.fasta",
                alignment_file=tmp_file,
                tmp_dir="/tmp/test_tmp",
            )
            config._check_required_files = lambda: None
            parser = EasyRbhParser(config)
            result = parser.to_list()
            self.assertIsInstance(result, list)
            self.assertEqual(len(result), 1)
            self.assertEqual(result[0]["query"], "q1")
            self.assertEqual(result[0]["target"], "t1")
        finally:
            os.unlink(tmp_file)

    def test_to_gen_yields_typed_values(self):
        """Verify to_gen yields dicts with correct types."""
        from pymmseqs.parsers.easy_rbh_parser import EasyRbhParser

        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            f.write("q1\tt1\t0.95\t100\t5\t0\t1\t100\t1\t100\t1E-50\t200\n")
            tmp_file = f.name

        try:
            config = EasyRbhConfig(
                query_fasta="/tmp/test_query.fasta",
                target_fasta_or_db="/tmp/test_target.fasta",
                alignment_file=tmp_file,
                tmp_dir="/tmp/test_tmp",
            )
            config._check_required_files = lambda: None
            parser = EasyRbhParser(config)
            rows = list(parser.to_gen())
            self.assertEqual(len(rows), 1)
            row = rows[0]
            self.assertIsInstance(row["fident"], float)
            self.assertIsInstance(row["alnlen"], int)
            self.assertIsInstance(row["evalue"], float)
            self.assertIsInstance(row["bits"], int)
        finally:
            os.unlink(tmp_file)


class TestConvert2Fasta(unittest.TestCase):
    """Test convert2fasta CLI argument generation and parser integration."""

    def test_convert2fasta_args_match_cli(self):
        """Verify generated CLI args match expected mmseqs command."""
        config = Convert2FastaConfig(
            sequence_db="/tmp/seqdb",
            fasta_file="/tmp/out.fasta",
            use_header_file=True,
        )
        config._check_required_files = lambda: None
        config._caller_dir = Path("/tmp")

        args = config._get_command_args("convert2fasta")

        self.assertEqual(args[0], "convert2fasta")
        self.assertEqual(args[1], "/tmp/seqdb")
        self.assertEqual(args[2], "/tmp/out.fasta")

        # Non-default bool appears as "1"
        self.assertIn("--use-header-file", args)
        self.assertEqual(args[args.index("--use-header-file") + 1], "1")

    def test_convert2fasta_end_to_end(self):
        """createdb -> convert2fasta should round-trip sequences back to FASTA."""
        from pymmseqs.commands.createdb import createdb
        from pymmseqs.commands.convert2fasta import convert2fasta

        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            fasta = tmp_path / "in.fasta"
            fasta.write_text(">s1\nACDEFGHIKL\n>s2\nMNPQRSTVWY\n")

            createdb(fasta, tmp_path / "seqdb")
            parser = convert2fasta(tmp_path / "seqdb", tmp_path / "out.fasta")

            df = parser.to_pandas()
            self.assertEqual(list(df.columns), ["header", "sequence"])
            self.assertEqual(len(parser), 2)
            self.assertEqual(set(df["sequence"]), {"ACDEFGHIKL", "MNPQRSTVWY"})

            # to_list / to_gen / to_path coverage
            self.assertEqual(len(parser.to_list()), 2)
            self.assertEqual(len(list(parser.to_gen())), 2)
            self.assertTrue(Path(parser.to_path()).exists())


class TestExtractOrfs(unittest.TestCase):
    """Test extractorfs CLI argument generation and parser integration."""

    # Contig with a clean forward ORF (ATG ... TAA)
    CONTIG = (
        "ATGGCAAAACGTTTAGCAGAAGAACTGGGCATTGAAGTGCAGGCACCGATTCTGAGCCGTGTT"
        "GGCGATGGCACCCAGGATAACCTGAGCGGCGCAGAAAAAGCAGTGCAGGTGAAAGTGAAAGCA"
        "CTGCCGGATGCACAGTTTGAAGTGTAA"
    )

    def test_extractorfs_args_match_cli(self):
        """Verify generated CLI args match expected mmseqs command."""
        config = ExtractOrfsConfig(
            sequence_db="/tmp/ntdb",
            orf_db="/tmp/orfdb",
            min_length=10,
            translate=True,
            orf_start_mode=0,
        )
        config._check_required_files = lambda: None
        config._caller_dir = Path("/tmp")

        args = config._get_command_args("extractorfs")

        self.assertEqual(args[0], "extractorfs")
        self.assertEqual(args[1], "/tmp/ntdb")
        self.assertEqual(args[2], "/tmp/orfdb")

        self.assertIn("--min-length", args)
        self.assertEqual(args[args.index("--min-length") + 1], "10")

        self.assertIn("--translate", args)
        self.assertEqual(args[args.index("--translate") + 1], "1")

        self.assertIn("--orf-start-mode", args)
        self.assertEqual(args[args.index("--orf-start-mode") + 1], "0")

        # Defaults are omitted
        self.assertNotIn("--max-length", args)
        self.assertNotIn("--contig-start-mode", args)

    def test_extractorfs_end_to_end(self):
        """createdb -> extractorfs should yield a coordinate-bearing ORF table."""
        from pymmseqs.commands.createdb import createdb
        from pymmseqs.commands.extractorfs import extractorfs

        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            nt = tmp_path / "nt.fasta"
            nt.write_text(f">contig1 example\n{self.CONTIG}\n")

            createdb(nt, tmp_path / "ntdb")
            parser = extractorfs(tmp_path / "ntdb", tmp_path / "orfdb", min_length=10)

            df = parser.to_pandas()
            for col in ("orf_id", "source_id", "source_name", "start", "end",
                        "strand", "length", "sequence"):
                self.assertIn(col, df.columns)

            self.assertGreaterEqual(len(parser), 1)
            # Source key resolved back to the original contig accession
            self.assertIn("contig1", set(df["source_name"]))
            # Strands are limited to +/-
            self.assertTrue(set(df["strand"].dropna()).issubset({"+", "-"}))

            # Coordinates are absolute positions on the contig: end>start on '+',
            # end<start on '-', and abs(end-start)+1 == length for every ORF.
            coord = df.dropna(subset=["start", "end", "strand"])
            for _, r in coord.iterrows():
                self.assertEqual(abs(int(r["end"]) - int(r["start"])) + 1, int(r["length"]))
                if r["strand"] == "+":
                    self.assertGreaterEqual(int(r["end"]), int(r["start"]))
                else:
                    self.assertLessEqual(int(r["end"]), int(r["start"]))

            # to_gen() and to_list() agree
            self.assertEqual(len(list(parser.to_gen())), len(parser.to_list()))

            stats = parser.summary()
            for key in ("total_orfs", "source_sequences", "forward_orfs",
                        "reverse_orfs", "mean_length"):
                self.assertIn(key, stats)
            self.assertEqual(stats["source_sequences"], 1)


if __name__ == "__main__":
    unittest.main()
