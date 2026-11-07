"""Regression checks for schema handling, missing data, and bounded reductions.

Run from the repository root: python3 -m unittest discover -s test -p 'test_gviz.py'
"""

import gzip
import math
from pathlib import Path
import tempfile
import unittest

import numpy as np

from gviz.aggregate import IntervalBins, WindowBins, summarize
from gviz.data import Selection, Table, detect_kind


class BinningTests(unittest.TestCase):
    def test_window_gaps_and_missing_estimates_are_not_zero_distances(self):
        bins = WindowBins(np.array([0, 10, 20, 30, 40]))
        bins.add(1, 10, 0.2)
        bins.add(1, 10, 0.4)
        bins.add(11, 20, math.nan)
        result = bins.finish()
        np.testing.assert_allclose(result['mean'][0], 0.3)
        self.assertTrue(np.isnan(result['mean'][1:]).all())
        np.testing.assert_allclose(result['missing'][:2], [0, 1])
        self.assertTrue(np.isnan(result['missing'][2:]).all())

    def test_intervals_respect_inclusive_boundaries_and_overlaps(self):
        bins = IntervalBins(np.array([0, 10, 20, 30]))
        bins.add(1, 10, 0.2)
        bins.add(6, 25, 0.4)
        bins.add(21, 30, math.nan)
        result = bins.finish()
        np.testing.assert_array_equal(result['count'], [2, 1, 2])
        np.testing.assert_allclose(result['mean'], [(10 * 0.2 + 5 * 0.4) / 15, 0.4, 0.4])
        np.testing.assert_allclose(result['overlap'], [1.5, 1, 1.5])

    def test_adjacent_calls_in_one_bin_do_not_invent_an_overlap(self):
        bins = IntervalBins(np.array([0, 10]))
        bins.add(1, 5, 0.2)
        bins.add(6, 10, math.nan)
        result = bins.finish()
        self.assertEqual(result['count'][0], 2)
        self.assertEqual(result['overlap'][0], 1)

    def test_interval_lanes_respect_inclusive_endpoints_and_bound_crowding(self):
        from gviz.render import pack_intervals
        packed = pack_intervals([(1, 10, 0.1), (10, 15, 0.2), (11, 20, 0.3)])
        self.assertEqual([row[3] for row in packed], [0, 1, 0])
        self.assertIsNone(pack_intervals([(1, 10, 0.1)] * 13))

    def test_display_names_do_not_merge_distinct_inputs(self):
        from gviz.render import display_names
        names = display_names(['ref.fna', 'ref.fasta.gz', 'other_contigs.fasta.gz'])
        self.assertEqual(names['ref.fna'], 'ref.fna')
        self.assertEqual(names['ref.fasta.gz'], 'ref.fasta.gz')
        self.assertEqual(names['other_contigs.fasta.gz'], 'other')

    def test_chromosome_long_intervals_use_fixed_arrays_and_clip_to_region(self):
        bins = IntervalBins(np.linspace(100_000_000, 300_000_000, 1001))
        bins.add(1, 3_000_000_000, 0.3)
        result = bins.finish()
        self.assertEqual(len(result['mean']), 1000)
        np.testing.assert_allclose(result['mean'], 0.3)
        np.testing.assert_array_equal(result['count'], 1)
        np.testing.assert_allclose(result['overlap'], 1)

    def test_sparse_windows_keep_their_positions_on_long_chromosomes(self):
        bins = WindowBins(np.array([0, 3_000_000_000]))
        bins.add(1, 100, 0.2)
        self.assertEqual(bins.finish()['positions'][0], 50)

    def test_missing_interval_after_finite_interval_stays_missing(self):
        bins = IntervalBins(np.linspace(0, 101, 100))
        bins.add(1, 40, 0.3)
        bins.add(80, 101, math.nan)
        result = bins.finish()
        self.assertTrue(np.isnan(result['mean'][-20:]).all())
        self.assertTrue((result['count'][-20:] == 1).all())

    def test_interval_reduction_matches_direct_overlap_weights(self):
        edges = np.linspace(13, 417, 39)
        bins = IntervalBins(edges)
        rng = np.random.default_rng(5)
        expected_count = np.zeros(38)
        expected_total = np.zeros(38)
        expected_weight = np.zeros(38)
        expected_coverage = np.zeros(38)
        for _ in range(100):
            start = int(rng.integers(1, 500))
            end = start + int(rng.integers(1, 100))
            value = float(rng.random()) if start % 3 else math.nan
            bins.add(start, end, value)
            overlap = np.maximum(0, np.minimum(edges[1:], end) - np.maximum(edges[:-1], start - 1))
            expected_count += overlap > 0
            expected_coverage += overlap
            if math.isfinite(value):
                expected_weight += overlap
                expected_total += value * overlap
        result = bins.finish()
        np.testing.assert_array_equal(result['count'], expected_count)
        np.testing.assert_allclose(result['mean'], expected_total / expected_weight)
        np.testing.assert_allclose(result['overlap'], expected_coverage / np.diff(edges))


class TableTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def table(self, text, zipped=False):
        path = self.root / ('input.tsv.gz' if zipped else 'input.tsv')
        if zipped:
            with gzip.open(path, 'wt') as handle:
                handle.write(text)
        else:
            path.write_text(text)
        return Table(path)

    def test_all_current_formats_are_distinguished(self):
        base = ['seq', 'start', 'end', 'reference', 'd']
        self.assertEqual(detect_kind(base), 'roll')
        self.assertEqual(detect_kind(base + ['alpha', 'side']), 'map')
        self.assertEqual(detect_kind(['seq', 'start', 'end', 'reference', 'd_fw', 'd_rc']), 'roll')
        self.assertEqual(detect_kind(['genome_a', 'genome_b', 'd']), 'dist')
        self.assertEqual(detect_kind(['genome_a', 'genome_b', 'd', 'seq', 'start', 'end', 'direction']), 'samples')
        with self.assertRaisesRegex(ValueError, 'Sampled-window output needs'):
            detect_kind(['genome_a', 'genome_b', 'd', 'direction'])

    def test_gzip_comments_and_exact_sequence_selection(self):
        table = self.table('# invocation: gdiff roll\nseq\tstart\tend\treference\td\n'
                           'chr1\t1\t10\tref\t0.2\nchr2\t1\t100\tref\t0.4\n', True)
        report = summarize(table, Selection(), 10)
        self.assertEqual(report.title, 'chr2')
        self.assertEqual(report.rows, 1)
        self.assertEqual(table.provenance, ['invocation: gdiff roll'])
        report = summarize(table, Selection('chr1', (1, 20)), 10)
        self.assertEqual(report.title, 'chr1')
        self.assertEqual(report.data['edges'][-1], 20)

    def test_both_roll_strands_remain_separate(self):
        table = self.table('seq\tstart\tend\treference\td_fw\td_rc\nchr1\t1\t10\tref\t0.1\t0.3\n')
        report = summarize(table, Selection(), 1)
        self.assertEqual(len(report.data['series']), 2)
        self.assertAlmostEqual(report.data['series'][('ref', 'forward')]['mean'][0], 0.1)
        self.assertAlmostEqual(report.data['series'][('ref', 'reverse complement')]['mean'][0], 0.3)

    def test_map_report_draws_actual_intervals_and_distinguishes_missing_distances(self):
        from gviz.render import figure, plt, GRAY
        from matplotlib.colors import to_rgba
        table = self.table('seq\tseq_len\tstart\tend\treference\td\tside\n'
                           'chr1\t1000\t1\t1000\tref\t0.2\tlow\n'
                           'chr1\t1000\t201\t800\tref\t0.1\tlow\n'
                           'chr1\t1000\t401\t600\tref\tNA\tlow\n')
        report = summarize(table, Selection(), 10)
        fig = figure(report)
        try:
            segments = fig.axes[0].collections[0].get_segments()
            self.assertEqual(len(segments), 3)
            np.testing.assert_allclose(segments[0][:, 0], [0, 1])
            np.testing.assert_allclose(segments[2][:, 0], [0.4, 0.6])
            np.testing.assert_allclose(fig.axes[0].collections[0].get_colors()[2], to_rgba(GRAY))
            self.assertFalse(any('Peak' in text.get_text() for text in fig.axes[0].texts))
        finally:
            plt.close(fig)

    def test_binning_does_not_clip_short_distance_peaks_out_of_the_figure(self):
        from gviz.render import figure, plt
        table = self.table('seq\tstart\tend\treference\td\n' + ''.join(
            f'chr1\t{i * 10 + 1}\t{i * 10 + 10}\tref\t{0.45 if i == 10 else 0.01}\n'
            for i in range(20)))
        fig = figure(summarize(table, Selection(), 1))
        try:
            self.assertGreaterEqual(fig.axes[0].get_ylim()[1], 0.45)
        finally:
            plt.close(fig)

    def test_map_strands_are_not_mixed(self):
        table = self.table('seq\tstart\tend\treference\td\tside\tis_rc\n'
                           'chr1\t1\t10\tref\t0.1\tlow\t0\n'
                           'chr1\t1\t10\tref\t0.3\tlow\t1\n')
        report = summarize(table, Selection(), 1)
        self.assertEqual(len(report.data['series']), 2)
        self.assertAlmostEqual(report.data['series'][('ref', 'low', 'is_rc=0')]['mean'][0], 0.1)

    def test_many_references_are_bounded_and_filterable(self):
        from gviz.render import figure, plt
        import warnings
        table = self.table('seq\tstart\tend\treference\td\n' + ''.join(
            f'chr1\t1\t10\tref{i}\t0.2\n' for i in range(20)))
        report = summarize(table, Selection(), 10)
        self.assertEqual(len(report.data['series']), 12)
        self.assertTrue(any('12 series' in note for note in report.notes))
        fig = figure(report)
        try:
            with warnings.catch_warnings():
                warnings.simplefilter('error', UserWarning)
                fig.canvas.draw()
        finally:
            plt.close(fig)
        report = summarize(table, Selection(references=('ref19',)), 10)
        self.assertEqual(list(report.data['series']), [('ref19', '')])

    def test_sample_directions_are_separate_and_missing_values_excluded(self):
        table = self.table('genome_a\tgenome_b\tseq\tstart\tend\tdirection\td\n'
                           'a\tb\tchr1\t1\t10\tab\t0.1\n'
                           'a\tb\tchr1\t11\t20\tba\t0.2\n'
                           'a\tb\tchr1\t21\t30\tab\tNA\n')
        report = summarize(table, Selection(), 10)
        self.assertEqual(report.data['histograms'][('a', 'b', 'gdiff')].sum(), 1)
        self.assertEqual(report.data['histograms'][('b', 'a', 'gdiff')].sum(), 1)
        self.assertEqual(report.data['missing'][('a', 'b', 'gdiff')], 1)

    def test_all_missing_samples_show_one_legend_entry_and_no_curve(self):
        from gviz.render import figure, plt
        table = self.table('genome_a\tgenome_b\tseq\tstart\tend\tdirection\td\n'
                           'a\tb\tchr1\t1\t10\tab\tNA\n')
        fig = figure(summarize(table, Selection()))
        try:
            self.assertEqual(len(fig.axes[0].get_legend().get_texts()), 1)
            self.assertTrue(any(text.get_text() == 'No distance estimates' for text in fig.axes[0].texts))
        finally:
            plt.close(fig)

    def test_absent_matrix_pairs_remain_missing_and_repeats_are_averaged(self):
        table = self.table('genome_a\tgenome_b\td\na\tb\t0.2\na\tb\t0.4\na\tc\tNA\n')
        report = summarize(table, Selection())
        matrix = report.data['matrix']
        self.assertEqual(matrix.shape, (1, 2))
        self.assertAlmostEqual(matrix[0, 0], 0.3)
        self.assertTrue(np.isnan(matrix[0, 1]))
        self.assertEqual(report.data['y_names'], ['a'])
        self.assertEqual(report.data['x_names'], ['b', 'c'])

    def test_complete_all_pairs_input_keeps_a_symmetric_matrix(self):
        table = self.table('genome_a\tgenome_b\td\na\tb\t0.2\na\tc\t0.3\nb\tc\t0.4\n')
        report = summarize(table, Selection())
        matrix = report.data['matrix']
        self.assertTrue(report.data['symmetric'])
        self.assertEqual(matrix.shape, (3, 3))
        np.testing.assert_allclose(matrix, matrix.T)
        self.assertTrue(np.isnan(np.diag(matrix)).all())

    def test_dense_map_tracks_keep_all_calls_in_a_bounded_summary(self):
        from gviz.render import figure, plt
        table = self.table('seq\tstart\tend\treference\td\tside\n' + ''.join(
            'chr1\t1\t100\tref\t0.2\tlow\n' for _ in range(300)))
        report = summarize(table, Selection(), 10)
        values = report.data['series'][('ref', 'low')]
        self.assertIsNone(values['intervals'])
        self.assertEqual(values['rows'], 300)
        np.testing.assert_allclose(values['overlap'], 300)
        fig = figure(report)
        try:
            self.assertEqual(fig.axes[0].get_ylabel(), 'Mean overlap')
        finally:
            plt.close(fig)

    def test_roll_panels_separate_references_on_the_same_distance_scale(self):
        from gviz.render import figure, plt
        table = self.table('seq\tstart\tend\treference\td\n'
                           'chr1\t1\t10\ta\t0.1\nchr1\t11\t20\ta\tNA\n'
                           'chr1\t1\t10\tb\t0.3\nchr1\t11\t20\tb\t0.2\n')
        report = summarize(table, Selection(), 1000)
        fig = figure(report)
        try:
            self.assertEqual(len(fig.axes), 2)
            self.assertEqual(fig.axes[0].get_ylim(), fig.axes[1].get_ylim())
            self.assertTrue(np.isnan(fig.axes[0].lines[0].get_ydata()[-1]))
        finally:
            plt.close(fig)

    def test_large_genome_sets_do_not_allocate_quadratic_matrices(self):
        from gviz.render import figure, plt
        table = self.table('genome_a\tgenome_b\td\n' + ''.join(f'q\tr{i}\t0.2\n' for i in range(100)))
        report = summarize(table, Selection())
        self.assertIsNone(report.data['matrix'])
        self.assertEqual(report.data['histogram'].sum(), 100)
        fig = figure(report)
        try:
            fig.canvas.draw()
            self.assertEqual(len(fig.axes[1].patches), 20)
            self.assertEqual(len(fig.axes[1].get_yticklabels()), 20)
        finally:
            plt.close(fig)

    def test_map_region_clips_coordinates_without_losing_overlapping_calls(self):
        table = self.table('seq\tseq_len\tstart\tend\treference\td\tside\n'
                           'chr1\t1000\t1\t1000\tref\t0.2\tlow\n')
        report = summarize(table, Selection('chr1', (101, 200)), 10)
        self.assertEqual(report.kind, 'map')
        np.testing.assert_allclose(report.data['series'][('ref', 'low')]['mean'], 0.2)
        np.testing.assert_array_equal(report.data['series'][('ref', 'low')]['count'], 1)

    def test_malformed_records_have_line_numbers(self):
        table = self.table('seq\tstart\tend\treference\td\nchr1\t1\t10\tref\twrong\n')
        with self.assertRaisesRegex(ValueError, 'input.tsv:2:'):
            list(table.rows())
        with self.assertRaisesRegex(ValueError, 'No rows match'):
            summarize(self.table('seq\tstart\tend\treference\td\nchr1\t1\t10\tref\t0.2\n'), Selection('absent'))

    def test_empty_and_duplicate_headers_are_errors(self):
        with self.assertRaisesRegex(ValueError, 'no TSV header'):
            self.table('# comments only\n')
        with self.assertRaisesRegex(ValueError, 'Duplicate column'):
            self.table('seq\tseq\tstart\tend\treference\td\n')


if __name__ == '__main__':
    unittest.main()
