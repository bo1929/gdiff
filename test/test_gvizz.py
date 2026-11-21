"""Regression checks for schema handling, exact interval tracks, and bounded outputs.

Run from the repository root: python3 -m unittest discover -s test -p 'test_gvizz.py'
"""

import gzip
import math
from pathlib import Path
import tempfile
import unittest

import numpy as np

from gvizz import aggregate
from gvizz.aggregate import summarize
from gvizz.data import Selection, Table, detect_kind

LABELS = {'$\\hat{D}$', 'lower tail', 'upper tail'}


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
        self.assertEqual(detect_kind(base + ['pi', 'side']), 'map')
        self.assertEqual(detect_kind(['seq', 'start', 'end', 'reference', 'd_fw', 'd_rc']), 'roll')
        self.assertEqual(detect_kind(['genome_a', 'genome_b', 'd']), 'dist')
        self.assertEqual(detect_kind(['genome_a', 'genome_b', 'd', 'seq', 'start', 'end', 'direction']), 'samples')
        with self.assertRaisesRegex(ValueError, 'Sampled-window output needs'):
            detect_kind(['genome_a', 'genome_b', 'd', 'direction'])

    def test_gzip_comments_and_exact_sequence_selection(self):
        table = self.table('# invocation: gdiff roll\nseq\tstart\tend\treference\td\n'
                           'chr1\t1\t10\tref\t0.2\nchr2\t1\t100\tref\t0.4\n', True)
        report = summarize(table, Selection())
        self.assertEqual(report.title, '2 sequences')  # no implicit selection
        self.assertEqual([s['seq'] for s in report.data['sequences']], ['chr1', 'chr2'])
        self.assertEqual(report.rows, 2)
        self.assertEqual(table.provenance, ['invocation: gdiff roll'])
        report = summarize(table, Selection('chr1'))
        self.assertEqual(report.title, 'chr1')
        self.assertEqual(report.data['sequences'][0]['extent'], (1, 10))

    def test_every_sequence_is_traversable_without_selection(self):
        table = self.table('seq\tseq_len\tstart\tend\treference\td\tside\n'
                           'long\t100000\t1\t10\tref\t0.2\tlow\n'
                           'short\t1000\t100\t600\tref\t0.1\tlow\n')
        report = summarize(table, Selection())
        self.assertEqual([s['seq'] for s in report.data['sequences']], ['long', 'short'])
        self.assertTrue(any('arrow links' in note for note in report.notes))
        report = summarize(table, Selection('short'))
        self.assertEqual([s['seq'] for s in report.data['sequences']], ['short'])

    def test_both_roll_strands_remain_separate(self):
        table = self.table('seq\tstart\tend\treference\td_fw\td_rc\nchr1\t1\t10\tref\t0.1\t0.3\n')
        report = summarize(table, Selection())
        tracks = report.data['sequences'][0]['tracks']
        self.assertEqual(len(tracks), 2)
        self.assertEqual(tracks[('ref', 'forward')], [(1, 10, 0.1)])
        self.assertEqual(tracks[('ref', 'reverse complement')], [(1, 10, 0.3)])

    def test_map_tracks_keep_exact_intervals_and_gray_missing(self):
        from gvizz.render import figure, plt, GRAY
        from matplotlib.colors import to_rgba
        table = self.table('seq\tseq_len\tstart\tend\treference\td\tside\n'
                           'chr1\t1000\t1\t1000\tref\t0.2\tlow\n'
                           'chr1\t1000\t201\t800\tref\t0.1\tlow\n'
                           'chr1\t1000\t401\t600\tref\tNA\tlow\n')
        report = summarize(table, Selection())
        self.assertEqual(report.data['sequences'][0]['tracks'][('ref', 'low')],
                         [(1, 1000, 0.2), (201, 800, 0.1), (401, 600, math.nan)])
        fig = figure(report)
        try:
            collection = fig.axes[0].collections[0]
            paths = collection.get_paths()
            self.assertEqual(len(paths), 3)
            xs = [sorted(set(np.round(p.vertices[:, 0], 6))) for p in paths]
            np.testing.assert_allclose(xs[0], [0, 1000])     # highest d first
            np.testing.assert_allclose(xs[1], [200, 800])
            np.testing.assert_allclose(xs[2], [400, 600])    # missing underneath
            np.testing.assert_allclose(collection.get_facecolors()[2], to_rgba(GRAY))
            ys = paths[0].vertices[:, 1]
            self.assertEqual(sorted(set(np.round(ys, 6))), [-0.45, 0.45])  # one fixed level
            self.assertEqual(len(fig.axes), 2)  # track + colorbar
        finally:
            plt.close(fig)

    def test_lower_distance_is_drawn_on_top(self):
        from gvizz.render import draw_order
        ordered = draw_order([(1, 10, 0.5), (2, 8, 0.1), (3, 9, math.nan), (4, 7, 0.3)])
        self.assertEqual([seg[0] for seg in ordered], [1, 4, 2, 3])

    def test_map_strands_are_not_mixed(self):
        table = self.table('seq\tstart\tend\treference\td\tside\tis_rc\n'
                           'chr1\t1\t10\tref\t0.1\tlow\t0\n'
                           'chr1\t1\t10\tref\t0.3\tlow\t1\n')
        report = summarize(table, Selection())
        tracks = report.data['sequences'][0]['tracks']
        self.assertEqual(len(tracks), 2)
        self.assertEqual(tracks[('ref', 'low', 'is_rc=0')], [(1, 10, 0.1)])

    def test_many_references_are_bounded_and_filterable(self):
        from gvizz.render import figure, plt
        import warnings
        table = self.table('seq\tstart\tend\treference\td\n' + ''.join(
            f'chr1\t1\t10\tref{i}\t0.2\n' for i in range(20)))
        report = summarize(table, Selection())
        self.assertEqual(len(report.data['sequences'][0]['tracks']), 12)
        self.assertTrue(any('12 tracks' in note for note in report.notes))
        fig = figure(report)
        try:
            with warnings.catch_warnings():
                warnings.simplefilter('error', UserWarning)
                fig.canvas.draw()
        finally:
            plt.close(fig)

    def test_track_cap_keeps_the_longest_intervals(self):
        original = aggregate.MAX_SEGMENTS
        aggregate.MAX_SEGMENTS = 3
        try:
            table = self.table('seq\tstart\tend\treference\td\n' + ''.join(
                f'chr1\t1\t{10 * i}\tref\t0.2\n' for i in range(1, 6)))
            report = summarize(table, Selection())
            self.assertEqual(len(report.data['sequences'][0]['tracks'][('ref', '')]), 3)
            self.assertEqual(report.data['sequences'][0]['tracks'][('ref', '')][-1][1], 50)
            self.assertTrue(any('shortest intervals omitted' in note for note in report.notes))
        finally:
            aggregate.MAX_SEGMENTS = original

    def test_overlapping_intervals_share_one_level(self):
        from gvizz.render import figure, plt
        table = self.table('seq\tstart\tend\treference\td\tside\n' + ''.join(
            'chr1\t1\t100\tref\t0.2\tlow\n' for _ in range(40)))
        report = summarize(table, Selection())
        fig = figure(report)
        try:
            self.assertEqual(len(fig.axes[0].collections[0].get_paths()), 40)  # one level, all drawn
        finally:
            plt.close(fig)

    def test_sample_directions_are_separate_and_missing_values_excluded(self):
        table = self.table('genome_a\tgenome_b\tseq\tstart\tend\tdirection\td\n'
                           'a\tb\tchr1\t1\t10\tab\t0.1\n'
                           'a\tb\tchr1\t11\t20\tba\t0.2\n'
                           'a\tb\tchr1\t21\t30\tab\tNA\n')
        report = summarize(table, Selection())
        np.testing.assert_allclose(report.data['ecdfs'][('a', 'b', 'gdiff')], [0.1])
        np.testing.assert_allclose(report.data['ecdfs'][('b', 'a', 'gdiff')], [0.2])
        self.assertEqual(report.data['missing'][('a', 'b', 'gdiff')], 1)

    def test_all_missing_samples_show_one_empty_panel(self):
        from gvizz.render import figure, plt
        table = self.table('genome_a\tgenome_b\tseq\tstart\tend\tdirection\td\n'
                           'a\tb\tchr1\t1\t10\tab\tNA\n')
        fig = figure(summarize(table, Selection()))
        try:
            self.assertIsNone(fig.axes[0].get_legend())
            self.assertTrue(any(text.get_text() == 'No distance estimates' for text in fig.axes[0].texts))
        finally:
            plt.close(fig)

    def test_absent_matrix_pairs_remain_missing_and_repeats_are_averaged(self):
        table = self.table('genome_a\tgenome_b\td\na\tb\t0.2\na\tb\t0.4\na\tc\tNA\n')
        report = summarize(table, Selection())
        matrix = report.data['matrix']
        self.assertEqual(matrix.shape, (3, 3))
        self.assertAlmostEqual(matrix[0, 1], 0.3)
        self.assertTrue(np.isnan(matrix[0, 2]))
        self.assertEqual(report.data['y_names'], ['a', 'b', 'c'])

    def test_complete_all_pairs_input_keeps_a_symmetric_ordered_matrix(self):
        table = self.table('genome_a\tgenome_b\td\na\tb\t0.2\na\tc\t0.3\nb\tc\t0.4\n')
        report = summarize(table, Selection())
        matrix = report.data['matrix']
        self.assertTrue(report.data['symmetric'])
        self.assertEqual(matrix.shape, (3, 3))
        np.testing.assert_allclose(matrix, matrix.T)
        self.assertTrue(np.isnan(np.diag(matrix)).all())
        self.assertEqual(report.data['y_names'], ['a', 'b', 'c'])

    def test_large_genome_sets_show_closest_pairs_instead_of_a_matrix(self):
        from gvizz.render import figure, plt
        table = self.table('genome_a\tgenome_b\td\n' + ''.join(f'q\tr{i}\t0.{i % 9 + 1}\n' for i in range(100)))
        report = summarize(table, Selection())
        self.assertIsNone(report.data['matrix'])
        self.assertEqual(len(report.data['closest']), 20)
        fig = figure(report)
        try:
            fig.canvas.draw()
            self.assertEqual(len(fig.axes[0].patches), 20)
        finally:
            plt.close(fig)

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


class ModeTests(unittest.TestCase):
    """Interval tracks, overlays, annotations, filters, traversal, and report output."""

    ROLL = ('seq\tstart\tend\treference\td\n'
            'chr1\t1\t10\tref\t0.2\n'
            'chr1\t11\t20\tref\t0.4\n')
    MAP = ('seq\tseq_len\tstart\tend\treference\td\tside\tpi\tq_interval\n'
           'chr1\t1000\t3\t8\tref\t0.1\tlow\t0.01\t0.02\n'
           'chr1\t1000\t12\t18\tref\tNA\thigh\t0.001\tNA\n'
           'chr1\t1000\t30\t40\tref\t0.3\tlow\t0.01\t0.4\n')

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def write(self, name, text):
        path = self.root / name
        path.write_text(text)
        return path

    def test_overlay_tracks_join_the_matching_reference(self):
        from gvizz.render import figure, plt
        roll = Table(self.write('roll.tsv', self.ROLL))
        suffixed = self.MAP.replace('\tref\t', '\tref_contigs.fasta.gz\t')
        report = summarize(roll, Selection(), overlay=Table(self.write('m.tsv', suffixed)))
        self.assertEqual(report.data['overlay_kind'], 'map')
        self.assertEqual({key[0] for key in report.data['sequences'][0]['overlay']}, {'ref'})
        fig = figure(report)
        try:
            tracks = [ax for ax in fig.axes if ax.get_ylabel() in LABELS]
            self.assertEqual(len(tracks), 2)  # roll track + map overlay track
            self.assertEqual({ax.get_ylabel() for ax in tracks}, {'$\\hat{D}$', 'lower tail'})
        finally:
            plt.close(fig)
        with self.assertRaisesRegex(ValueError, 'different kind'):
            summarize(roll, Selection(), overlay=Table(self.write('r2.tsv', self.ROLL)))

    def test_map_filters_are_gone_and_rows_pass_through(self):
        table = Table(self.write('m.tsv', self.MAP))
        self.assertEqual(summarize(table, Selection()).rows, 3)  # no row filters remain

    def test_annotations_tsv_gff3_and_sequence_filter(self):
        from gvizz.data import read_annotations
        from gvizz.render import figure, plt
        tsv = self.write('genes.tsv', 'seq\tstart\tend\tstrand\tgene\nchr1\t3\t8\t-\tub\nchr2\t1\t5\t+\telsewhere\n')
        features, truncated = read_annotations(tsv)
        self.assertFalse(truncated)
        self.assertEqual([f.name for f in features], ['ub', 'elsewhere'])
        gff = self.write('genes.gff3', 'chr1\tsrc\tgene\t3\t8\t.\t-\t.\tID=g1;Name=ub\n')
        self.assertEqual(read_annotations(gff)[0][0].name, 'ub')
        features, truncated = read_annotations(tsv, limit=1)
        self.assertTrue(truncated)
        report = summarize(Table(self.write('r.tsv', self.ROLL)), Selection('chr1'),
                           annotations=read_annotations(tsv))
        self.assertEqual([f.name for f in report.data['sequences'][0]['genes']], ['ub'])
        fig = figure(report)
        try:
            self.assertEqual(fig.axes[-2].get_ylabel(), 'genes')  # before the colorbar axis
        finally:
            plt.close(fig)

    def test_dist_matrix_orders_genomes_by_mean_distance(self):
        from gvizz.render import figure, plt
        table = Table(self.write('d.tsv', 'genome_a\tgenome_b\td\n'
                                          'b\tc\t0.4\n'
                                          'a\tb\t0.1\n'
                                          'a\tc\t0.2\n'))
        report = summarize(table, Selection())
        self.assertTrue(report.data['symmetric'])
        self.assertEqual(report.data['y_names'], ['a', 'b', 'c'])  # ordered by mean distance
        self.assertAlmostEqual(report.data['matrix'][0, 1], 0.1)
        fig = figure(report)
        try:
            fig.canvas.draw()
        finally:
            plt.close(fig)

    def test_heading_and_stats(self):
        report = summarize(Table(self.write('r.tsv', self.ROLL)), Selection(), heading='Custom')
        self.assertEqual(report.heading, 'Custom')
        self.assertTrue(any(label == 'Selected rows' for label, _ in report.stats))

    def test_html_viewer_embeds_data_controls_and_no_stat_cards(self):
        from gvizz.render import write_report, viewer_payload
        table = Table(self.write('r.tsv', '# invocation: gdiff roll\n'
                                          'seq\tstart\tend\treference\td\n'
                                          'chr1\t1\t10\tref\t0.2\nchr2\t1\t10\tref\t0.3\n'))
        report = summarize(table, Selection())
        payload = viewer_payload(report)
        self.assertEqual([s['seq'] for s in payload['sequences']], ['chr1', 'chr2'])
        self.assertEqual(payload['sequences'][0]['tracks'][0]['segs'], [[1, 10, 0.2]])
        self.assertGreater(payload['vmax'], 0)
        out = self.root / 'report.html'
        write_report(report, table, out)
        html = out.read_text()
        for needle in ('gvizz-data', 'id="prev"', 'id="next"', 'id="filters"', 'id="plot"'):
            self.assertIn(needle, html)
        self.assertNotIn('class="card"', html)
        self.assertIn('invocation: gdiff roll', html)
        self.assertIn('lower distance is drawn on top', html)

    def test_cli_requires_seq_for_figure_output_and_writes_html_without_it(self):
        from gvizz.cli import main
        table = self.write('c.tsv', self.ROLL)
        with self.assertRaises(SystemExit):
            main([str(table), '-o', str(self.root / 'out.png')])
        self.assertEqual(main([str(table), '-o', str(self.root / 'out.html')]), 0)
        self.assertEqual(main([str(table), '--seq', 'chr1', '-o', str(self.root / 'out.png')]), 0)
        self.assertTrue((self.root / 'out.png').stat().st_size > 0)

    def test_cli_overlay_annotations_and_dist_seq_guard(self):
        from gvizz.cli import main
        roll = self.write('roll.tsv', self.ROLL)
        mapped = self.write('map.tsv', self.MAP)
        genes = self.write('genes.tsv', 'seq\tstart\tend\tstrand\tgene\nchr1\t3\t8\t-\tub\n')
        out = self.root / 'overlay.html'
        self.assertEqual(main([str(roll), str(mapped), '--annotations', str(genes),
                               '-o', str(out)]), 0)
        html = out.read_text()
        self.assertIn('gvizz-data', html)
        self.assertIn('lower tail', html)
        dist = self.write('dist.tsv', 'genome_a\tgenome_b\td\na\tb\t0.2\n')
        with self.assertRaises(SystemExit):
            main([str(dist), '--seq', 'chr1', '-o', str(self.root / 'dist.html')])


if __name__ == '__main__':
    unittest.main()
