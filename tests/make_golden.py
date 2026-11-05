#!/usr/bin/env python
"""Record what GlinkV2 outputs, for the GoogleTest suite to compare GlinkV3 against.

    python tests/make_golden.py    # with GlinkV2's environment active

This is the only file that touches GlinkV2. The tests read just what it writes, so they need neither
Python nor GlinkV2. Rerun it only when the fixtures or the cases change, and review what moves:
GlinkV2 is the reference, so a changed value here means the comparison changed.

--out, --sites and --samples record a larger graph elsewhere; GLINK_GOLDEN_DIR points the tests at it.

Layout: JSON for structure and small arrays (floats as float.hex, so exact), raw little-endian files
for large ones. Every array comes from GlinkV2's own command line or functions.
"""
import argparse
import json
import os
import random
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

from glink import _engine, conditioning, coordinates, core, dosage_ld, matching
from glink.args import node_filter

HERE = Path(__file__).resolve().parent
V2_BIN = Path(sys.executable).parent
DEFAULT_VG = 'vg'


def hexf(x):
    return float(x).hex()


def dump(path, value):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(value, separators=(',', ':')) + '\n')


def upper_triangle(matrix):
    """Both LD matrices are symmetric, so the files hold rows i, columns i.. only."""
    assert np.array_equal(matrix, matrix.T, equal_nan=True)
    return matrix[np.triu_indices(len(matrix))]


def check(condition, message):
    """A fixture that does not exercise its feature would make the tests prove nothing."""
    if not condition:
        raise SystemExit(f'make_golden: fixture too weak: {message}')


def run(*command, cwd):
    return subprocess.run([str(c) for c in command], cwd=cwd, capture_output=True, text=True, errors='replace')


# ---- the graph ------------------------------------------------------------------------------


def graph_golden(out, vg):
    ids, names, bits = core.read_gfa(out / 'graph.gfa')
    bits = np.asarray(bits, dtype=np.uint8)
    strand_ids, strand_names, _, (forward, reverse), summary = core.read_gfa(out / 'graph.gfa', oriented=True)
    assert np.array_equal(ids, strand_ids) and names == strand_names
    check(summary['flipped'] > 0, 'no walk was written on the opposite strand')
    check(summary['unanchored'] > 0, 'every walk touches the reference')
    check(any(n.endswith('#3') for n in names), 'no third haplotype')
    check(not any(n.startswith(('CHM13', 'GRCh38')) for n in names), 'a reference walk became a column')
    (out / 'graph.bits.u8').write_bytes(bits.tobytes())
    (out / 'graph.forward.u8').write_bytes(np.asarray(forward, dtype=np.uint8).tobytes())
    (out / 'graph.reverse.u8').write_bytes(np.asarray(reverse, dtype=np.uint8).tobytes())
    tips = {}
    for backbone in ('CHM13', 'HG000', 'HG001#2#ctg1', 'NOBODY'):
        tips[backbone] = [int(t) for t in conditioning.read_topology(out / 'graph.gfa', backbone)[2]]
    check(tips['CHM13'] and tips['HG001#2#ctg1'] and not tips['NOBODY'], 'backbone tips')
    dump(out / 'graph.json', dict(ids=ids.tolist(), names=list(names), rows=bits.shape[0], width=bits.shape[1],
                                  walks=summary['walks'], flipped=summary['flipped'],
                                  unanchored=summary['unanchored'], tips=tips))


# ---- coordinates ----------------------------------------------------------------------------

# Each small graph takes a route through the coordinate reader the generated graph does not.
SMALL_GRAPHS = {
    'walk_offset': 'S\t1\tAAAA\nS\t2\tCC\nS\t3\tGGG\nS\t4\tTT\nL\t1\t+\t2\t+\t0M\nL\t1\t+\t3\t+\t0M\nL\t2\t+\t4\t+\t0M\nL\t3\t+\t4\t+\t0M\n'
                   'W\tCHM13\t0\tchr1\t100\t108\t>1>2>4\nW\tA\t1\tc\t0\t9\t>1>3>4\n',
    'two_chromosomes': 'S\t1\tAAAA\nS\t2\tCC\nS\t3\tGGG\nS\t4\tTT\nS\t5\tA\nL\t1\t+\t2\t+\t0M\nL\t3\t+\t4\t+\t0M\nL\t2\t+\t5\t+\t0M\nL\t5\t+\t3\t+\t0M\n'
                       'W\tCHM13\t0\tchr1\t0\t6\t>1>2\nW\tCHM13\t0\tchr2\t0\t5\t>3>4\nW\tA\t1\tc\t0\t12\t>1>2>5>3>4\n',
    'length_tags': 'S\t1\t*\tLN:i:40\nS\t2\t*\tLN:i:5\nS\t3\tACGT\tLN:i:4\nS\t4\t*\tLN:i:7\nL\t1\t+\t2\t+\t0M\nL\t2\t+\t4\t+\t0M\nL\t1\t+\t3\t+\t0M\nL\t3\t+\t4\t+\t0M\n'
                   'W\tCHM13\t0\tchr1\t0\t52\t>1>2>4\nW\tA\t1\tc\t0\t51\t>1>3>4\n',
    'path_reference': 'S\t1\tAAAA\nS\t2\tCC\nS\t3\tGGG\nS\t4\tTT\nL\t1\t+\t2\t+\t0M\nL\t2\t+\t4\t+\t0M\nL\t1\t+\t3\t+\t0M\nL\t3\t+\t4\t+\t0M\n'
                      'P\tCHM13#0#chr1\t1+,2+,4+\t*\nP\tA#1#c\t1+,3+,4+\t0M,0M\n',
    'path_subrange': 'S\t1\tAAAA\nS\t2\tCC\nS\t3\tGGG\nS\t4\tTT\nL\t1\t+\t2\t+\t0M\nL\t2\t+\t4\t+\t0M\nL\t1\t+\t3\t+\t0M\nL\t3\t+\t4\t+\t0M\n'
                     'P\tCHM13#0#chr1[50-58]\t1+,2+,4+\t*\nP\tA#1#c\t1+,3+,4+\t*\n',
    'reversed_reference': 'S\t1\tAAAA\nS\t2\tCC\nS\t3\tGGG\nS\t4\tTT\nL\t4\t-\t2\t-\t0M\nL\t2\t-\t1\t-\t0M\nL\t1\t+\t3\t+\t0M\nL\t3\t+\t4\t+\t0M\n'
                          'W\tCHM13\t0\tchr1\t0\t8\t<4<2<1\nW\tA\t1\tc\t0\t9\t>1>3>4\n',
    'rgfa': 'S\t1\tAAAA\tSN:Z:chr1\tSO:i:0\tSR:i:0\nS\t2\tCC\tSN:Z:chr1\tSO:i:4\tSR:i:0\nS\t3\tGGG\tSN:Z:alt\tSO:i:0\tSR:i:1\nS\t4\tTT\tSN:Z:chr1\tSO:i:6\tSR:i:0\n'
            'L\t1\t+\t2\t+\t0M\nL\t2\t+\t4\t+\t0M\nL\t1\t+\t3\t+\t0M\nL\t3\t+\t4\t+\t0M\nP\tA#1#c\t1+,3+,4+\t*\n',
    'unplaced': 'S\t1\tAAAA\nS\t2\tCC\nS\t3\tGGG\nS\t9\tTT\nL\t1\t+\t2\t+\t0M\nL\t1\t+\t3\t+\t0M\nW\tCHM13\t0\tchr1\t0\t6\t>1>2\nW\tA\t1\tc\t0\t7\t>1>3\nW\tB\t1\tc\t0\t2\t>9\n',
    'star_overlap_link': 'S\t1\tAAAA\nS\t2\tCC\nL\t1\t+\t2\t+\t*\nW\tCHM13\t0\tchr1\t0\t6\t>1>2\nW\tA\t1\tc\t0\t4\t>1\n',
    'walk_span_differs': 'S\t1\tAAAA\nS\t2\tCC\nL\t1\t+\t2\t+\t0M\nW\tCHM13\t0\tchr1\t0\t99\t>1>2\nW\tA\t1\tc\t0\t4\t>1\n',
    'bad_walk_start': 'S\t1\tAAAA\nW\tCHM13\t0\tchr1\tx\t4\t>1\nW\tA\t1\tc\t0\t4\t>1\n',
    'bad_subrange': 'S\t1\tAAAA\nP\tCHM13#0#chr1[50]\t1+\t*\nW\tA\t1\tc\t0\t4\t>1\n',
    'undefined_link': 'S\t1\tAAAA\nL\t1\t+\t7\t+\t0M\nW\tCHM13\t0\tchr1\t0\t4\t>1\nW\tA\t1\tc\t0\t4\t>1\n',
    'undefined_step': 'S\t1\tAAAA\nW\tCHM13\t0\tchr1\t0\t4\t>7\nW\tA\t1\tc\t0\t4\t>1\n',
    'overlap_link': 'S\t1\tAAAA\nS\t2\tCC\nL\t1\t+\t2\t+\t5M\nW\tCHM13\t0\tchr1\t0\t6\t>1>2\nW\tA\t1\tc\t0\t4\t>1\n',
    'overlap_path': 'S\t1\tAAAA\nS\t2\tCC\nL\t1\t+\t2\t+\t0M\nP\tCHM13#0#chr1\t1+,2+\t5M\nW\tA\t1\tc\t0\t4\t>1\n',
    'length_mismatch': 'S\t1\tAAAA\tLN:i:9\nW\tCHM13\t0\tchr1\t0\t4\t>1\nW\tA\t1\tc\t0\t4\t>1\n',
    'no_length': 'S\t1\t*\nW\tCHM13\t0\tchr1\t0\t4\t>1\nW\tA\t1\tc\t0\t4\t>1\n',
    'duplicate_node': 'S\t1\tAAAA\nS\t1\tAAAA\nW\tCHM13\t0\tchr1\t0\t4\t>1\nW\tA\t1\tc\t0\t4\t>1\n',
    'no_reference': 'S\t1\tAAAA\nS\t2\tCC\nL\t1\t+\t2\t+\t0M\nW\tA\t1\tc\t0\t6\t>1>2\n',
}

# (graph, reference, source)
COORDINATE_RUNS = [
    ('graph.gfa', 'CHM13', 'auto'), ('graph.gfa', 'CHM13', 'walk'), ('graph.gfa', 'CHM13#0', 'auto'),
    ('graph.gfa', 'GRCh38', 'auto'), ('graph.gfa', 'HG001#2#ctg1', 'auto'), ('graph.gfa', 'CHM13', 'rgfa'),
    ('graph.gfa', 'NOBODY', 'auto'),
    *[(f'{name}.gfa', 'CHM13', 'auto') for name in SMALL_GRAPHS],
    ('walk_offset.gfa', 'CHM13', 'walk'), ('rgfa.gfa', 'CHM13', 'walk'), ('rgfa.gfa', 'CHM13', 'rgfa'),
    ('path_reference.gfa', 'CHM13#0', 'walk'), ('walk_offset.gfa', 'A', 'auto'),
]


def coordinate_golden(out, work):
    runs = []
    for name, text in SMALL_GRAPHS.items():
        (out / 'small').mkdir(exist_ok=True)
        (out / 'small' / f'{name}.gfa').write_text(text)
    for k, (graph, reference, source) in enumerate(COORDINATE_RUNS):
        gfa = out / graph if graph == 'graph.gfa' else out / 'small' / graph
        entry = dict(graph=str(gfa.relative_to(out)), reference=reference, source=source)
        try:
            built, _ = coordinates.prepare_coordinates(gfa, work / f'coords{k}.npz', reference=reference,
                                                       reference_source=source, rebuild=True)
            entry['arrays'] = {key: np.asarray(value).tolist() for key, value in built.arrays.items()}
            entry['summary'] = built.metadata['summary']
            for key in ('links', 'reference_nodes', 'reference_source', 'nodes'):
                entry[key] = built.metadata[key]
        except ValueError as error:
            entry['error'] = str(error)
        runs.append(entry)
    placed = [r for r in runs if 'arrays' in r]
    failed = [r for r in runs if 'error' in r]
    check(len(placed) > 8 and len(failed) >= 5, f'coordinates: {len(placed)} built, {len(failed)} failed')
    check(any(r['summary']['ambiguous'] for r in placed), 'no ambiguous node')
    check(any(r['summary']['reversed_anchors'] for r in placed), 'no reversed anchor')
    check(any(r['reference_source'] == 'rgfa_SO' for r in placed), 'no rGFA run')
    dump(out / 'coordinates.json', runs)


# ---- the conditioning tree ------------------------------------------------------------------


def tree_arrays(gfa, records, backbone='CHM13'):
    ids, adjacency, tips, _ = conditioning.read_topology(gfa, backbone)
    fields = conditioning.build_snarl_tree(records, ids)
    owner = conditioning.assign_nodes(fields, ids, adjacency)
    return dict(start_node=fields['start_node'].tolist(), end_node=fields['end_node'].tolist(),
                parent=fields['parent'].tolist(), depth=fields['depth'].tolist(), node_id=ids.tolist(),
                node_snarl=owner.tolist())


def tree_golden(out, vg, work):
    """vg runs once and its output is kept, so the comparison does not hang on the vg version."""
    gfa = out / 'graph.gfa'
    tips = conditioning.read_topology(gfa, 'CHM13')[2]
    records, _ = conditioning.run_vg_snarls(gfa, 'CHM13', tips, vg, 1)
    (out / 'graph.snarls.json').write_text(''.join(json.dumps(r) + '\n' for r in records))
    with open(out / 'snarls.bin', 'wb') as binary:
        command = [vg, 'snarls', '-n', '-P', 'CHM13', '-t', '1']
        for tip in tips:
            command += ['-w', str(tip)]
        subprocess.run(command + [str(gfa)], stdout=binary, check=True)
    tree = tree_arrays(gfa, records)
    check(max(tree['depth'], default=-1) >= 1, 'no nested snarl')
    check(any(s >= 0 for s in tree['node_snarl']), 'no node has a parent')
    starts = np.asarray(tree['start_node'])
    owner = np.asarray(tree['node_snarl'])
    parents = np.full(len(owner), -1, dtype=np.int64)
    parents[owner >= 0] = starts[owner[owner >= 0]]
    tree['parent_node_id'] = parents.tolist()
    dump(out / 'tree.json', tree)

    plain = out / 'plain.gfa'
    records, _ = conditioning.run_vg_snarls(plain, 'CHM13', conditioning.read_topology(plain, 'CHM13')[2], vg, 1)
    (out / 'plain.snarls.json').write_text(''.join(json.dumps(r) + '\n' for r in records))

    # On a chain the backbone changes nothing, so this graph forks: each backbone ends on its own tip.
    forked = out / 'small' / 'forked.gfa'
    forked.write_text('S\t1\tAA\nS\t2\tC\nS\t3\tG\nS\t4\tAA\nS\t5\tC\nS\t6\tG\nS\t7\tAA\nS\t8\tTT\nS\t9\tC\nS\t10\tG\nS\t11\tTT\n'
                      'L\t1\t+\t2\t+\t0M\nL\t1\t+\t3\t+\t0M\nL\t2\t+\t4\t+\t0M\nL\t3\t+\t4\t+\t0M\nL\t4\t+\t5\t+\t0M\nL\t4\t+\t6\t+\t0M\n'
                      'L\t5\t+\t7\t+\t0M\nL\t6\t+\t7\t+\t0M\nL\t4\t+\t8\t+\t0M\nL\t8\t+\t9\t+\t0M\nL\t8\t+\t10\t+\t0M\nL\t9\t+\t11\t+\t0M\nL\t10\t+\t11\t+\t0M\n'
                      'W\tCHM13\t0\tchr1\t0\t8\t>1>2>4>5>7\nW\tOTHER\t1\tc\t0\t9\t>1>3>4>8>9>11\nW\tA\t1\tc\t0\t8\t>1>2>4>6>7\nW\tA\t2\tc\t0\t9\t>1>3>4>8>10>11\n')
    backbones = {}
    for backbone in ('CHM13', 'OTHER', 'A'):
        fields, ids, owner, _, _ = conditioning.build(forked, backbone, vg, 1)
        backbones[backbone] = dict(start_node=fields['start_node'].tolist(), end_node=fields['end_node'].tolist(),
                                   parent=fields['parent'].tolist(), node_id=ids.tolist(), node_snarl=owner.tolist())
    check(backbones['CHM13'] != backbones['OTHER'], f'the backbone changes nothing: {backbones}')
    dump(out / 'tree_backbones.json', backbones)

    # --vg-graph: the same steps through a GBZ and its --no-translation export.
    subprocess.run([vg, 'gbwt', '-G', str(gfa), '--gbz-format', '-g', str(out / 'graph.gbz')], check=True, capture_output=True)
    with open(out / 'graph.gbz.gfa', 'w') as export:
        subprocess.run([vg, 'convert', '-f', '--no-translation', str(out / 'graph.gbz')], stdout=export, check=True)
    gbz_gfa = out / 'graph.gbz.gfa'
    tips = conditioning.read_topology(gbz_gfa, 'CHM13')[2]
    check(tips, 'the GBZ export lost the CHM13 walk')
    fields, ids, owner, _, _ = conditioning.build(gbz_gfa, 'CHM13', vg, 1, out / 'graph.gbz')
    dump(out / 'tree_gbz.json', dict(start_node=fields['start_node'].tolist(), end_node=fields['end_node'].tolist(),
                                     parent=fields['parent'].tolist(), node_id=ids.tolist(), node_snarl=owner.tolist()))


# ---- whole runs -----------------------------------------------------------------------------


def case(matrix_type='phased', filters=(), collapse=True, assume_complete=False, tree=False, coords=True, condition=False,
         collapse_backbone=False, broken=6, split=False, drop=False, samples=None, min_pairs=2, min_allele_count=0,
         blocks=False, max_kb=200, neutral=False, r2=True, alias=None, graph='graph'):
    return dict(graph=graph, matrix_type=matrix_type, filters=list(filters), collapse=collapse, assume_complete=assume_complete, tree=tree,
                coords=coords, condition=condition, collapse_backbone=collapse_backbone, broken=broken, split=split, drop=drop,
                samples=samples, min_pairs=min_pairs, min_allele_count=min_allele_count, blocks=blocks, max_kb=max_kb,
                neutral=neutral, r2=r2, alias=alias)


# Conditioned on parents means with the tree.
def conditioned(**k):
    return case(tree=True, condition=True, **k)


PLACED = '@nodes.txt'
CASES = {
    # matrix types and absence models
    'phased_complete': case(assume_complete=True),
    'phased_observed': case(),
    'phased_alias': case(assume_complete=True, alias='haplotype'),
    'dosage_complete': case('dosage', assume_complete=True),
    'dosage_observed': case('dosage'),
    'presence_complete': case('presence', assume_complete=True),
    'presence_observed': case('presence'),
    # filters
    'band': case('dosage', ['0.1:0.9'], assume_complete=True),
    'band_stack': case('presence', ['0.5:0.95', 'minor:0.05'], False, True),
    'majority': case('dosage', ['majority'], assume_complete=True),
    'polymorphic': case('phased', ['polymorphic'], assume_complete=True),
    'all_explicit': case('phased', ['all'], assume_complete=True),
    'minor': case('phased', ['minor:0.1'], False, True),
    'minor_dosage': case('dosage', ['minor:0.1'], False, True),
    'list_then_minor': case('phased', [PLACED, 'minor:0.05'], False, True),
    'two_lists': case('phased', [PLACED, '@few.txt'], False, True),
    'common_complete': case('phased', ['common'], assume_complete=True),
    'common_presence': case('presence', ['common'], assume_complete=True),
    # samples
    'samples_phased': case('phased', ['minor:0.05'], False, True, samples='samples.txt'),
    'samples_dosage': case('dosage', ['minor:0.05'], False, True, samples='samples.txt'),
    'samples_observed': case('dosage', samples='samples.txt'),
    # parents
    'parents_phased': conditioned(filters=['minor:0.05'], collapse=False),
    'parents_dosage': conditioned(matrix_type='dosage', filters=['minor:0.05'], collapse=False),
    'parents_presence': conditioned(matrix_type='presence', filters=['0.1:0.9']),
    'parents_samples': conditioned(filters=[PLACED, 'minor:0.05'], collapse=False, samples='samples.txt'),
    'tree_not_conditioned': case('dosage', assume_complete=True, tree=True),
    # collapsing parents
    'backbone_phased': conditioned(collapse_backbone=True, filters=['minor:0.05'], collapse=False),
    'backbone_dosage': conditioned(matrix_type='dosage', collapse_backbone=True, filters=['common'], collapse=False),
    'broken_1': conditioned(collapse_backbone=True, broken=1, collapse=False),
    'broken_2': conditioned(collapse_backbone=True, broken=2, collapse=False),
    'broken_100': conditioned(collapse_backbone=True, broken=100, collapse=False),
    # on a graph with no truncated haplotypes, so chain nodes are parents everyone visits
    'plain_conditioned': conditioned(graph='plain', collapse=False),
    'plain_backbone': conditioned(graph='plain', collapse_backbone=True, collapse=False),
    'plain_backbone_dosage': conditioned(graph='plain', matrix_type='dosage', collapse_backbone=True, filters=['common'], collapse=False),
    'plain_gabriel': conditioned(graph='plain', collapse_backbone=True, filters=['common'], collapse=False, blocks=True, neutral=True),
    # complements
    'drop_phased': case('phased', ['minor:0.02'], assume_complete=True, drop=True),
    'drop_dosage': case('dosage', ['minor:0.02'], assume_complete=True, drop=True),
    'drop_plain': case('phased', ['minor:0.02'], assume_complete=True, drop=True, graph='plain'),
    'drop_plain_dosage': case('dosage', ['minor:0.02'], False, True, drop=True, graph='plain'),
    'split_drop_plain': case('phased', ['minor:0.03'], assume_complete=True, split=True, drop=True, graph='plain'),
    'neutral_plain': case('phased', ['common'], False, True, blocks=True, neutral=True, graph='plain'),
    'neutral_plain_off': case('phased', ['common'], False, True, blocks=True, graph='plain'),
    'drop_parents': conditioned(filters=['minor:0.02'], drop=True),
    'drop_parents_dosage': conditioned(matrix_type='dosage', filters=['minor:0.02'], drop=True, collapse=False),
    'drop_samples': conditioned(filters=['minor:0.02'], drop=True, samples='samples.txt'),
    # inversions
    'split': case('phased', assume_complete=True, split=True),
    'split_drop': case('phased', ['minor:0.03'], assume_complete=True, split=True, drop=True),
    'split_dosage': case('dosage', ['minor:0.03'], False, True, split=True),
    'split_parents': conditioned(split=True, filters=['minor:0.03'], collapse=False),
    'split_observed': case('phased', split=True),
    # LD support
    'strict_pairs': conditioned(filters=['minor:0.05'], collapse=False, min_pairs=60, min_allele_count=3),
    'pairs_at_columns': case('phased', ['minor:0.05'], False, True, min_pairs=121),
    'pairs_over_columns': case('phased', ['minor:0.05'], False, True, min_pairs=122),
    'allele_count': case('dosage', ['minor:0.02'], False, True, min_allele_count=8),
    'no_coordinates': case('dosage', ['minor:0.05'], False, True, coords=False),
    # blocks
    'gabriel_phased': case('phased', [PLACED, 'minor:0.05'], False, True, blocks=True),
    'gabriel_dosage': case('dosage', [PLACED, 'minor:0.05'], False, True, blocks=True),
    'gabriel_collapsed': case('phased', [PLACED, 'minor:0.05'], True, True, blocks=True),
    'gabriel_span': case('phased', [PLACED, 'minor:0.05'], False, True, blocks=True, max_kb=2),
    'gabriel_span_dosage': case('dosage', [PLACED, 'minor:0.05'], False, True, blocks=True, max_kb=0.5),
    'gabriel_no_r2': case('dosage', [PLACED, 'minor:0.05'], False, True, blocks=True, r2=False),
    'gabriel_parents': conditioned(filters=[PLACED, 'minor:0.05'], collapse=False, blocks=True),
    'gabriel_parents_dosage': conditioned(matrix_type='dosage', filters=[PLACED, 'minor:0.05'], collapse=False, blocks=True),
    'gabriel_samples': conditioned(filters=[PLACED, 'minor:0.05'], collapse=False, samples='samples.txt', blocks=True, max_kb=5),
    'neutral_phased': conditioned(filters=[PLACED, 'minor:0.05'], drop=False, collapse=False, blocks=True, neutral=True),
    'neutral_dosage': case('dosage', [PLACED, 'minor:0.05'], False, True, blocks=True, neutral=True, max_kb=5),
    'neutral_complements': case('phased', [PLACED, 'minor:0.02'], False, True, blocks=True, neutral=True),
    'paper_phased': conditioned(filters=[PLACED, 'common'], collapse=False, collapse_backbone=True, blocks=True, neutral=True),
    'paper_dosage': conditioned(matrix_type='dosage', filters=[PLACED, 'common'], collapse=False, collapse_backbone=True,
                                min_pairs=50, min_allele_count=5, blocks=True, neutral=True, max_kb=100),
}


def v2_command(spec, directory):
    """GlinkV2's spelling of a case, run from the golden directory."""
    gfa, snarls = f"{spec['graph']}.gfa", f"{spec['graph']}.snarls.json"
    command = [V2_BIN / 'glink', '--gfa', gfa, '--matrix-type', spec['alias'] or spec['matrix_type'],
               '--matrix-out', f'{directory}/matrix.npz']
    if spec['coords']:
        command += ['--node-coordinates', f'{directory}/coords.npz']
    else:
        command += ['--no-node-coordinates']
    if spec['tree']:
        command += ['--node-parents', f'{directory}/parents.npz', '--snarls', snarls]
    else:
        command += ['--no-node-parents']
    for spec_filter in spec['filters']:
        command += ['--filter', spec_filter]
    command += ['--collapse' if spec['collapse'] else '--no-collapse']
    for flag, key in (('--assume-complete', 'assume_complete'), ('--condition-on-parents', 'condition'),
                      ('--collapse-backbone-parents', 'collapse_backbone'), ('--split-inversions', 'split'),
                      ('--drop-complements', 'drop'), ('--neutral-complements', 'neutral')):
        if spec[key]:
            command += [flag]
    command += ['--broken-parent-min-haplotypes', spec['broken'], '--min-pairs', spec['min_pairs'],
                '--min-allele-count', spec['min_allele_count']]
    if spec['samples']:
        command += ['--samples', spec['samples']]
    if spec['r2']:
        command += ['--r2-out', f'{directory}/r2.npy']
    if spec['blocks']:
        command += ['--blocks', 'gabriel', '--blocks-max-kb', spec['max_kb'], '--blocks-out', f'{directory}/blocks.tsv']
    return command


def console_lines(text):
    """The lines of GlinkV2's report that say something about the results, without file paths."""
    return [line.split(' -> ')[0] for line in text.splitlines() if line.strip()]


def case_golden(name, spec, out, work, stats):
    directory = work / name
    directory.mkdir()
    done = run(*v2_command(spec, directory), cwd=out)
    if done.returncode:
        raise SystemExit(f'GlinkV2 failed on case {name}: {" ".join(map(str, v2_command(spec, directory)))}\n{done.stderr[-2000:]}')
    m = core.load_matrix(directory / 'matrix.npz')
    values = np.asarray(m.values, dtype=np.int8)
    check(values.size, f'{name}: no rows')
    destination = out / 'cases' / name
    destination.mkdir(parents=True)
    (destination / 'values.i8').write_bytes(values.tobytes())
    if spec['r2']:
        r2 = np.load(directory / 'r2.npy')
        assert r2.dtype == np.float32
        (destination / 'r2.f32').write_bytes(upper_triangle(r2).tobytes())
    d = dosage_ld.dprime(values) if m.matrix_type == 'dosage' else core.dprime(values)
    (destination / 'dprime.f64').write_bytes(upper_triangle(np.asarray(d, dtype=np.float64)).tobytes())
    if spec['blocks']:
        shutil.copy(directory / 'blocks.tsv', destination / 'blocks.tsv')
        labels = [line.split('\t')[4] for line in (directory / 'blocks.tsv').read_text().splitlines()[1:]]
        stats['blocks'] = stats.get('blocks', 0) + any(int(size) > 1 for size in labels)
    record = dict(rows=values.shape[0], cols=values.shape[1], matrix_type=m.matrix_type, columns=[str(c) for c in m.columns],
                  node_ids=np.asarray(m.node_ids).tolist(), input_node_ids=np.asarray(m.input_node_ids).tolist(),
                  labels=[str(label) for label in core.row_labels(m)], density=hexf(m.density()),
                  console=console_lines(done.stdout))
    if m.pattern_of is not None:
        record['pattern_of'] = np.asarray(m.pattern_of).tolist()
    if m.parent_node_ids is not None:
        record['parent_node_ids'] = np.asarray(m.parent_node_ids).tolist()
    if m.node_orientation is not None:
        record['node_orientation'] = np.asarray(m.node_orientation).tolist()
        record['input_node_orientation'] = np.asarray(m.input_node_orientation).tolist()
    if m.complement_pairs is not None:
        record['complement_pairs'] = np.asarray(m.complement_pairs).reshape(-1).tolist()
    if m.coordinates is not None:
        record['coordinates'] = {key: np.asarray(value).tolist() for key, value in m.coordinates.items()}
    callability = {k: v for k, v in m.callability.items() if k in ('absence_model', 'backbone_collapse')}
    if 'complements' in m.callability:
        callability['complements'] = dict(dropped=m.callability['complements']['dropped'])
    if 'orientation' in m.callability:
        callability['orientation'] = {k: v for k, v in m.callability['orientation'].items() if k != 'rule'}
    record['callability'] = callability
    dump(destination / 'matrix.json', record)
    return values


def require_coverage(out):
    """Each feature the cases claim to test has to change the result somewhere."""
    def load(name):
        return json.loads((out / 'cases' / name / 'matrix.json').read_text())

    collapse = {name: load(name)['callability']['backbone_collapse'] for name in ('broken_1', 'broken_2', 'broken_100')}
    check(collapse['broken_1']['broken_parents'] > collapse['broken_100']['broken_parents'], f'broken parents {collapse}')
    plain = load('plain_backbone')['callability']['backbone_collapse']
    check(plain['backbone_parents'] > 0 and plain['moved_nodes'] > 0, f'backbone parents {plain}')
    moved = load('plain_backbone')['parent_node_ids'] != load('plain_conditioned')['parent_node_ids']
    check(moved, 'collapsing backbone parents moved no parent')
    check(load('split')['callability']['orientation']['inverted_nodes'] > 0, 'no inverted node')
    for name in ('drop_plain', 'drop_plain_dosage', 'drop_parents', 'split_drop_plain'):
        check(load(name)['callability']['complements']['dropped'] > 0, f'{name}: no complement dropped')
    neutral = [(out / 'cases' / name / 'blocks.tsv').read_text() for name in ('neutral_plain', 'neutral_plain_off')]
    check(neutral[0] != neutral[1], '--neutral-complements changed no block')
    check(load('phased_complete')['rows'] < len(load('phased_complete')['input_node_ids']), 'collapse merged nothing')
    no_r2 = (out / 'cases' / 'gabriel_no_r2' / 'blocks.tsv').read_text().splitlines()[1].split('\t')
    check(no_r2[8] == '', 'blocks without r² carry statistics')
    for name, array in (('strict_pairs', 'r2.f32'), ('pairs_over_columns', 'r2.f32'), ('allele_count', 'r2.f32')):
        values = np.frombuffer((out / 'cases' / name / array).read_bytes(), dtype=np.float32)
        check(np.isnan(values).any(), f'{name}: every pair passed')
    all_nan = np.frombuffer((out / 'cases' / 'pairs_over_columns' / 'r2.f32').read_bytes(), dtype=np.float32)
    check(np.isnan(all_nan).all(), 'pairs_over_columns: expected every r² undefined')
    at_limit = np.frombuffer((out / 'cases' / 'pairs_at_columns' / 'r2.f32').read_bytes(), dtype=np.float32)
    check(np.isfinite(at_limit).any(), 'pairs_at_columns: every r² undefined')
    undefined_d = np.frombuffer((out / 'cases' / 'samples_dosage' / 'dprime.f64').read_bytes())
    check(np.isnan(undefined_d).all(), 'samples_dosage: 40 people should leave dosage D′ undefined')
    defined_d = np.frombuffer((out / 'cases' / 'gabriel_dosage' / 'dprime.f64').read_bytes())
    check(np.isfinite(defined_d).any(), 'gabriel_dosage: no defined D′')


# ---- engines --------------------------------------------------------------------------------


def count_tables(rng, ploidy):
    """Edge cases, then random tables of every density, around the 50-person support rule."""
    size = 2 if ploidy == 1 else 3

    def table(cells):
        t = np.zeros((3, 3), dtype=np.uint16)
        t[:size, :size] = np.asarray(cells).reshape(size, size)
        return t

    tables = [table(np.zeros(size * size, dtype=int))]
    for scale in (1, 5, 49, 50, 51, 500, 7000):
        for k in range(size * size):
            cells = np.zeros(size * size, dtype=int)
            cells[k] = scale
            tables.append(table(cells))
        tables.append(table(np.full(size * size, scale)))
        if ploidy == 1:
            tables += [table([scale, 0, 0, scale]), table([0, scale, scale, 0]), table([scale, 1, 1, scale]),
                       table([scale, scale, 0, 1])]
        else:
            tables += [table([scale, 0, 0, 0, scale, 0, 0, 0, scale]), table([0, 0, scale, 0, scale, 0, scale, 0, 0]),
                       table([scale, 1, 0, 1, scale, 1, 0, 1, scale])]
    for _ in range(1500):
        total = int(rng.choice([8, 30, 49, 50, 60, 120, 300, 2000, 60000]))
        cells = rng.multinomial(total, rng.dirichlet(np.full(size * size, float(rng.choice([.1, .5, 1, 5])))))
        tables.append(table(cells))
    return tables


def interval_row(table, result):
    return [table.reshape(9).tolist(), hexf(result[0]), hexf(result[1]), hexf(result[2]), int(result[3])]


def panel(rng, markers, columns, ploidy, missing, complements=False):
    """Markers in LD: mosaics of a few founder rows, with noise, and some exact complements."""
    founders = rng.random((6, columns))
    frequency = rng.uniform(.1, .9, markers)
    values = np.zeros((markers, columns), dtype=np.int8)
    f = 0
    for m in range(markers):
        if rng.random() < .1:
            f = int(rng.integers(6))
        calls = (founders[f] < frequency[m]).astype(np.int8)
        if ploidy == 2:
            calls = calls + (rng.permutation(founders[f]) < frequency[m]).astype(np.int8)
        noise = rng.random(columns) < .02
        calls[noise] = rng.integers(0, ploidy + 1, int(noise.sum()))
        values[m] = calls
    if complements:
        for m in range(1, markers, 4):
            values[m] = ploidy - values[m - 1]
    if missing:
        values[rng.random(values.shape) < missing] = -1
    return values


def packed(values):
    return np.stack([np.packbits(values == g, axis=1) for g in range(3)], axis=1)


def engine_golden(out):
    rng = np.random.default_rng(29)
    tables = dict(dosage=[], gabriel_diploid=[], gabriel_haploid=[])
    for t in count_tables(rng, 2):
        tables['dosage'].append(interval_row(t, _engine.dosage_confidence(t.reshape(9))))
        tables['gabriel_diploid'].append(interval_row(t, _engine.gabriel_confidence(t.reshape(9), 2)))
    for t in count_tables(rng, 1):
        tables['gabriel_haploid'].append(interval_row(t, _engine.gabriel_confidence(t.reshape(9), 1)))
    panels = []
    for name, (markers, columns, ploidy, missing, complements) in dict(
        diploid=(32, 70, 2, 0, False), diploid_missing=(32, 70, 2, .05, False), diploid_complements=(32, 70, 2, 0, True),
        haploid=(32, 120, 1, 0, False), haploid_missing=(32, 120, 1, .05, False), haploid_complements=(32, 120, 1, 0, True),
        wide=(16, 300, 2, .02, False),
    ).items():
        values = panel(rng, markers, columns, ploidy, missing, complements)
        planes = packed(values)
        positions = np.cumsum(rng.integers(1, 1500, markers)).astype(np.int64)
        entry = dict(name=name, markers=markers, width=planes.shape[2], ploidy=ploidy, packed=planes.tobytes().hex(),
                     positions=positions.tolist(), spans=[])
        if ploidy == 2:
            entry['classes'] = _engine.dosage_classes(planes, 1).reshape(-1).tolist()
            entry['dprime'] = [hexf(v) for v in _engine.dosage_dprime(planes, 1).reshape(-1)]
        for span in (0, 2000, 10000, 200000):
            for neutral in (False, True):
                bands = _engine.gabriel_bands(planes, positions, span, ploidy, 1, neutral)
                selected = _engine.gabriel_select(bands, positions, span)
                blocks = _engine.gabriel_blocks(planes, positions, span, ploidy, 1, neutral)
                assert np.array_equal(selected, blocks), name
                entry['spans'].append(dict(span=span, neutral=neutral, bands=bands.tobytes().hex(), blocks=blocks.tolist()))
        panels.append(entry)
    check(any(len(s['blocks']) > 1 for p in panels for s in p['spans']), 'no panel made two blocks')
    bands_seen = {b for p in panels for s in p['spans'] for b in bytes.fromhex(s['bands'])}
    check({0, 1, 2, 3, 4, 5, 6, 7} <= bands_seen, f'panels reach only bands {sorted(bands_seen)}')
    dump(out / 'engine.json', dict(tables=tables, panels=panels))


# ---- filters, block comparison, errors ------------------------------------------------------


def filter_golden(out):
    specs = ['all', 'majority', 'polymorphic', 'common', 'minor:0.05', 'minor:0', 'minor:0.49', 'minor:0.5', 'minor:-0.1',
             'minor:x', 'minor:', 'minor: 0.1', 'minor:1e-3', '0.1:0.9', '0:1', '0.9:0.1', '0.5:0.5', '-1:0.5', '0:2', 'a:b', ':',
             '0.1:', ' 0.1:0.9', '0.1:0.9 ', 'nonsense', '', '@nodes.txt', '@', 'minor:nan', 'nan:nan', 'minor:inf', '1e-1:1']
    rows = []
    for spec in specs:
        try:
            kind, value = node_filter(spec)
        except argparse.ArgumentTypeError as error:
            rows.append(dict(spec=spec, error=str(error)))
            continue
        row = dict(spec=spec, kind=kind)
        if kind == 'band':
            row.update(low=hexf(value[0]), high=hexf(value[1]))
        elif kind == 'minor':
            row.update(value=hexf(value))
        elif kind == 'common':
            row.update(value=hexf(value[0]), people=value[1])
        elif kind == 'list':
            row.update(path=str(value))
        rows.append(row)
    dump(out / 'filters.json', rows)


def block_table(rng, plink, count, wide):
    lefts = np.sort(rng.integers(0, 2000 if wide else 40, count))
    spans = rng.integers(0, 300 if wide else 6, count)
    if plink:
        rows = [f' 1 {left + 1} {left + span} {span / 1000} {rng.integers(2, 9)} a|b\n' for left, span in zip(lefts, spans)]
        return ' CHR BP1 BP2 KB NSNPS SNPS\n' + ''.join(rows)
    header = 'feature\tid\tposition\tblock\tblock_size\tblock_start\tblock_end\tblock_span\tblock_mean_r2\tblock_min_r2\n'
    rows, feature = [], 0
    for label, (left, span) in enumerate(zip(lefts, spans)):
        size = int(rng.integers(1, 4))
        for _ in range(size):
            rows.append(f'{feature}\t{feature}\t{left + 1}\t{label}\t{size}\t{left + 1}\t{left + span}\t{span}\t0.5\t0.5\n')
            feature += 1
    return header + ''.join(rows)


def compare_golden(out, work):
    rng = np.random.default_rng(11)
    rows = []
    for trial in range(60):
        wide = bool(trial % 2)
        a = block_table(rng, trial % 3 != 0, int(rng.integers(0, 30)), wide)
        b = block_table(rng, trial % 4 != 0, int(rng.integers(0, 30)), wide)
        (work / 'a').write_text(a)
        (work / 'b').write_text(b)
        region = ['chr1:11-900', 'chr1:1-2300', 'chr1:500-500', 'chrX:1000-1200'][trial % 4]
        done = run(V2_BIN / 'glink-compare-blocks', 'a', 'b', '--region', region, cwd=work)
        if done.returncode:
            raise SystemExit(f'glink-compare-blocks failed: {done.stderr[-500:]}')
        scores = json.loads(done.stdout)
        rows.append(dict(a=a, b=b, **{k: scores[k] for k in ('region', 'region_bp', 'a_blocks', 'b_blocks', 'a_bp', 'b_bp', 'shared_bp')}))
    check(any(r['shared_bp'] for r in rows) and any(r['a_bp'] != r['shared_bp'] for r in rows), 'block tables never overlap')
    bad_tables = {
        'empty': '',
        'unknown_header': 'what is this\n1 2 3\n',
        'det_non_numeric': ' CHR BP1 BP2 KB NSNPS SNPS\n 1 x y 0 2 a|b\n',
        'det_too_few_columns': ' CHR BP1 BP2 KB NSNPS SNPS\n 1\n',
        'glink_short_row': 'feature\tid\tposition\tblock\tblock_size\tblock_start\tblock_end\n0\t1\n',
        'glink_non_numeric': 'feature\tid\tposition\tblock\tblock_size\tblock_start\tblock_end\n0\t1\t2\t0\tx\ty\tz\n',
        'glink_empty_block': 'feature\tid\tposition\tblock\tblock_size\tblock_start\tblock_end\n0\t1\t\t\t\t\t\n',
    }
    bad = []
    good = ' CHR BP1 BP2 KB NSNPS SNPS\n 1 1 100 0.1 2 a|b\n'
    for name, text in bad_tables.items():
        (work / 'bad').write_text(text)
        (work / 'good').write_text(good)
        done = run(V2_BIN / 'glink-compare-blocks', 'bad', 'good', '--region', 'chr1:1-200', cwd=work)
        bad.append(dict(name=name, text=text, good=good, exit=done.returncode))
    regions = []
    for region in ('nonsense', 'chr1', 'chr1:', 'chr1:abc-def', 'chr1:5', 'chr1:100-50', 'chr1:0-0'):
        (work / 'good').write_text(good)
        done = run(V2_BIN / 'glink-compare-blocks', 'good', 'good', '--region', region, cwd=work)
        regions.append(dict(region=region, exit=done.returncode))
    dump(out / 'compare.json', dict(tables=rows, bad_tables=bad, bad_regions=regions))


BAD_GFAS = {
    'binary_garbage': bytes(range(256)) * 40,
    'walk_too_few_fields': b'H\tVN:Z:1.1\nS\t1\tA\nW\tS1\n',
    'walk_unknown_node': b'S\t1\tA\nW\tS1\t1\tc\t0\t1\t>99\n',
    'walk_not_a_path': b'S\t1\tA\nW\tS1\t1\tc\t0\t1\tnonsense\n',
    'segment_without_id': b'S\nW\tS1\t1\tc\t0\t1\t>1\n',
    'segment_id_not_a_number': b'S\tabc\tA\nW\tS1\t1\tc\t0\t1\t>abc\n',
    'segment_id_overflow': b'S\t99999999999999999999999\tA\nW\tS1\t1\tc\t0\t1\t>99999999999999999999999\n',
    'negative_segment_id': b'S\t-5\tA\nW\tS1\t1\tc\t0\t1\t>-5\n',
    'path_without_segments': b'S\t1\tA\nP\tname\n',
    'path_empty_step': b'S\t1\tA\nP\tname\t1+,,\t*\n',
    'no_haplotypes': b'H\tVN:Z:1.1\nS\t1\tA\nS\t2\tC\n',
    'only_reference': b'S\t1\tA\nW\tCHM13\t0\tchr1\t0\t1\t>1\n',
    'no_segments': b'H\tVN:Z:1.1\n',
    'empty': b'',
    'crlf': b'S\t1\tA\r\nS\t2\tC\r\nL\t1\t+\t2\t+\t0M\r\nW\tS1\t1\tc\t0\t2\t>1>2\r\nW\tS1\t2\tc\t0\t1\t>1\r\n',
    'no_final_newline': b'S\t1\tA\nS\t2\tC\nW\tS1\t1\tc\t0\t2\t>1>2\nW\tS1\t2\tc\t0\t1\t>1',
    'duplicate_segment': b'S\t1\tA\nS\t1\tA\nW\tS1\t1\tc\t0\t1\t>1\nW\tS1\t2\tc\t0\t1\t>1\n',
}


def error_golden(out, work):
    """GlinkV2's verdict on bad GFAs and bad command lines. V3 spells a few options differently."""
    bad = out / 'bad'
    bad.mkdir()
    gfas = []
    skeleton = ['--no-node-parents', '--no-node-coordinates', '--assume-complete', '--matrix-type', 'phased',
                '--matrix-out', str(work / 'bad.npz')]
    for name, content in BAD_GFAS.items():
        (bad / f'{name}.gfa').write_bytes(content)
        done = run(V2_BIN / 'glink', '--gfa', f'bad/{name}.gfa', *skeleton, cwd=out)
        message = done.stderr.strip().splitlines()[-1] if done.returncode and done.stderr.strip() else ''
        gfas.append(dict(name=name, exit=done.returncode, message=message))
    check(any(g['exit'] == 0 for g in gfas) and any(g['exit'] for g in gfas), 'bad GFAs: all one verdict')

    # In V3's spelling; `v2` is GlinkV2's where it differs. Each runs from the golden directory.
    no_graph_steps = ['--no-conditioning-tree', '--no-node-coordinates', '--assume-complete']
    no_graph_steps_v2 = ['--no-node-parents', '--no-node-coordinates', '--assume-complete']
    gfa = ['--gfa', 'graph.gfa']
    matrix = ['--matrix-out', 'OUT/m']
    r2 = ['--r2-out', 'OUT/r2.ld']
    r2_v2 = ['--r2-out', 'OUT/r2.npy']
    lines = dict(
        help=(['--help'], None),
        equals_form=(gfa + no_graph_steps + matrix + ['--threads=2', '--filter=minor:0.1'], gfa + no_graph_steps_v2 + matrix + ['--threads=2', '--filter=minor:0.1']),
        flag_with_value=(gfa + no_graph_steps + matrix + ['--assume-complete=1'], gfa + no_graph_steps_v2 + matrix + ['--assume-complete=1']),
        abbreviated_option=(gfa + no_graph_steps + matrix + ['--matrix-typ', 'dosage'], gfa + no_graph_steps_v2 + matrix + ['--matrix-typ', 'dosage']),
        no_args=([], None),
        unknown_option=(['--no-such-option'], None),
        gfa_without_value=(['--gfa'], None),
        missing_gfa=(['--gfa', 'absent.gfa'] + no_graph_steps + matrix, ['--gfa', 'absent.gfa'] + no_graph_steps_v2 + matrix),
        no_output=(gfa + no_graph_steps, gfa + no_graph_steps_v2),
        threads_zero=(gfa + no_graph_steps + matrix + ['--threads', '0'], gfa + no_graph_steps_v2 + matrix + ['--threads', '0']),
        threads_word=(gfa + no_graph_steps + matrix + ['--threads', 'many'], gfa + no_graph_steps_v2 + matrix + ['--threads', 'many']),
        bad_filter=(gfa + no_graph_steps + matrix + ['--filter', 'nonsense'], gfa + no_graph_steps_v2 + matrix + ['--filter', 'nonsense']),
        bad_minor=(gfa + no_graph_steps + matrix + ['--filter', 'minor:0.9'], gfa + no_graph_steps_v2 + matrix + ['--filter', 'minor:0.9']),
        bad_type=(gfa + no_graph_steps + matrix + ['--matrix-type', 'nonsense'], gfa + no_graph_steps_v2 + matrix + ['--matrix-type', 'nonsense']),
        absent_samples=(gfa + no_graph_steps + matrix + ['--samples', 'absent.txt'], gfa + no_graph_steps_v2 + matrix + ['--samples', 'absent.txt']),
        absent_list=(gfa + no_graph_steps + matrix + ['--filter', '@absent.txt'], gfa + no_graph_steps_v2 + matrix + ['--filter', '@absent.txt']),
        list_keeps_none=(gfa + no_graph_steps + matrix + ['--filter', '@bad/no_nodes.txt'], gfa + no_graph_steps_v2 + matrix + ['--filter', '@bad/no_nodes.txt']),
        list_not_numbers=(gfa + no_graph_steps + matrix + ['--filter', '@bad/words.txt'], gfa + no_graph_steps_v2 + matrix + ['--filter', '@bad/words.txt']),
        filter_keeps_none=(gfa + no_graph_steps + r2 + ['--filter', 'minor:0.4999', '--filter', '@bad/one_fixed_node.txt'],
                           gfa + no_graph_steps_v2 + r2_v2 + ['--filter', 'minor:0.4999', '--filter', '@bad/one_fixed_node.txt']),
        duplicate_samples=(gfa + no_graph_steps + matrix + ['--samples', 'bad/duplicate_samples.txt'], gfa + no_graph_steps_v2 + matrix + ['--samples', 'bad/duplicate_samples.txt']),
        empty_samples=(gfa + no_graph_steps + matrix + ['--samples', 'bad/empty_samples.txt'], gfa + no_graph_steps_v2 + matrix + ['--samples', 'bad/empty_samples.txt']),
        unknown_samples=(gfa + no_graph_steps + matrix + ['--samples', 'bad/unknown_samples.txt'], gfa + no_graph_steps_v2 + matrix + ['--samples', 'bad/unknown_samples.txt']),
        pairs_too_few=(gfa + no_graph_steps + matrix + ['--min-pairs', '1'], gfa + no_graph_steps_v2 + matrix + ['--min-pairs', '1']),
        negative_allele_count=(gfa + no_graph_steps + matrix + ['--min-allele-count', '-1'], gfa + no_graph_steps_v2 + matrix + ['--min-allele-count', '-1']),
        condition_without_tree=(gfa + no_graph_steps + matrix[:] + ['--condition-on-parents'], gfa + no_graph_steps_v2 + matrix + ['--condition-on-parents']),
        condition_and_complete=(gfa + ['--no-node-coordinates', '--assume-complete', '--condition-on-parents'] + matrix,
                                gfa + ['--no-node-coordinates', '--assume-complete', '--condition-on-parents'] + matrix),
        collapse_without_condition=(gfa + no_graph_steps + matrix + ['--collapse-backbone-parents'], gfa + no_graph_steps_v2 + matrix + ['--collapse-backbone-parents']),
        broken_zero=(gfa + no_graph_steps + matrix + ['--broken-parent-min-haplotypes', '0'], gfa + no_graph_steps_v2 + matrix + ['--broken-parent-min-haplotypes', '0']),
        drop_presence=(gfa + no_graph_steps + matrix + ['--matrix-type', 'presence', '--drop-complements'],
                       gfa + no_graph_steps_v2 + matrix + ['--matrix-type', 'presence', '--drop-complements']),
        split_without_reference=(['--gfa', 'small/no_reference.gfa'] + no_graph_steps + matrix + ['--split-inversions'],
                                 ['--gfa', 'small/no_reference.gfa'] + no_graph_steps_v2 + matrix + ['--split-inversions']),
        blocks_without_coordinates=(gfa + no_graph_steps + r2 + ['--blocks'], gfa + no_graph_steps_v2 + r2_v2 + ['--blocks', 'gabriel']),
        blocks_without_output=(gfa + no_graph_steps + matrix + ['--blocks'], gfa + no_graph_steps_v2 + matrix + ['--blocks', 'gabriel']),
        blocks_unplaced=(gfa + ['--no-conditioning-tree', '--node-coordinates', 'OUT/c', '--assume-complete'] + r2 + ['--blocks', '--filter', '@bad/unplaced_nodes.txt'],
                         gfa + ['--no-node-parents', '--node-coordinates', 'OUT/c.npz', '--assume-complete'] + r2_v2 + ['--blocks', 'gabriel', '--filter', '@bad/unplaced_nodes.txt']),
        blocks_presence=(gfa + ['--no-conditioning-tree', '--node-coordinates', 'OUT/c', '--assume-complete', '--matrix-type', 'presence'] + r2 + ['--blocks', '--filter', '@nodes.txt'],
                         gfa + ['--no-node-parents', '--node-coordinates', 'OUT/c.npz', '--assume-complete', '--matrix-type', 'presence'] + r2_v2 + ['--blocks', 'gabriel', '--filter', '@nodes.txt']),
        no_reference=(['--gfa', 'small/no_reference.gfa', '--no-conditioning-tree', '--assume-complete', '--node-coordinates', 'OUT/c'] + matrix,
                      ['--gfa', 'small/no_reference.gfa', '--no-node-parents', '--assume-complete', '--node-coordinates', 'OUT/c.npz'] + matrix),
        bad_reference_source=(gfa + ['--no-conditioning-tree', '--assume-complete', '--coordinate-reference-source', 'nope'] + matrix,
                              gfa + ['--no-node-parents', '--assume-complete', '--coordinate-reference-source', 'nope'] + matrix),
        no_backbone=(gfa + ['--no-node-coordinates', '--assume-complete', '--backbone', 'NOBODY', '--conditioning-tree', 'OUT/t'] + matrix,
                     gfa + ['--no-node-coordinates', '--assume-complete', '--backbone', 'NOBODY', '--node-parents', 'OUT/t.npz'] + matrix),
        two_snarl_sources=(gfa + ['--no-node-coordinates', '--assume-complete', '--snarls', 'graph.snarls.json', '--vg-graph', 'graph.gbz'] + matrix,
                           gfa + ['--no-node-coordinates', '--assume-complete', '--snarls', 'graph.snarls.json', '--vg-graph', 'graph.gbz'] + matrix),
        no_coordinates_and_cache=(gfa + no_graph_steps + ['--node-coordinates', 'OUT/c'] + matrix, gfa + no_graph_steps_v2 + ['--node-coordinates', 'OUT/c.npz'] + matrix),
        no_tree_and_input=(gfa + no_graph_steps + ['--snarls', 'graph.snarls.json'] + matrix, gfa + no_graph_steps_v2 + ['--snarls', 'graph.snarls.json'] + matrix),
        snarls_bad_json=(gfa + ['--no-node-coordinates', '--assume-complete', '--snarls', 'bad/snarls_not_json.txt', '--snarls-format', 'json', '--conditioning-tree', 'OUT/t'] + matrix,
                         gfa + ['--no-node-coordinates', '--assume-complete', '--snarls', 'bad/snarls_not_json.txt', '--snarls-format', 'json', '--node-parents', 'OUT/t.npz'] + matrix),
        snarls_unknown_node=(gfa + ['--no-node-coordinates', '--assume-complete', '--snarls', 'bad/snarls_unknown_node.json', '--conditioning-tree', 'OUT/t'] + matrix,
                             gfa + ['--no-node-coordinates', '--assume-complete', '--snarls', 'bad/snarls_unknown_node.json', '--node-parents', 'OUT/t.npz'] + matrix),
        snarls_twice=(gfa + ['--no-node-coordinates', '--assume-complete', '--snarls', 'bad/snarls_twice.json', '--conditioning-tree', 'OUT/t'] + matrix,
                      gfa + ['--no-node-coordinates', '--assume-complete', '--snarls', 'bad/snarls_twice.json', '--node-parents', 'OUT/t.npz'] + matrix),
        snarls_orphan=(gfa + ['--no-node-coordinates', '--assume-complete', '--snarls', 'bad/snarls_orphan.json', '--conditioning-tree', 'OUT/t'] + matrix,
                       gfa + ['--no-node-coordinates', '--assume-complete', '--snarls', 'bad/snarls_orphan.json', '--node-parents', 'OUT/t.npz'] + matrix),
    )
    (bad / 'no_nodes.txt').write_text('999999999\n')
    (bad / 'words.txt').write_text('one two\n')
    records = json.loads((out / 'graph.json').read_text())
    (bad / 'duplicate_samples.txt').write_text('HG000\nHG000\n')
    (bad / 'empty_samples.txt').write_text('\n')
    (bad / 'unknown_samples.txt').write_text('HG000\nNOBODY\n')
    coords = json.loads((out / 'coordinates.json').read_text())[0]['arrays']
    position = np.asarray(coords['position'])
    ids = np.asarray(coords['node_id'])
    check((position == 0).any(), 'every node of the graph has a coordinate')
    (bad / 'unplaced_nodes.txt').write_text(' '.join(str(i) for i in ids[position == 0]) + '\n')
    # A node on every haplotype: no minor allele, so the frequency filter leaves nothing.
    (bad / 'one_fixed_node.txt').write_text(f"{records['ids'][0]}\n")
    snarls = (out / 'graph.snarls.json').read_text().splitlines()
    (bad / 'snarls_not_json.txt').write_text('not json\n')
    (bad / 'snarls_unknown_node.json').write_text(json.dumps(dict(start=dict(name='999999'), end=dict(name='1'))) + '\n')
    (bad / 'snarls_twice.json').write_text(snarls[0] + '\n' + snarls[0] + '\n')
    orphan = json.loads(snarls[0])
    orphan['parent'] = dict(start=dict(name='1'), end=dict(name='999998'))
    (bad / 'snarls_orphan.json').write_text(json.dumps(orphan) + '\n')

    verdicts = []
    for name, (args, v2) in lines.items():
        scratch = work / f'line_{name}'
        scratch.mkdir()
        command = [str(a).replace('OUT', str(scratch)) for a in (v2 if v2 is not None else args)]
        done = run(V2_BIN / 'glink', *command, cwd=out)
        message = done.stderr.strip().splitlines()[-1] if done.returncode and done.stderr.strip() else ''
        verdicts.append(dict(name=name, args=args, exit=done.returncode, message=message))
    dump(out / 'errors.json', dict(gfas=gfas, command_lines=verdicts))


# ---- main -----------------------------------------------------------------------------------


def write_inputs(out, seed):
    ids, names, _ = core.read_gfa(out / 'graph.gfa')
    coords, _ = coordinates.prepare_coordinates(out / 'graph.gfa', out / '_c.npz', rebuild=True)
    (out / '_c.npz').unlink()
    Path(str(out / '_c.npz') + '.lock').unlink(missing_ok=True)
    array = coords.arrays
    rng = random.Random(seed)
    placed = [int(i) for i, p, a in zip(array['node_id'], array['position'], array['ambiguous']) if p > 0 and a == 0]
    chosen = sorted(rng.sample(placed, int(len(placed) * .8)))
    (out / 'nodes.txt').write_text('\n'.join(map(str, chosen)) + '\n')
    (out / 'few.txt').write_text('\n'.join(map(str, chosen[::2])) + '\n')
    samples = sorted({n.split('#')[0] for n in names if n.endswith('#2')})[:-1]
    rng.shuffle(samples)  # the sample file's order is the column order
    (out / 'samples.txt').write_text('\n'.join(samples[:40]) + '\n')


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--out', type=Path, default=HERE / 'golden')
    p.add_argument('--sites', type=int, default=30)
    p.add_argument('--samples', type=int, default=60)
    p.add_argument('--seed', type=int, default=7)
    p.add_argument('--vg', default=os.environ.get('VG', DEFAULT_VG))
    a = p.parse_args()
    final = a.out.resolve()
    final.parent.mkdir(parents=True, exist_ok=True)
    fresh = Path(tempfile.mkdtemp(dir=final.parent, prefix=f'.{final.name}.'))
    work = Path(tempfile.mkdtemp(prefix='glink-golden-'))
    try:
        subprocess.run([sys.executable, str(HERE / 'make_gfa.py'), str(fresh / 'graph.gfa'), '--sites', str(a.sites),
                        '--samples', str(a.samples), '--seed', str(a.seed)], check=True)
        subprocess.run([sys.executable, str(HERE / 'make_gfa.py'), str(fresh / 'plain.gfa'), '--sites', str(a.sites),
                        '--samples', str(a.samples), '--seed', str(a.seed), '--plain'], check=True)
        write_inputs(fresh, a.seed)
        graph_golden(fresh, a.vg)
        coordinate_golden(fresh, work)
        tree_golden(fresh, a.vg, work)
        filter_golden(fresh)
        engine_golden(fresh)
        compare_golden(fresh, work)
        error_golden(fresh, work)
        stats = {}
        for name, spec in CASES.items():
            case_golden(name, spec, fresh, work, stats)
        check(stats.get('blocks', 0) >= 5, f'only {stats.get("blocks", 0)} cases made a multi-marker block')
        require_coverage(fresh)
        v2_root = Path(core.__file__).resolve().parents[2]
        commit = run('git', 'rev-parse', 'HEAD', cwd=v2_root).stdout.strip()
        version = run(a.vg, 'version', cwd=work).stdout.splitlines()[0].strip()
        dump(fresh / 'manifest.json', dict(v2_commit=commit, vg_version=version, sites=a.sites, samples=a.samples, seed=a.seed,
                                           cases={name: spec for name, spec in CASES.items()}))
        for leftover in fresh.rglob('*.lock'):
            leftover.unlink()
        if final.exists():
            shutil.rmtree(final)
        fresh.rename(final)
    finally:
        shutil.rmtree(fresh, ignore_errors=True)
        shutil.rmtree(work, ignore_errors=True)
    size = sum(f.stat().st_size for f in final.rglob('*') if f.is_file())
    print(f'wrote {len(CASES)} cases, {size / 1e6:.1f} MB, to {final}')


if __name__ == '__main__':
    main()
