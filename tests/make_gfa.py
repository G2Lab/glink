"""Write a random but realistic test GFA for comparing GlinkV2 and GlinkV3.

A chain of bubbles on a CHM13 backbone: each site has a reference allele node, an
alternative node, a deletion edge, an inversion of the reference node, or a nested bubble
(so that a node's parent can itself have a parent). Haplotypes are mosaics of a few
founders, so there is real LD and real blocks.

The awkward cases are deliberate: some haplotypes are truncated (unknown beyond their
ends), a few start inside a nested bubble (visiting a child while skipping its parent), a
few are written on the opposite strand, one sample has a single haplotype and one a third,
one third haplotype touches no reference node, one haplotype is a P line, a GRCh38 walk
must not become a column, and S lines are shuffled.

--plain leaves all of that out but the shuffling, so that every haplotype spans the graph and the
nodes on the chain are visited by everyone.

    python tests/make_gfa.py OUT.gfa [--sites 1500] [--samples 60] [--seed 1] [--plain]
"""
import argparse

import numpy as np


def main():
    p = argparse.ArgumentParser()
    p.add_argument('out')
    p.add_argument('--sites', type=int, default=1500)
    p.add_argument('--samples', type=int, default=60)
    p.add_argument('--seed', type=int, default=1)
    p.add_argument('--plain', action='store_true')
    a = p.parse_args()
    rng = np.random.default_rng(a.seed)
    n = a.sites
    kind = rng.choice(['snp', 'del', 'inv', 'multi', 'nested'], size=n, p=[.5, .18, .08, .1, .14])
    kind[n // 2] = 'nested'  # always at least one, for the haplotypes starting inside it
    kind[n // 3] = 'inv'     # and one inversion
    shared = [1 + 3 * k for k in range(n + 1)]           # S_k
    ref = [2 + 3 * k for k in range(n)]                  # R_k
    alt = [3 + 3 * k for k in range(n)]                  # A_k
    extra = {k: 10 * n + k for k in range(n) if kind[k] == 'multi'}   # a third allele
    inner = {k: tuple(20 * n + 4 * k + j for j in range(1, 4)) for k in range(n) if kind[k] == 'nested'}  # X, Y, close
    length = {}
    for node in shared:
        length[node] = int(rng.integers(20, 400))
    for k in range(n):
        length[ref[k]] = int(rng.integers(1, 60))
        length[alt[k]] = int(rng.integers(1, 60))
        if k in extra:
            length[extra[k]] = int(rng.integers(1, 60))
        for node in inner.get(k, ()):
            length[node] = int(rng.integers(1, 60))

    founders = rng.random((8, n))
    freq = rng.uniform(0.02, 0.6, size=n)
    freq[n // 3] = 0.5  # common, so that haplotypes read it both ways

    def allele(f, k):
        u = founders[f, k]
        if kind[k] in ('multi', 'nested'):  # nested: 1 and 2 pick the inner branch
            return 0 if u > freq[k] else (1 if u > freq[k] / 3 else 2)
        return int(u < freq[k])

    edges = set()
    W = []

    def walk(alleles, lo, hi):
        steps = [('>', shared[lo])]
        for k in range(lo, hi):
            al = alleles[k]
            if al == 0:
                steps.append(('>', ref[k]))
            elif kind[k] == 'del':
                pass
            elif kind[k] == 'inv':
                steps.append(('<', ref[k]))
            elif kind[k] == 'nested':
                x, y, close = inner[k]
                steps += [('>', alt[k]), ('>', x if al == 1 else y), ('>', close)]
            elif al == 2:
                steps.append(('>', extra[k]))
            else:
                steps.append(('>', alt[k]))
            steps.append(('>', shared[k + 1]))
        for (o1, a1), (o2, a2) in zip(steps, steps[1:]):
            edges.add((a1, '+' if o1 == '>' else '-', a2, '+' if o2 == '>' else '-'))
        return steps

    def text(steps):
        return ''.join(o + str(x) for o, x in steps)

    def reversed_text(steps):
        return ''.join(('<' if o == '>' else '>') + str(x) for o, x in reversed(steps))

    def link(steps):
        for (o1, a1), (o2, a2) in zip(steps, steps[1:]):
            edges.add((a1, '+' if o1 == '>' else '-', a2, '+' if o2 == '>' else '-'))

    reference = walk([0] * n, 0, n)
    total = sum(length[x] for _, x in reference)
    W.append(f'W\tCHM13\t0\tchr21\t0\t{total}\t{text(reference)}')
    W.append(f'W\tGRCh38\t0\tchr21\t0\t{total}\t{text(reference)}')  # also not a haplotype
    haps = []
    for s in range(a.samples):
        for h in (1, 2):
            if s == a.samples - 1 and h == 2 and not a.plain:
                continue  # a sample with one haplotype
            # mosaic of founders with switches
            alleles, f = [], int(rng.integers(8))
            for k in range(n):
                if rng.random() < 0.004:
                    f = int(rng.integers(8))
                al = allele(f, k)
                if rng.random() < 0.01:
                    al = 1 - min(al, 1)
                alleles.append(al)
            lo, hi = 0, n
            if rng.random() < 0.15 and not a.plain:
                lo = int(rng.integers(0, n // 3))
            if rng.random() < 0.15 and not a.plain:
                hi = int(rng.integers(2 * n // 3, n))
            steps = walk(alleles, lo, hi)
            if len(haps) in range(10, 18) and not a.plain:  # starts inside the nested bubble, skipping its parent
                k = n // 2
                alleles = alleles[:k] + [1 + len(haps) % 2] + alleles[k + 1:]
                steps = walk(alleles, k, hi)[2:]
            haps.append((f'HG{s:03d}', h, steps))
    # Third haplotypes: one on the backbone, one touching no reference node.
    if not a.plain:
        haps.append(('HG000', 3, walk([0] * n, 0, n // 3)))
        haps.append(('HG001', 3, [('>', alt[k]) for k in range(n) if kind[k] == 'snp'][:1]))
    for _, _, steps in haps:
        link(steps)
    lines = ['H\tVN:Z:1.1\tRS:Z:CHM13']
    nodes = sorted(length)
    order = rng.permutation(len(nodes))
    for i in order:
        lines.append(f'S\t{nodes[i]}\t{"ACGT" * (length[nodes[i]] // 4) + "A" * (length[nodes[i]] % 4)}')
    for a1, o1, a2, o2 in sorted(edges):
        lines.append(f'L\t{a1}\t{o1}\t{a2}\t{o2}\t0M')
    lines += W
    for i, (sample, h, steps) in enumerate(haps):
        if i == 3:
            segs = ','.join(f'{x}{"+" if o == ">" else "-"}' for o, x in steps)
            lines.append(f'P\t{sample}#{h}#ctg1\t{segs}\t*')
        else:
            size = sum(length[x] for _, x in steps)
            written = reversed_text(steps) if i % 13 == 5 and not a.plain else text(steps)  # the opposite strand
            lines.append(f'W\t{sample}\t{h}\tctg{i}\t0\t{size}\t{written}')
    with open(a.out, 'w') as handle:
        handle.write('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main()
