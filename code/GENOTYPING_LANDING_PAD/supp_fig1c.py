# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 10:08:05 2026

@author: ig8
"""

from dna_features_viewer import GraphicFeature, GraphicRecord
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec


features = [
    #GraphicFeature(0, 22, strand=1, label="Mock_CLYBL_FW", color="#d6b295"),
    GraphicFeature(468, 489, strand=1, label="IG59", color="#f8d3a9"),
    #GraphicFeature(500, 520, strand=1, label="VM19", color="#f7977a"),
    #GraphicFeature(513, 533, strand=1, label="VM17", color="#f7977a"),
    GraphicFeature(602, 1434, strand=1, label="CLYBL homology arm", color="#c6c9d1"),
    #GraphicFeature(1413, 1434, strand=1, label="IG53", color="#c6c9d1"),
    GraphicFeature(1442, 1476, strand=1, label="loxP", color="#31849b"),
    GraphicFeature(1476, 1502, strand=1, label="SA", color="#b4abac"),
    GraphicFeature(1525, 1579, strand=1, label="T2A", color="#ff9ccd"),
    GraphicFeature(1582, 2374, strand=1, label="NeoR/KanR", color="#b7e6d7"),
    #GraphicFeature(1624, 1644, strand=-1, label="VM20", color="#f7977a"),
    GraphicFeature(2375, 2497, strand=1, label="SV40 poly(A) signal", color="#b4abac"),
    GraphicFeature(2498, 2532, strand=1, label="loxP", color="#31849b"),
    GraphicFeature(2534, 2914, strand=1, label="CMV enhancer", color="#c6c9d1"),
    GraphicFeature(2915, 3193, strand=1, label="chicken β-actin promoter", color="#c6c9d1"),
    GraphicFeature(3194, 4212, strand=1, label="chimeric intron", color="#a187be"),
    #GraphicFeature(3233, 3284, strand=1, label="Illumina", color="#9eafd2"),
    #GraphicFeature(3926, 3946, strand=1, label="VM100", color="#75c6a9"),
    #GraphicFeature(4169, 4189, strand=1, label="CAG-sF (from Edinburgh)", color="#b4abac"),
    #GraphicFeature(4208, 4247, strand=1, label="Illumina", color="#9eafd2"),
    GraphicFeature(4217, 4271, strand=1, label="Bxb1 attR", color="#3366cc"),
    #GraphicFeature(4222, 4242, strand=1, label="P", color="#ffef86"),
    #GraphicFeature(4242, 4271, strand=1, label="right", color="#c6c9d1"),
    #GraphicFeature(4271, 4296, strand=1, label="VM77", color="#f7977a"),
    #GraphicFeature(4272, 4298, strand=1, label="SEQUENCING_DONOR_FORWARD", color="#f58a5e"),
    #GraphicFeature(4299, 4307, strand=1, label="BARCODE", color="#b4abac"),
    #GraphicFeature(4307, 4317, strand=1, label="Kozak sequence", color="#b4abac"),
    GraphicFeature(4319, 4346, strand=1, label="c-myc NLS", color="#993366"),
    GraphicFeature(4361, 4388, strand=1, label="c-myc NLS", color="#993366"),
    GraphicFeature(4403, 5117, strand=1, label="EGFP OPTIMISED", color="#b1ff67"),
    GraphicFeature(5147, 5201, strand=1, label="T2A", color="#ff9ccd"),
    #GraphicFeature(5161, 5167, strand=-1, label="AAUAAA motif", color="#ff9ccd"),
    #GraphicFeature(5204, 5210, strand=1, label="DNA BLOCK", color="#ff00ff"),
    #GraphicFeature(5219, 5239, strand=-1, label="SEQUENCING_REVERSE (=VM103)", color="#b4abac"),
    GraphicFeature(5219, 5819, strand=1, label="PuroR", color="#b7e6d7"),
    GraphicFeature(5825, 5879, strand=1, label="Bxb1 attL", color="#3366cc"),
    #GraphicFeature(5846, 5851, strand=1, label="P", color="#ffef86"),
    #GraphicFeature(5846, 5879, strand=1, label="right", color="#c7b0e3"),
    #GraphicFeature(5849, 5874, strand=-1, label="P'", color="#ffef86"),
    GraphicFeature(5879, 6104, strand=1, label="bGH poly(A) signal", color="#b4abac"),
    #GraphicFeature(6049, 6069, strand=-1, label="VM101", color="#d59687"),
    GraphicFeature(6104, 6790, strand=1, label="Short UCOE", color="#ffef86"),
    #GraphicFeature(6694, 6714, strand=1, label="VM21", color="#f7977a"),
    GraphicFeature(6790, 7588, strand=1, label="CLYBL homology arm", color="#c6c9d1"),
    GraphicFeature(7939, 7959, strand=-1, label="IG60", color="#c6c9d1"),
    #GraphicFeature(7970, 7992, strand=-1, label="VM22", color="#f7977a"),
    #GraphicFeature(8476, 8492, strand=-1, label="Mock_CLYBL_REV", color="#f8d3a9"),
]

record = GraphicRecord(
    sequence_length=8492,
    features=features,
)


fig, ax = plt.subplots(figsize=(18, 3.5))

record.plot(
    ax=ax,
    #feature_level_height=0.4,
    plot_sequence=False,
)

plt.tight_layout()
plt.savefig("CLYBL_egfp_map.png", dpi=600, bbox_inches="tight")
plt.show()



# Load samtools depth TSV (position, coverage)
# If your file has no header:
depth = pd.read_csv(
    r".\egfp_7k_20_coverage.tsv",
    sep="\t",
    names=["pos", "depth"]
)

# If your file has a header instead, use:
# depth = pd.read_csv("coverage.tsv", sep="\t")

fig = plt.figure(figsize=(18, 5))
gs = GridSpec(2, 1, height_ratios=[1, 2], hspace=0.05)

ax_cov = fig.add_subplot(gs[0])
ax_map = fig.add_subplot(gs[1], sharex=ax_cov)

# Coverage track
ax_cov.fill_between(
    depth["pos"],
    depth["depth"],
    color="steelblue",
    alpha=0.6
)
ax_cov.plot(
    depth["pos"],
    depth["depth"],
    color="steelblue",
    linewidth=0.8
)

ax_cov.set_ylabel("Depth")
ax_cov.set_xlim(1, record.sequence_length)
ax_cov.spines["top"].set_visible(False)
ax_cov.spines["right"].set_visible(False)
ax_cov.tick_params(axis="x", labelbottom=False)

# DNA Features Viewer map
record.plot(ax=ax_map)
ax_map.set_xlim(1, record.sequence_length)

plt.tight_layout()
plt.savefig("egfp_coverage_with_features_20.png", dpi=600, bbox_inches="tight")
plt.show()
