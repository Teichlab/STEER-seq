# -*- coding: utf-8 -*-
"""
Created on Sun Jul  5 11:26:44 2026

@author: ig8
"""

from dna_features_viewer import GraphicFeature, GraphicRecord

features = [
    #GraphicFeature(start=0, end=22, strand=1, label="Mock_CLYBL_FW", color="#d6b295"),
    GraphicFeature(start=468, end=489, strand=1, label="IG59", color="#f8d3a9"),
    GraphicFeature(start=500, end=520, strand=1, label="VM19", color="#f7977a"),
    GraphicFeature(start=513, end=533, strand=1, label="VM17", color="#f7977a"),
    GraphicFeature(start=602, end=1434, strand=-1, label="HA-L CLYBL", color="#d59687"),
    #GraphicFeature(start=1144, end=1164, strand=1, label="VM38", color="#f7977a"),
    #GraphicFeature(start=1272, end=1295, strand=1, label="IG57", color="#d6b295"),
    #GraphicFeature(start=1366, end=1437, strand=1, label="Homology_Arm", color="#b1ff67"),
    GraphicFeature(start=1413, end=1434, strand=1, label="IG53", color="#c6c9d1"),
    GraphicFeature(start=1434, end=1454, strand=-1, label="hCLYBLmutant_sgRNA", color="#b4abac"),
    GraphicFeature(start=1434, end=1454, strand=-1, label="hCLYBL_sgRNA", color="#9eafd2"),
    GraphicFeature(start=1436, end=2234, strand=-1, label="HA-R CLYBL", color="#b7e6d7"),
    #GraphicFeature(start=1437, end=1508, strand=1, label="Homology_Arm", color="#9eafd2"),
    GraphicFeature(start=1566, end=1589, strand=-1, label="IG58", color="#d59687"),
    GraphicFeature(start=1732, end=1752, strand=-1, label="VM39", color="#f7977a"),
    GraphicFeature(start=2585, end=2605, strand=-1, label="IG60", color="#c6c9d1"),
    GraphicFeature(start=2616, end=2638, strand=-1, label="VM22", color="#f7977a"),
    #GraphicFeature(start=3122, end=3138, strand=-1, label="Mock_CLYBL_REV", color="#f8d3a9"),
]

record = GraphicRecord(
    sequence_length=3138,
    features=features,
)
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(14, 3))
record.plot(
    ax=ax,
    #feature_level_height=0.4,
    plot_sequence=False,
)

plt.tight_layout()
plt.savefig("CLYBL_record.png", dpi=600, bbox_inches="tight")
plt.show()



from dna_features_viewer import GraphicFeature, GraphicRecord
import matplotlib.pyplot as plt

features = [
    #GraphicFeature(0, 22, strand=1, label="Mock_CLYBL_FW", color="#d6b295"),
    GraphicFeature(271, 291, strand=1, label="VM85", color="#f7977a"),
    GraphicFeature(468, 489, strand=1, label="IG59", color="#f8d3a9"),
    GraphicFeature(500, 520, strand=1, label="VM19", color="#f7977a"),
    GraphicFeature(513, 533, strand=1, label="VM17", color="#f7977a"),

    GraphicFeature(602, 1434, strand=1, label="CLYBL homology arm", color="#c6c9d1"),
    GraphicFeature(1144, 1164, strand=1, label="VM38", color="#f7977a"),

    GraphicFeature(1442, 1476, strand=1, label="loxP", color="#31849b"),
    GraphicFeature(1476, 1502, strand=1, label="SA", color="#b4abac"),

    GraphicFeature(1525, 1579, strand=1, label="T2A", color="#ff9ccd"),
    GraphicFeature(1582, 2374, strand=1, label="NeoR/KanR", color="#b7e6d7"),
    GraphicFeature(1591, 1611, strand=-1, label="VM87", color="#f7977a"),
    GraphicFeature(1624, 1644, strand=-1, label="VM20", color="#f7977a"),

    GraphicFeature(2375, 2497, strand=1, label="SV40 poly(A)", color="#b4abac"),
    GraphicFeature(2498, 2532, strand=1, label="loxP", color="#31849b"),

    GraphicFeature(2534, 2914, strand=1, label="CMV enhancer", color="#c6c9d1"),
    GraphicFeature(2915, 3193, strand=1, label="CAG promoter", color="#c6c9d1"),
    #GraphicFeature(3166, 3172, strand=1, label="TATA", color="#f58a5e"),
    GraphicFeature(3194, 4212, strand=1, label="Chimeric intron", color="#a187be"),

    GraphicFeature(3926, 3946, strand=1, label="VM100", color="#75c6a9"),

    GraphicFeature(4217, 4275, strand=1, label="Bxb1 attP", color="#3366cc"),
    #GraphicFeature(4222, 4247, strand=1, label="P", color="#ffef86"),
    #GraphicFeature(4245, 4270, strand=-1, label="P'", color="#ffef86"),

    GraphicFeature(4279, 4286, strand=1, label="Kozak", color="#b4abac"),
    GraphicFeature(4283, 5804, strand=1, label="Bxb1 integrase", color="#b1ff67"),
    GraphicFeature(5783, 5804, strand=1, label="SV40 NLS", color="#c6c9d1"),

    GraphicFeature(5813, 5870, strand=1, label="P2A", color="#ff9ccd"),
    GraphicFeature(5876, 6296, strand=1, label="Bsd", color="#b7e6d7"),

    GraphicFeature(6311, 6365, strand=1, label="T2A", color="#ff9ccd"),
    GraphicFeature(6371, 6398, strand=1, label="c-myc NLS", color="#993366"),
    GraphicFeature(6413, 6440, strand=1, label="c-myc NLS", color="#993366"),

    #GraphicFeature(6488, 6498, strand=1, label="Kozak", color="#b4abac"),
    GraphicFeature(6494, 7205, strand=1, label="mApple", color="#f58a5e"),
    #GraphicFeature(6497, 6517, strand=1, label="VM34", color="#f7977a"),

    #GraphicFeature(7198, 7207, strand=1, label="splice?", color="#75c6a9"),
    GraphicFeature(7205, 7263, strand=1, label="Bxb1 attP", color="#3366cc"),
    #GraphicFeature(7210, 7235, strand=1, label="P", color="#ffef86"),
    #GraphicFeature(7233, 7258, strand=-1, label="P'", color="#ffef86"),

    GraphicFeature(7263, 7488, strand=1, label="bGH poly(A)", color="#b4abac"),
    GraphicFeature(7433, 7453, strand=-1, label="VM101", color="#d59687"),

    GraphicFeature(7488, 8174, strand=1, label="Short UCOE", color="#ffef86"),
    GraphicFeature(8078, 8098, strand=1, label="VM21", color="#f7977a"),

    GraphicFeature(8174, 8972, strand=1, label="CLYBL homology arm", color="#c6c9d1"),

    GraphicFeature(8470, 8490, strand=-1, label="VM39", color="#f7977a"),
    GraphicFeature(9279, 9300, strand=-1, label="VM24", color="#f7977a"),
    GraphicFeature(9323, 9343, strand=-1, label="IG60", color="#c6c9d1"),
    GraphicFeature(9354, 9376, strand=-1, label="VM22", color="#f7977a"),
    #GraphicFeature(9860, 9876, strand=-1, label="Mock_CLYBL_REV", color="#f8d3a9"),
]

record = GraphicRecord(
    sequence_length=9876,
    features=features,
)

fig, ax = plt.subplots(figsize=(18, 3.5))

record.plot(
    ax=ax,
    #feature_level_height=0.4,
    plot_sequence=False,
)

plt.tight_layout()
plt.savefig("CLYBL_donor_map.png", dpi=600, bbox_inches="tight")
plt.show()

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# Load samtools depth TSV (position, coverage)
# If your file has no header:
depth = pd.read_csv(
    ".\het_f9_lp_allele.tsv",
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
plt.savefig("coverage_with_features.png", dpi=600, bbox_inches="tight")
plt.show()



import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# Load samtools depth TSV (position, coverage)
# If your file has no header:
depth = pd.read_csv(
    ".\het_f9_wt_allele.tsv",
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
plt.savefig("coverage_with_features.png", dpi=600, bbox_inches="tight")
plt.show()