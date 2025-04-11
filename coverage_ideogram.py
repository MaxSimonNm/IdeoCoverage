import argparse
import pandas as pd
from pyfaidx import Fasta
import pybedtools
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
import matplotlib.ticker as ticker


def get_chrom_lengths(fasta_file):
    """Get chromosome lengths from a FASTA file."""
    fasta = Fasta(fasta_file)
    return {seq.name: len(seq) for seq in fasta}


def get_bed_coverage(bed_file):
    """Get covered regions from a BED file."""
    bed = pybedtools.BedTool(bed_file).sort().merge()
    chr_covered = {}
    for feature in bed:
        chrom = feature.chrom
        start, end = int(feature.start), int(feature.end)
        chr_covered.setdefault(chrom, []).append((start, end))
    return chr_covered, bed


def get_cytoband_info(cytoband_file):
    """Get centromere and telomere positions from a cytoband file."""
    df = pd.read_csv(cytoband_file, sep='\t', header=None, names=["chr", "start", "end", "name", "gieStain"])

    # Centromeres
    centromeres = df[df["gieStain"] == "acen"]
    centromere_positions = {}
    for chrom in centromeres["chr"].unique():
        subdf = centromeres[centromeres["chr"] == chrom]
        centromere_start = subdf["start"].min()
        centromere_end = subdf["end"].max()
        centromere_positions[chrom] = (centromere_start, centromere_end)

    # Telomeres
    telomere_positions = {}
    for chrom in df["chr"].unique():
        subdf = df[df["chr"] == chrom]
        if subdf.empty:
            continue
        start_tel = (subdf.iloc[0]["start"], subdf.iloc[0]["end"])
        end_tel = (subdf.iloc[-1]["start"], subdf.iloc[-1]["end"])
        telomere_positions[chrom] = (start_tel, end_tel)

    return centromere_positions, telomere_positions


def check_overlap(bedtool, region):
    """Check if any BED region overlaps a given (chrom, start, end) region."""
    chrom, start, end = region
    region_bed = pybedtools.BedTool(f"{chrom}\t{start}\t{end}\n", from_string=True)
    return len(bedtool.intersect(region_bed, u=True)) > 0


def plot_ideogram(chrom_lengths, bed_coverage_dict, bedtool, centromeres, telomeres, output_file):
    """Plot the ideogram with coverage, centromeres, and telomeres."""
    chromosomes = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    chrom_width = 0.4
    spacing = 1.0

    fig, ax = plt.subplots(figsize=(12, 12))

    for idx, chrom in enumerate(chromosomes):
        if chrom not in chrom_lengths:
            continue
        x = idx * spacing
        chr_len = chrom_lengths[chrom]

        # Draw chromosome body
        ax.add_patch(Rectangle((x, 0), chrom_width, chr_len, color='lightgray'))

        # Draw covered regions
        for start, end in bed_coverage_dict.get(chrom, []):
            ax.add_patch(Rectangle((x, start), chrom_width, end - start, color='navy'))

        # Draw centromere
        if chrom in centromeres:
            cstart, cend = centromeres[chrom]
            covered = check_overlap(bedtool, (chrom, cstart, cend))
            centromere_color = 'red' if covered else '#ffcccc'  # lighter red if uncovered
            ax.add_patch(Rectangle((x, cstart), chrom_width, cend - cstart, color=centromere_color))

        # Draw telomeres
        if chrom in telomeres:
            (tel_start_start, tel_start_end), (tel_end_start, tel_end_end) = telomeres[chrom]
            # Start telomere
            covered_start = check_overlap(bedtool, (chrom, tel_start_start, tel_start_end))
            color_start = 'green' if covered_start else '#ccffcc'  # lighter green if uncovered
            ax.add_patch(Rectangle((x, tel_start_start), chrom_width, tel_start_end - tel_start_start, color=color_start))
            # End telomere
            covered_end = check_overlap(bedtool, (chrom, tel_end_start, tel_end_end))
            color_end = 'green' if covered_end else '#ccffcc'
            ax.add_patch(Rectangle((x, tel_end_start), chrom_width, tel_end_end - tel_end_start, color=color_end))

        # Chromosome label
        ax.text(x + chrom_width / 2, chr_len + 1e7, chrom, ha='center', va='bottom', fontsize=9)

    # Y-axis: genomic coordinate ticks
    ax.set_ylim(0, max(chrom_lengths.values()) * 1.1)
    ax.set_xlim(-0.5, len(chromosomes))
    ax.invert_yaxis()  # Flip ideogram
    ax.set_xticks([])
    ax.set_ylabel('Genomic Position (bp)')
    ax.set_title('BED Coverage Ideogram with Centromeres and Telomeres')

    # Add minor ticks every 50Mb
    ax.yaxis.set_major_locator(ticker.MultipleLocator(50_000_000))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(10_000_000))

    # Add legend
    legend_elements = [
        Patch(facecolor='lightgray', edgecolor='black', label='Chromosome Body'),
        Patch(facecolor='navy', label='Covered Region'),
        Patch(facecolor='red', label='Covered Centromere'),
        Patch(facecolor='#ffcccc', label='Uncovered Centromere'),
        Patch(facecolor='green', label='Covered Telomere'),
        Patch(facecolor='#ccffcc', label='Uncovered Telomere')
    ]
    ax.legend(handles=legend_elements, loc='lower left', bbox_to_anchor=(1, 0.5))

    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    print(f"✅ Saved plot to {output_file}")


def main():
    """Main function to parse arguments and create the ideogram plot."""
    parser = argparse.ArgumentParser(description="Create BED coverage ideogram with centromeres and telomeres.")
    parser.add_argument("--fasta", required=True, help="Reference FASTA file (indexed with .fai or use pyfaidx).")
    parser.add_argument("--bed", required=True, help="BED file with regions of interest.")
    parser.add_argument("--cytoband", required=True, help="Cytoband file (e.g., UCSC cytoBand.txt).")
    parser.add_argument("--output", required=True, help="Output image file (e.g., output.svg).")
    args = parser.parse_args()

    chrom_lengths = get_chrom_lengths(args.fasta)
    bed_coverage_dict, bedtool = get_bed_coverage(args.bed)
    centromeres, telomeres = get_cytoband_info(args.cytoband)

    plot_ideogram(chrom_lengths, bed_coverage_dict, bedtool, centromeres, telomeres, args.output)


if __name__ == "__main__":
    main()
