import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import numpy as np
import pandas as pd

def infer_chrom_lengths_from_vcf(vcf_path):
    """
    Infer chromosome lengths from a VCF file.


    Returns
    -------
    chrom_lengths : dict
        {chrom: length}
    source : str
        "contig" or "max_snp"
    """
    from collections import defaultdict
    try:
        import pysam
    except ImportError:
        raise ImportError("pysam is required to read chromosome lengths from a VCF")

    vcf = pysam.VariantFile(vcf_path)

    # look in vcf contig lengths
    chrom_lengths = {
        contig: vcf.header.contigs[contig].length
        for contig in vcf.header.contigs
        if vcf.header.contigs[contig].length is not None
    }

    if chrom_lengths:
        return chrom_lengths, "contig"

    # fallback if contigs in vcf-file not provided: max SNP position
    max_pos = defaultdict(int)
    for rec in vcf.fetch():
        max_pos[rec.chrom] = max(max_pos[rec.chrom], rec.pos)

    if not max_pos:
        raise ValueError("VCF contains no variant records")

    return dict(max_pos), "max_snp"


def plot_introgression_support_exact(
    bed_df,
    chromosomes=None,
    count_mode="haplotype",       
    cmap="Reds",
    vcf_file=None,
    haplotype=None,               
    backbone_linewidth=10,
    introgression_linewidth=16,
    title=None
):
    """
    Plot exact introgression support over all individuals / haplotypes and chromosomes

    """

    df = bed_df.copy()
    df[['individual', 'haplotype']] = (
        df['individual_haplotype'].str.rsplit('_', n=1, expand=True)
    )

    # --- filter by haplotype if requested ---
    if haplotype is not None:
        if str(haplotype) not in ["1", "2"]:
            raise ValueError("haplotype must be 1 or 2")
        df = df[df['haplotype'] == str(haplotype)]

    # ---- filter chromosomes if provided ----
    if chromosomes is not None:
        df = df[df['chrom'].isin(chromosomes)]

    chroms = sorted(df['chrom'].unique())

    # ---- infer chromosome lengths from VCF if provided ----
    chrom_lengths = None
    length_source = None

    if vcf_file is not None:
        chrom_lengths, length_source = infer_chrom_lengths_from_vcf(vcf_file)

    fig, ax = plt.subplots(figsize=(16, 0.6 * len(chroms)))

    segments = []
    max_count = 0
    chrom_ranges = {}

    for chrom in chroms:
        cdf = df[df['chrom'] == chrom]

        if chrom_lengths is not None and chrom in chrom_lengths:
            chrom_end = chrom_lengths[chrom]
        else:
            chrom_end = cdf['end'].max()

        chrom_ranges[chrom] = (0, chrom_end)

        # get all breakpoints
        breakpoints = np.sort(
            np.unique(np.concatenate([cdf['start'].values, cdf['end'].values]))
        )

        for i in range(len(breakpoints) - 1):
            start, end = breakpoints[i], breakpoints[i + 1]

            overlapping = cdf[
                (cdf['start'] < end) & (cdf['end'] > start)
            ]

            if overlapping.empty:
                continue

            if count_mode == "individual":
                count = overlapping['individual'].nunique()
            else:
                count = len(overlapping)

            segments.append((chrom, start, end, count))
            max_count = max(max_count, count)

    # color mapping
    norm = Normalize(vmin=0, vmax=max_count)
    sm = ScalarMappable(norm=norm, cmap=cmap)

    chrom_y = {chrom: i for i, chrom in enumerate(chroms)}

    for chrom, (start, end) in chrom_ranges.items():
        y = chrom_y[chrom]
        ax.hlines(
            y=y,
            xmin=start,
            xmax=end,
            color="lightgray",
            linewidth=backbone_linewidth,
            zorder=1
        )

    # ---- draw introgressed segments ----
    for chrom, start, end, count in segments:
        y = chrom_y[chrom]
        ax.hlines(
            y=y,
            xmin=start,
            xmax=end,
            color=sm.to_rgba(count),
            linewidth=introgression_linewidth,
            zorder=2
        )

    # ---- axes & labels ----
    ax.set_yticks(range(len(chroms)))
    ax.set_yticklabels(chroms, fontsize=14)

    ax.set_xlabel("Genome position", fontsize=15)
    ax.tick_params(axis="x", labelsize=13)
    ax.set_xlim(left=0)

    # ---- figure title ----
    if not title:
        title = "Introgression"
    if haplotype is not None:
        title += f" — haplotype {haplotype}"


    ax.set_title(title, fontsize=17)

    cbar = fig.colorbar(sm, ax=ax, pad=0.01)
    cbar.set_label(
        "Number of haplotypes" if count_mode == "haplotype" else "Number of individuals",
        fontsize=14
    )
    cbar.ax.tick_params(labelsize=13)

    plt.tight_layout()
    plt.show()




def plot_introgression_from_bed(
    bed_df,
    chromosome=None,
    individuals_to_plot=None,
    hap_height=1.0,
    run_line_width=6,
    fragment_color='red',
    genome_bg_color='lightgray'
):
    """
    Plot introgressed fragments directly from a BED-like file/dataframe.

    Parameters
    ----------
    bed_df : pd.DataFrame
        Columns: chrom, start, end, individual_haplotype
    chromosome : str, optional
        If provided, filter by chromosome
    individuals_to_plot : list[str], optional
        If provided, filter by individual (ignores haplotype suffix)
    hap_height : float
        Vertical spacing per haplotype
    run_line_width : float
        Width of the fragment lines
    fragment_color : str
        Color for introgressed fragments
    genome_bg_color : str
        Color for genome background lines
    """

    df = bed_df.copy()
    # split individual_haplotype into individual and hap
    df[['individual', 'haplotype']] = df['individual_haplotype'].str.rsplit('_', n=1, expand=True)
    df['haplotype'] = df['haplotype'].astype(int)

    if chromosome is not None:
        df = df[df['chrom'] == chromosome]
    
    if individuals_to_plot is not None:
        df = df[df['individual'].isin(individuals_to_plot)]
    
    if df.empty:
        raise ValueError("No data to plot after filtering")
    
    # Prepare y positions
    individuals = sorted(df['individual'].unique())
    y_base = {ind: i*2*hap_height for i, ind in enumerate(individuals)}
    
    # Genome-wide range
    genome_min = df['start'].min()
    genome_max = df['end'].max()
    
    fig, ax = plt.subplots(figsize=(15, len(individuals)*1.5))
    
    for _, row in df.iterrows():
        y = y_base[row['individual']] + (row['haplotype']-1)*hap_height
        
        # genome background line
        ax.hlines(y=y, xmin=genome_min, xmax=genome_max, color=genome_bg_color,
                  linewidth=run_line_width, zorder=0)
        
        # introgressed fragment
        ax.hlines(y=y, xmin=row['start'], xmax=row['end'], color=fragment_color,
                  linewidth=run_line_width, zorder=1)
    
    ax.set_xlabel("Genome position")
    ax.set_ylabel("Individual (haplotype 1/2)")
    ax.set_yticks([y_base[ind]+hap_height/2 for ind in individuals])
    ax.set_yticklabels(individuals)
    title = "Introgressed fragments"
    if chromosome is not None:
        title += f" on chromosome {chromosome}"
    ax.set_title(title)
    
    plt.tight_layout()
    plt.show()



def plot_introgression_bed_genome_subplots(
    bed_df,
    chromosomes=None,
    individuals_to_plot=None,
    hap_height=1.0,
    run_line_width=6,
    fragment_color='red',
    genome_bg_color='lightgray',
    chromosome_bg_colors=('white', 'whitesmoke')
):
    """
    Plot introgressed fragments with one subplot per chromosome.

    Parameters
    ----------
    bed_df : pd.DataFrame
        Columns: chrom, start, end, individual_haplotype 
    chromosomes : list[str], optional
        Only plot these chromosomes (default: all present in bed_df)
    individuals_to_plot : list[str], optional
        Filter by individual (ignores haplotype suffix)
    hap_height : float
        Vertical spacing per haplotype
    run_line_width : float
        Width of fragment lines
    fragment_color : str
        Color of introgressed fragments
    genome_bg_color : str
        Background line color
    chromosome_bg_colors : tuple
        Alternating background colors for chromosomes
    """

    df = bed_df.copy()
    df[['individual', 'haplotype']] = df['individual_haplotype'].str.rsplit('_', n=1, expand=True)
    df['haplotype'] = df['haplotype'].astype(int)

    if chromosomes is not None:
        df = df[df['chrom'].isin(chromosomes)]
    if individuals_to_plot is not None:
        df = df[df['individual'].isin(individuals_to_plot)]
    if df.empty:
        raise ValueError("No data to plot after filtering")

    individuals = sorted(df['individual'].unique())
    y_base = {ind: i*2*hap_height for i, ind in enumerate(individuals)}
    
    chroms = sorted(df['chrom'].unique())
    fig_height = max(4, len(individuals) * 1.5 * len(chroms))
    
    fig, axes = plt.subplots(len(chroms), 1, figsize=(15, fig_height), sharex=False)
    if len(chroms) == 1:
        axes = [axes]  # make iterable

    for i, (ax, chrom) in enumerate(zip(axes, chroms)):
        chrom_df = df[df['chrom'] == chrom]
        genome_min = chrom_df['start'].min()
        genome_max = chrom_df['end'].max()

        # optional alternating background for readability
        ax.set_facecolor(chromosome_bg_colors[i % len(chromosome_bg_colors)])

        for _, row in chrom_df.iterrows():
            y = y_base[row['individual']] + (row['haplotype']-1)*hap_height
            # genome background line
            ax.hlines(y=y, xmin=genome_min, xmax=genome_max, color=genome_bg_color,
                      linewidth=run_line_width, zorder=0)
            # introgressed fragment
            ax.hlines(y=y, xmin=row['start'], xmax=row['end'], color=fragment_color,
                      linewidth=run_line_width, zorder=1)

        ax.set_ylabel("Individuals")
        ax.set_yticks([y_base[ind]+hap_height/2 for ind in individuals])
        ax.set_yticklabels(individuals)
        ax.set_title(f"Chromosome {chrom}")
        ax.grid(axis='x', linestyle='--', alpha=0.3)

    axes[-1].set_xlabel("Genome position")
    plt.tight_layout()
    plt.show()



def plot_introgression_bed_genome_with_true(
    bed_df,
    true_tracts_bed=None,
    chromosomes=None,
    individuals_to_plot=None,
    hap_height=1.0,
    run_line_width=6,
    fragment_color='red',
    genome_bg_color='lightgray',
    true_tract_color='orange',
    true_tract_offset=0.3
):
    """
    Plot introgressed fragments from a BED-like file across one or multiple chromosomes,
    with optional true tracts overlay.

    Parameters
    ----------
    bed_df : pd.DataFrame
        Columns: chrom, start, end, individual_haplotype 
    true_tracts_bed : pd.DataFrame, optional
        Same columns as bed_df, plotted above the haplotype lines
    chromosomes : list[str], optional
        Only plot these chromosomes (default: all present in bed_df)
    individuals_to_plot : list[str], optional
        Filter by individual (ignores haplotype suffix)
    hap_height : float
        Vertical spacing per haplotype
    run_line_width : float
        Width of fragment lines
    fragment_color : str
        Color of introgressed fragments
    genome_bg_color : str
        Background line color
    true_tract_color : str
        Color for true tracts overlay
    true_tract_offset : float
        Vertical offset for true tracts relative to hap_height
    """

    df = bed_df.copy()
    df[['individual', 'haplotype']] = df['individual_haplotype'].str.rsplit('_', n=1, expand=True)
    df['haplotype'] = df['haplotype'].astype(int)

    if chromosomes is not None:
        df = df[df['chrom'].isin(chromosomes)]
    
    if individuals_to_plot is not None:
        df = df[df['individual'].isin(individuals_to_plot)]
    
    if df.empty:
        raise ValueError("No data to plot after filtering")
    
    individuals = sorted(df['individual'].unique())
    y_base = {ind: i*2*hap_height for i, ind in enumerate(individuals)}
    
    fig_height = max(6, len(individuals) * 1.5)
    fig, ax = plt.subplots(figsize=(15, fig_height))
    
    chroms = sorted(df['chrom'].unique())
    for chrom in chroms:
        chrom_df = df[df['chrom'] == chrom]
        genome_min = chrom_df['start'].min()
        genome_max = chrom_df['end'].max()

        # Plot called fragments
        for _, row in chrom_df.iterrows():
            y = y_base[row['individual']] + (row['haplotype']-1)*hap_height
            ax.hlines(y=y, xmin=genome_min, xmax=genome_max, color=genome_bg_color,
                      linewidth=run_line_width, zorder=0)
            # fragment
            ax.hlines(y=y, xmin=row['start'], xmax=row['end'], color=fragment_color,
                      linewidth=run_line_width, zorder=1)

    if true_tracts_bed is not None:
        bed_df_true = true_tracts_bed.copy()
        bed_df_true[['individual', 'haplotype']] = bed_df_true['individual_haplotype'].str.rsplit('_', n=1, expand=True)
        bed_df_true['haplotype'] = bed_df_true['haplotype'].astype(int)
        if chromosomes is not None:
            bed_df_true = bed_df_true[bed_df_true['chrom'].isin(chromosomes)]
        if individuals_to_plot is not None:
            bed_df_true = bed_df_true[bed_df_true['individual'].isin(individuals_to_plot)]
        
        for _, row in bed_df_true.iterrows():
            if row['individual'] not in y_base:
                continue
            y = y_base[row['individual']] + (row['haplotype']-1)*hap_height + true_tract_offset
            ax.hlines(y=y, xmin=row['start'], xmax=row['end'], color=true_tract_color,
                      linewidth=run_line_width*0.6, alpha=0.7, zorder=2)

    ax.set_xlabel("Genome position")
    ax.set_ylabel("Individual (haplotype 1/2)")
    ax.set_yticks([y_base[ind]+hap_height/2 for ind in individuals])
    ax.set_yticklabels(individuals)
    
    title = "Introgressed fragments with true tracts"
    if chromosomes is not None:
        title += " (selected chromosomes)"
    ax.set_title(title)
    
    plt.tight_layout()
    plt.show()
