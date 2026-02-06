import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import numpy as np
import pandas as pd


def plot_introgression_fragments(
    df,
    chromosome=None,
    individuals_to_plot=None,
    segments_to_plot=None,
    genome_length=None,
    hap_height=1.2,                 # vertical height of haplotype bars
    inter_hap_spacing=0.1,          # spacing between haplotypes of same individual
    inter_ind_spacing=0.4,          # extra spacing between individuals
    snp_marker_size=18,             
    run_line_width=10,              
    fragment_color='red',
    non_fragment_color='blue',
    uninformative_color='green',
    genome_bg_color='lightgray',
    label_fragments=True,
    plot_uninformative_snps=False,
    true_tracts_bed=None, 
    true_tract_color='mediumseagreen',  
    true_tract_offset=0.25,
    title_fontsize=18,              
    axis_label_fontsize=12,         
    tick_label_fontsize=12,
    title=None          
):
    """
    Plot introgressed fragments with stacked haplotypes, individual labels centered,
    and adjustable spacing.
    """

    plot_df = df.copy()

    # Apply filters
    if chromosome is not None:
        plot_df = plot_df[plot_df['chrom'].astype(int) == int(chromosome)]
    if individuals_to_plot is not None:
        plot_df = plot_df[plot_df['individual'].isin(individuals_to_plot)]
    if segments_to_plot is not None:
        plot_df = plot_df[plot_df['segment'].isin(segments_to_plot)]

    if plot_df.empty:
        raise ValueError("No data to plot after filtering")

    individuals = sorted(plot_df['individual'].unique())

    # Compute base y-positions
    y_base = {}
    current_y = 0
    for ind in individuals:
        y_base[ind] = current_y
        # Each individual has 2 haplotypes stacked with small spacing
        current_y += 2 * hap_height + inter_ind_spacing

    fig, ax = plt.subplots(
        figsize=(15, current_y * 0.6)
    )

    # Genome-wide plotting range
    genome_min = 0
    if genome_length is not None:
        genome_max = genome_length
    else:
        genome_max = max(max(s) for s in plot_df['segment_snps'])

    # Plot fragments and SNPs
    for _, row in plot_df.iterrows():
        # haplotype offset within individual
        hap_offset = (row['haplotype'] - 1) * (hap_height + inter_hap_spacing)
        y = y_base[row['individual']] + hap_offset

        # Background genome line
        ax.hlines(
            y=y,
            xmin=genome_min,
            xmax=genome_max,
            color=genome_bg_color,
            linewidth=run_line_width,
            zorder=0
        )

        # Fragment SNPs
        for pos in row['fragment_positions']:
            ax.plot(
                pos, y,
                marker='|',
                markersize=snp_marker_size,
                markeredgewidth=2.5,
                color=fragment_color,
                zorder=2
            )

        # Non-fragment SNPs
        non_fragment_snps = [
            pos for pos in row['segment_snps']
            if pos not in row['fragment_positions']
        ]
        for pos in non_fragment_snps:
            ax.plot(
                pos, y,
                marker='|',
                markersize=snp_marker_size,
                markeredgewidth=2.0,
                color=non_fragment_color,
                zorder=1.5
            )

        # Uninformative SNPs
        if plot_uninformative_snps and 'uninformative_snps' in row:
            unin_only = [
                pos for pos in row['uninformative_snps']
                if pos not in row['fragment_positions']
            ]
            for pos in unin_only:
                ax.plot(
                    pos, y,
                    marker='|',
                    markersize=snp_marker_size,
                    markeredgewidth=1.8,
                    color=uninformative_color,
                    zorder=1
                )

        # Fragment label
        if label_fragments and row['fragment_positions']:
            frag_center = (
                min(row['fragment_positions']) +
                max(row['fragment_positions'])
            ) / 2
            ax.text(
                frag_center,
                y + hap_height * 0.45,
                str(row['segment']),
                ha='center',
                va='bottom',
                fontsize=10,
                zorder=3
            )

    # Plot true tracts
    if true_tracts_bed is not None:
        bed_df = true_tracts_bed.copy()
        if chromosome is not None and 'chromosome' in bed_df.columns:
            bed_df = bed_df[bed_df['chromosome'] == chromosome]

        for _, row in bed_df.iterrows():
            ind, hap = row['individual_haplotype'].rsplit('_', 1)
            hap = int(hap)

            if ind not in y_base:
                continue

            hap_offset = (hap - 1) * (hap_height + inter_hap_spacing)
            y = y_base[ind] + hap_offset + true_tract_offset

            ax.hlines(
                y=y,
                xmin=row['start'],
                xmax=row['end'],
                color=true_tract_color,
                linewidth=run_line_width *1.5, #* 0.8,
                alpha=0.9,
                zorder=4
            )

    # Create y-ticks for haplotypes
    yticks = []
    ylabels_hap = []
    for ind in individuals:
        for hap in [1, 2]:
            hap_y = y_base[ind] + (hap - 1) * (hap_height + inter_hap_spacing)
            yticks.append(hap_y)
            ylabels_hap.append(f"{hap}")

    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels_hap, fontsize=tick_label_fontsize)

    # Add individual labels centered between haplotypes
    for ind in individuals:
        center_y = y_base[ind] + (hap_height + inter_hap_spacing)/2
        ax.text(
            -0.02 * genome_max,  # slight left offset outside the genome
            center_y,
            ind,
            ha='right',
            va='center',
            fontsize=axis_label_fontsize + 1,
            #fontweight='bold'
        )

    # Axes formatting
    ax.set_xlim(genome_min, genome_max)
    ax.set_xlabel("Genome position", fontsize=axis_label_fontsize)
    ax.set_ylabel("Haplotypes", fontsize=axis_label_fontsize)
    ax.tick_params(axis='x', labelsize=tick_label_fontsize)
    
    if not title:
        title = "Introgressed fragments"
    
    ax.set_title(
        title
        + (f" on chromosome {chromosome}" if chromosome else ""),
        fontsize=title_fontsize
    )

    plt.tight_layout()
    plt.show()


def plot_introgression_full_genome(
    df,
    chromosome=None,
    individuals_to_plot=None,
    segments_to_plot=None,
    hap_height=1.0,
    snp_marker_size=10,
    run_line_width=6,
    intro_color='red',
    non_intro_color='blue',
    uninformative_color='green',
    genome_bg_color='lightgray',
    exclude_empty_segments=False,
    introgressed_fraction_threshold=None,
    label_segments=True,
    plot_uninformative_snps=True,
    true_tracts_bed=None, 
    true_tract_color='orange',  
    true_tract_offset=0.3  # vertical offset relative to hap_height
):
    """
    Plot haplotypes across a chromosome/genome showing:
    - Gray line for genome outside segments
    - SPrime SNPs colored by introgressed/non-introgressed
    - Optional uninformative (VCF-only) SNPs colored differently (green by default)
    - Optional labels above segments showing segment ID and introgressed fraction
    - Optional BED of true tracts plotted above/below haplotypes
    """

    plot_df = df.copy()
    
    # Apply filters to df
    if chromosome is not None:
        plot_df = plot_df[plot_df['chrom'] == chromosome]
    if individuals_to_plot is not None:
        plot_df = plot_df[plot_df['individual'].isin(individuals_to_plot)]
    if segments_to_plot is not None:
        plot_df = plot_df[plot_df['segment'].isin(segments_to_plot)]
    if exclude_empty_segments:
        plot_df = plot_df[plot_df['n_introgressed'] > 0]
    if introgressed_fraction_threshold is not None:
        plot_df = plot_df[plot_df['introgressed_fraction'] >= introgressed_fraction_threshold]
    
    if plot_df.empty:
        raise ValueError("No data to plot after filtering")
    
    fig, ax = plt.subplots(figsize=(15, len(plot_df['individual'].unique())*1.5))
    
    individuals = sorted(plot_df['individual'].unique())
    y_base = {ind: i*2*hap_height for i, ind in enumerate(individuals)}  # two haplines per individual
    
    # Genome-wide range for plotting background
    genome_min = plot_df['segment_snps'].apply(min).min()
    genome_max = plot_df['segment_snps'].apply(max).max()
    
    # Plot haplotypes from df
    for _, row in plot_df.iterrows():
        y = y_base[row['individual']] + (row['haplotype']-1 if row['haplotype'] else 0)*hap_height
        
        # Genome background line
        ax.hlines(y=y, xmin=genome_min, xmax=genome_max, color=genome_bg_color, linewidth=run_line_width, zorder=0)
        
        # SPrime SNPs
        for pos, f in zip(row['segment_snps'], row['introgressed_flags']):
            color = intro_color if f==1 else non_intro_color
            ax.plot(pos, y, marker='|', markersize=snp_marker_size, color=color, zorder=1)
        
        # Uninformative SNPs
        if plot_uninformative_snps and 'uninformative_snps' in row:
            for pos in row['uninformative_snps']:
                ax.plot(pos, y, marker='|', markersize=snp_marker_size, color=uninformative_color, zorder=0.5)
        
        # Segment label
        if label_segments:
            seg_center = (min(row['segment_snps']) + max(row['segment_snps'])) / 2
            label = f"{row['segment']} ({row['introgressed_fraction']:.2f})"
            ax.text(seg_center, y + hap_height*0.4, label, ha='center', va='bottom', fontsize=8, zorder=2)
    
    # Plot true tracts if provided
    if true_tracts_bed is not None:
        bed_df = true_tracts_bed.copy()
        if chromosome is not None and 'chromosome' in bed_df.columns:
            bed_df = bed_df[bed_df['chromosome'] == chromosome]
        for _, row in bed_df.iterrows():
            ind, hap = row['individual_haplotype'].split('_')  # e.g., "tsk_50_1"
            hap = int(hap)
            if ind not in y_base:
                continue
            y = y_base[ind] + (hap-1)*hap_height + true_tract_offset  # slightly above
            ax.hlines(y=y, xmin=row['start'], xmax=row['end'], color=true_tract_color,
                      linewidth=run_line_width*0.6, alpha=0.7, zorder=3)
    
    ax.set_xlabel("Genome position")
    ax.set_ylabel("Individual (two lines = two haplotypes)")
    ax.set_yticks([y_base[ind]+hap_height/2 for ind in individuals])
    ax.set_yticklabels(individuals)
    
    title = "Introgressed fragments"
    if chromosome is not None:
        title += f" on chromosome {chromosome}"
    ax.set_title(title)
    
    plt.tight_layout()
    plt.show()





def plot_introgression_structure_by_chrom(
    df,
    chromosome=None,
    individuals_to_plot=None,
    segments_to_plot=None,
    hap_height=1.0,
    snp_marker_size=10,
    run_line_width=6,
    intro_color='red',
    non_intro_color='blue',
    segment_line_color='lightgray',
    segment_line_padding=1
):
    """
    Plot introgressed and non-introgressed SNPs along genome segments,
    one line per haplotype (two per individual for phased data).
    
    Parameters
    ----------
    df : pd.DataFrame
        Output from map_sprime_segments_with_full_structure()
    chromosome : str, optional
        Only plot this chromosome. If None, plot all chromosomes.
    individuals_to_plot : list[str], optional
        Only plot these individuals. If None, all are plotted.
    segments_to_plot : list[int], optional
        Only plot these segments. If None, all are plotted.
    hap_height : float
        Vertical spacing between haplotype lines.
    snp_marker_size : int
        Marker size for SNPs.
    run_line_width : float
        Width of segment line.
    intro_color : str
        Color for introgressed SNPs.
    non_intro_color : str
        Color for non-introgressed SNPs.
    segment_line_color : str
        Color for background segment line.
    segment_line_padding : int
        Extra padding (bps) on each side of segment line to make it more visible.
    """
    
    plot_df = df.copy()
    
    # Filter chromosome
    if chromosome is not None:
        plot_df = plot_df[plot_df['chrom'] == chromosome]
    
    # Filter individuals and segments
    if individuals_to_plot is not None:
        plot_df = plot_df[plot_df['individual'].isin(individuals_to_plot)]
    if segments_to_plot is not None:
        plot_df = plot_df[plot_df['segment'].isin(segments_to_plot)]
    
    if plot_df.empty:
        raise ValueError("No data to plot after filtering.")
    
    fig, ax = plt.subplots(figsize=(15, len(plot_df['individual'].unique()) * 1.5))
    
    # Assign y positions
    individuals = sorted(plot_df['individual'].unique())
    y_base = {ind: i*2*hap_height for i, ind in enumerate(individuals)}  # two lines per individual
    
    for _, row in plot_df.iterrows():
        y = y_base[row['individual']] + (row['haplotype']-1 if row['haplotype'] else 0)*hap_height
        snps = row['segment_snps']
        flags = row['introgressed_flags']
        
        # Draw segment line (with optional padding)
        xmin = min(snps) - segment_line_padding
        xmax = max(snps) + segment_line_padding
        ax.hlines(
            y=y,
            xmin=xmin,
            xmax=xmax,
            color=segment_line_color,
            linewidth=run_line_width,
            zorder=0
        )
        
        # Draw SNPs on top
        for pos, f in zip(snps, flags):
            color = intro_color if f==1 else non_intro_color
            ax.plot(pos, y, marker='|', markersize=snp_marker_size, color=color, zorder=1)
    
    ax.set_xlabel("Genome position")
    ax.set_ylabel("Individuals (two lines = two haplotypes)")
    ax.set_yticks([y_base[ind]+hap_height/2 for ind in individuals])
    ax.set_yticklabels(individuals)
    title = f"Introgressed fragments"
    if chromosome is not None:
        title += f" on chromosome {chromosome}"
    ax.set_title(title)
    
    plt.tight_layout()
    plt.show()
