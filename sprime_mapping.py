import pandas as pd
import pysam
import pybedtools
from collections import defaultdict


def split_archaic_fragments_from_segment(segment_snps, hap_hits):
    """
    Split a segment into fragments of contiguous archaic SNPs.

    Parameters
    ----------
    segment_snps : list[int]
        Sorted positions of all SNPs in the segment (archaic + non-archaic)
    hap_hits : set[int]
        Positions of SNPs in this haplotype that carry the archaic allele

    Returns
    -------
    List of fragments: each fragment is a list of SNP positions
    optional: List of all snps in fragments: each fragment is a list of SNP positions, including non-archaic snps
    """

    fragments = []
    current_frag = []
    for pos in segment_snps:


        if pos in hap_hits:
            current_frag.append(pos)
        else:
            if current_frag:
                fragments.append(current_frag)
                current_frag = []


    if current_frag:
        fragments.append(current_frag)

    return fragments




def map_sprime_segments(
    score_file,
    vcf_file,
    out_file=None,
    target_individuals=None,
    segment_fraction=None,
    min_snps=None,
    phased=True,
    merge_distance=0,
    only_tract_output=True,
    return_full_records=False,
    single_snp_bed=False, 
):
    """
    Takes a SPrime score output file as input and returns a bed file (or a dataframe) with introgressed tracts per individual / haplotype.
    Addtionally, some filtering functions are implemented.

    Parameters
    ----------
    score_file : str, Filepath to the SPrime output file
    vcf_file: str, filepath to the vcf-file which was used as input for SPrime
    out_file : (optional) str, filepath for writing the bed-file with the introgressed tracts per individual; if None, no output file is written and the dataframe is returned
    target_indivdiuals: list, indicating the individuals in the target panel, i.e. for which introgressed fragments shouls be detected (has to correspond to the ids in the vcf_file)
    segment_fraction: (optional) float, if specified, indicates the minimum fraction of SNPs in the target which correspond to the SNPs in the segment provided by SPrime so that it is marked as introgressed in this individual
    min_snps: (optional) int: within a segment, the minimum number of adjacent SNPs so that this region of a segment is marked as introgressed in the specific individual
    phased: (bool) if phased, the inidividual names in the final bed-file are suffixed by '_1' for the first haplotype, '_2' for the second haplotype; if not phased, tracts for both haplotypes are merged. Currently only diploid input is supported. Since SPrime takes diploid phased data as input, one usually will use Phased=True,
    default: True
    merge_distance: (optional) int, for the final merging of tracts, an int > 0 indicates that also regions with a gap should be merged in case that the distance is smaller than merge_distance in base pairs. Can also (preferably) done in a later step.
    default: 0
    only_tract_output: (optional) True, in case of false, an extra column with 'haplotype' is generated in the bedfile, otherwise the individual names are suffixed as described in phased and no extra haplotype-column is generated
    default: True
    return_full_records: (optional) False: instead of returning the bed-file, a larger dataframe with the number of snps and the segment fraction for each intrgressed region is returned. Overwrites only_tract_output.
    default: False
    """

    score = pd.read_csv(score_file, sep="\t")
    segments = defaultdict(list)
    for _, r in score.iterrows():
        segments[r["SEGMENT"]].append(r)

    segment_chrom = {k: v[0]["CHROM"] for k, v in segments.items()}
    segment_positions = {k: [r["POS"] for r in v] for k, v in segments.items()}
    segment_sizes = {k: len(v) for k, v in segments.items()}

    # Load VCF
    vcf = pysam.VariantFile(vcf_file)
    samples = list(vcf.header.samples)

    if not target_individuals:
        target_individuals = samples

    # Collect archaic SNPs per haplotype
    hap_hits = defaultdict(lambda: defaultdict(lambda: {0: set(), 1: set()}))
    for seg_id, snps in segments.items():
        for snp in snps:
            chrom = str(snp["CHROM"])
            pos = int(snp["POS"])
            intro_is_alt = snp["ALLELE"] == 1

            recs = list(vcf.fetch(chrom, pos - 1, pos))
            if not recs:
                continue
            rec = recs[0]

            for sample in samples:
                if sample in target_individuals:
                    gt = rec.samples[sample]["GT"]
                    if gt is None or None in gt:
                        continue
                    if (intro_is_alt and gt[0] == 1) or (not intro_is_alt and gt[0] == 0):
                        hap_hits[seg_id][sample][0].add(pos)
                    if (intro_is_alt and gt[1] == 1) or (not intro_is_alt and gt[1] == 0):
                        hap_hits[seg_id][sample][1].add(pos)

    # Build fragments / tracts
    records = []
    for seg_id, sample_dict in hap_hits.items():
        chrom = segment_chrom[seg_id]
        seg_size = segment_sizes[seg_id]
        seg_snps = segment_positions[seg_id]

        for sample, haps in sample_dict.items():
            for hap in (0, 1):
                hits = sorted(haps[hap])
                if not hits:
                    continue

                # segment filter: a haploitype / individual has to have at least this fraction of snps which are introgressed according to SPrime
                frac = len(hits) / seg_size
                if segment_fraction is not None and frac < segment_fraction:
                    continue


                if single_snp_bed:
                    # one record per snp - in case on wants to count snps per individual
                    for snp in hits:
                        records.append({
                            "chrom": chrom,
                            "start": snp,
                            "end": snp + 1, # +1 because of half-open intervals [)
                            "individual": sample,
                            "haplotype": hap + 1,
                            "nsnps": 1,
                            "segment_fraction": 1/seg_size,
                            "segment": seg_id,
                            "archaic_snps_present": [snp]
                        })
                    continue

                # min snps filter within segment: only regions of continguous archaic snps of the size given by min_snps are added
                if min_snps is not None:
                    fragments = split_archaic_fragments_from_segment(seg_snps, set(hits))
                    for ifrag, frag in enumerate(fragments):
                        if len(frag) < min_snps:
                            continue

                        records.append({
                            "chrom": chrom,
                            "start": min(frag),
                            "end": max(frag) +1, # +1 because of half-open intervals [)
                            "individual": sample,
                            "haplotype": hap + 1,
                            "nsnps": len(frag),
                            "segment_fraction": frac,
                            "segment": seg_id,
                            "archaic_snps_present": frag,
                        })
                else:
                    # No fragment filtering by min_snps: take full span (i.e. full segment, except for the borders, in which the first and last matching SNP is used)

                    records.append({
                        "chrom": chrom,
                        "start": min(hits),
                        "end": max(hits) +1, # +1 because of half-open intervals [)
                        "individual": sample,
                        "haplotype": hap + 1,
                        "nsnps": len(hits),
                        "segment_fraction": frac,
                        "segment": seg_id,
                        "all_snps_present": seg_snps,
                        "archaic_snps_present": hits,
                    })

    df = pd.DataFrame(records)

    if return_full_records:
        return df


    if df.empty:
        if out_file:
            if only_tract_output:
                pd.DataFrame(columns=["chrom","start","end","name"]).to_csv(
                    out_file, sep="\t", index=False, header=False
                )
            else:
                df.to_csv(out_file, sep="\t", index=False)
        return df

    # Merge overlapping fragments
    group_cols = ["individual","haplotype","chrom"] if phased else ["individual","chrom"]
    merged_records = []

    for keys, sub in df.groupby(group_cols):
        sub = sub.sort_values("start")
        bed = pybedtools.BedTool.from_dataframe(sub[["chrom","start","end"]])
        merged_bed = bed.merge(d=merge_distance)
        for iv in merged_bed:
            rec = {"chrom": iv.chrom, "start": int(iv.start), "end": int(iv.end)}
            if phased:
                individual, hap, _ = keys
                rec["individual"] = individual
                rec["haplotype"] = hap
            else:
                individual, _ = keys
                rec["individual"] = individual
            merged_records.append(rec)

    df = pd.DataFrame(merged_records)

    #write output
    if only_tract_output:
        bed_rows = []
        for _, row in df.iterrows():
            name = f"{row['individual']}_{row['haplotype']}" if phased else row["individual"]
            bed_rows.append([row["chrom"], row["start"], row["end"], name])
        df = pd.DataFrame(bed_rows)
        if out_file:
            df.to_csv(out_file, sep="\t", index=False, header=False)
    else:
        if out_file:
            df.to_csv(out_file, sep="\t", index=False)
    return df

