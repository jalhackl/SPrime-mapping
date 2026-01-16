SPrime mapping function

sprime_mapping.py contains the function which does the mapping of SPrime output to a vcf-file.
SPrime provides allele-based output: Candidate introgressed segments for the whole population can be found in the score-file. 
For every SNP in each segment, it is indicated wheter the REF or the ALT allele is considered the introgressed one (in the ALLELE column, 0-REF is introgressed,
1-ALT is introgressed).
In case of only 1 target individual in the vcf-file 'mapping' is trivial, because the prediction as introgressed was done solely on this individual.
In case of more individuals, however, it is not trivial - the SNPs of each individual have to be compared to the SNPs marked as introgressed in the 
score-file segments.

sprime_info preprocesses the obtained dataframe so that the total length of the introgressed tracts, average segments per 
individual/haplotype etc. can be easily computed.


In sprime_mapping_example.ipynb an example is shown for  matching and the computation of averages.

sprime_mapping_example_filtering.ipynb uses the same data and includes different filtering options (fraction of SNPs which has to be present
in a segment, min nr. of contiguous SNPs within a fragment in a particular segment, extraction of all SNPs).

The necessary input is a SPrime-score file (e.g. sprime.1src.out.100000.score) and the corresponding vcf-file (archie.1_biallelic.vcf.gz).
Furthermore, a file with the target individuals can be provided (archie.1.tgt.ind.list), otherwise for all individuals (incl. reference) the mapping is done.

In principle, the notebook should work for other files by only changing the filenames at the beginning.

The example input is a phased vcf-file with a sequence length of 25MB (50 ref / 10 tgt individuals).
One possibly could improve performance, but for this set-up the mapping takes only 1 sec, so it should be fine.



Here the parameters (can also be found in the function header):

    Parameters
    ----------
    score_file : str, Filepath to the SPrime output file
    vcf_file: str, filepath to the vcf-file which was used as input for SPrime
    out_file : (optional) str, filepath for writing the bed-file with the introgressed tracts per individual; if None, no output file is written and the dataframe is returned
    target_indivdiuals: list, indicating the individuals in the target panel, i.e. for which introgressed fragments shouls be detected (has to correspond to the ids in the vcf_file)
    segment_fraction: (optional) float, if specified, indicates the minimum fraction of SNPs in the target which correspond to the SNPs in the segment provided by SPrime so that it is marked as introgressed in this individual
    min_snps: (optional) int: within a segment, the minimum number of adjacent SNPs so that this region of a segment is marked as introgressed in the specific individual
    phased: (bool) if phased, the inidividual names in the final bed-file are suffixed by '_1' for the first haplotype, '_2' for the second haplotype; if not phased, tracts for both haplotypes are merged. Since SPrime takes diploid phased data as input, one usually will use Phased=True. Currently only diploid input is supported.
    default: True
    merge_distance: (optional) int, for the final merging of tracts, an int > 0 indicates that also regions with a gap should be merged in case that the distance is smaller than merge_distance in base pairs.
    default: 0
    only_tract_output: (optional) True, in case of false, an extra column with 'haplotype' is generated in the bedfile, otherwise the individual names are suffixed as described in phased and no extra haplotype-column is generated
    default: True
    return_full_records: (optional) False: instead of returning the bed-file, a larger dataframe with the number of snps and the segment fraction for each intrgressed region is returned. Overwrites only_tract_output.
    default: False
    
For the generated bed-file representing the introgressed fragments, the structure is: chromosome/start/end/individual

Pandas, Numpy and Pysam have to be installed (should be no problem via 'pip install') as well as pybedtools (perhaps works via 'pip install', but maybe conda-installation via 'conda install' is necessary, in my experience this package sometimes causes problems, and even with conda it does not work on Windows).
It was tested on Python 3.9.19.