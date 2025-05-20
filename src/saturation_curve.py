import click
import numpy as np

from typing import List, Tuple

"""
Name: saturation_curve
Author: Eugenio Mattei
Affiliation: BroadInstitute of MIT and Harvard
Version: 1.0.0
Description: This script provides functionality to generate a saturation curve from scATAC fragment file.
"""


def sample_counts(
        duplicates: int,
        downsample_probabilities: List = np.arange(0.1, 1, 0.1)
        ) -> Tuple:
    """
    Calculate the number of duplicates and total counts for each downsample probability.

    Parameters:
    - duplicates: Number of duplicates to sample.
    - downsample_rate: The rate at which to downsample the fragments.

    Returns:
    - Tuple of arrays containing total counts and duplicate counts.
    """
    # Calculate the total counts and duplicate counts based on the downsample probabilities.
    # The total counts are sampled from a binomial distribution.
    total_counts = np.random.binomial(duplicates, downsample_probabilities)
    # The duplicate counts are calculated as the total counts minus 1.
    # Ensure the duplicate counts are non-negative.
    duplicate_counts = np.maximum(total_counts - 1, 0)

    return total_counts, duplicate_counts


def saturation_from_list(
        duplicates: List,
        downsample_probabilities: List = np.arange(0.1, 1, 0.1)
        ) -> Tuple:
    """
    Downsample the fragments in the input file and write to the output file.

    Parameters:
    - duplicates:
        List containing the number of duplicates to sample per fragment.
        This is the 5th column of a fragment file.
    - downsample_rate:
        The rate at which to downsample the fragments.
        This should be a range or list of probabilities between 0 and 1.

    Returns:
    - Tuple of arrays containing total counts and duplicate counts.
    """
    # Placeholder for downsampling logic
    print(f"Downsampling with rates {downsample_probabilities}")
    # Initialize the sampled counts and duplicates
    sampled_counts = np.zeros(len(downsample_probabilities))
    sampled_duplicates = np.zeros(len(downsample_probabilities))
    # Iterate over the downsample rates
    for dup in duplicates:
        # Sample counts and duplicates for the given downsample rate
        temp_counts, temp_duplicates = sample_counts(dup, downsample_probabilities)
        # Accumulate the results
        sampled_counts += temp_counts
        sampled_duplicates += temp_duplicates

    return sampled_counts, sampled_duplicates


@click.command()
@click.option('--input', type=click.File('r'), help="Input file name or '-' to read from stdin.", required=True)
@click.option('--downsample_rate', type=float, default=0.1, help="Downsample rate.")
@click.option('--output_prefix', type=str, default='saturation_curve', help="Output prefix.")
@click.option('--plot', is_flag=True, help="Generate a plot of the saturation curve.")
def main(input, downsample_rate, output_prefix, plot):
    """
    This function takes in input a fragment file, takes the 5th column with the duplicate numbers
    and each downsample rate, it calculates the number of duplicates and total counts.

    Args:\n
        input (file): Input file containing the fragment file. Use '-' to read from stdin. (Required).\n
        output_prefix (str): Prefix for the output files (Optional).\n
        downsample_rate (float): The rate at which to downsample the fragments. (Optional).\n
        plot (bool): If True, generates a plot of the saturation curve. (Optional).\n

    Returns:\n
        saturation_curve (File): A file containing the saturation curve data.\n


    Output format:\n
        The output is a tab-separated file called `saturation_curve.txt` with the following columns:\n
        1. downsample probability\n
        2. total counts\n
        3. duplicate counts\n
        4. Percentage of duplicates\n

        If the `-p` flag is set, a plot of the saturatin curve will be generated and saved as `saturation_curve.png`.


    Usage:\n
    python saturation_curve.py -i fragment_file.tsv --plot\n
    or\n
    gzip -dc fragment_file.tsv.gz | python saturation_curve.py --input - --downsample_rate 0.1 -o output_prefix --plot
    """
    rate = downsample_rate
    downsample_probabilities = np.arange(0.1, 1, rate)
    output_size = len(downsample_probabilities)

    # Define the two vectors holding the results
    total_counts = np.zeros(output_size)
    total_duplicates = np.zeros(output_size)
    #genomic_regions = set()
    unique_fragments = 0
    # Iterate over the fragment file
    for line in input:
        content = line.strip().split("\t")
        # Extract the number of duplicates from the 5th column
        dup = int(content[4])
        # Sample counts and duplicates for the given downsample rate
        sampled_counts, sampled_duplicates = sample_counts(dup, downsample_probabilities)
        # Accumulate the results
        total_counts += sampled_counts
        total_duplicates += sampled_duplicates
        unique_fragments += total_counts > 0

    duplicate_rate = (total_duplicates - unique_fragments) / unique_fragments
    effective_library_size = total_counts / (1 + duplicate_rate)

    
    # Save the results to a file
    output_file = f"{output_prefix}_saturation_curve.txt"
    with open(output_file, 'w') as f:
        # Write the header
        f.write("#downsample_probability\ttotal_counts\tduplicate_counts\tpercentage_duplicates\tduplicate_rate\teffective_library_size\n")
        for i,_ in enumerate(downsample_probabilities):
            prob = downsample_probabilities[i]
            total = total_counts[i]
            dup = total_duplicates[i]
            # Calculate the percentage of duplicates
            percentage_duplicates = (dup / total) * 100 if total > 0 else 0
            # Retrieve duplicate rate and effective library size using index
            dup_rate = duplicate_rate[i]
            eff_lib_size = effective_library_size[i]
            f.write(f"{prob:.1f}\t{total}\t{dup}\t{percentage_duplicates:.2f}\t{dup_rate:.2f}\t{eff_lib_size:.2f}\n")


if __name__ == "__main__":
    main()
