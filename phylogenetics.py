import os
import subprocess
import glob
import argparse

def run_command(command):
    """Runs a command and exits if it fails."""
    print(f"Executing: {command}")
    try:
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error executing command: {command}\n{e}")
        exit(1)

def main():
    """Main function to run the entire pipeline."""
    parser = argparse.ArgumentParser(description="Automated viral genomics and phylogenetics pipeline.")
    parser.add_argument(
        '--ref_id',
        type=str,
        default='NC_045512.2',
        help='The NCBI accession number for the reference genome.'
    )
    parser.add_argument(
        'sample_ids',
        type=str,
        nargs='+',
        help='One or more SRA Run IDs (e.g., SRR34928089) to analyze.'
    )
    args = parser.parse_args()

    REF_GENOME_ID = args.ref_id
    sample_ids = args.sample_ids
    REF_GENOME_FILE = "reference.fa"

    if not os.path.exists(REF_GENOME_FILE):
        print("--- Downloading Reference Genome ---")
        command = f"efetch -db nuccore -id {REF_GENOME_ID} -format fasta > {REF_GENOME_FILE}"
        run_command(command)
    else:
        print("--- Reference Genome already exists. Skipping download. ---")

    os.makedirs("all_genomes", exist_ok=True)

    print("\n--- Starting Sample Processing ---")
    for sample_id in sample_ids:
        print(f"\n--- Processing Sample: {sample_id} ---")

        fastq1 = f"{sample_id}_1.fastq.gz"
        fastq2 = f"{sample_id}_2.fastq.gz"
        snippy_dir = f"{sample_id}_snippy"

        if not os.path.exists(snippy_dir):
            command = f"fasterq-dump --split-files {sample_id}"
            run_command(command)
            command = (f"snippy --force --outdir {snippy_dir} "
                       f"--ref {REF_GENOME_FILE} "
                       f"--R1 {fastq1} "
                       f"--R2 {fastq2}")
            run_command(command)
            if os.path.exists(fastq1): os.remove(fastq1)
            if os.path.exists(fastq2): os.remove(fastq2)
        else:
            print(f"Snippy directory '{snippy_dir}' already exists. Skipping analysis.")

        consensus_file = os.path.join(snippy_dir, "snps.consensus.fa")
        destination_file = os.path.join("all_genomes", f"{sample_id}.fa")
        if os.path.exists(consensus_file):
            with open(consensus_file, 'r') as original, open(destination_file, 'w') as destination:
                header = original.readline()
                destination.write(f">{sample_id}\n")
                for line in original:
                    destination.write(line)
        else:
            print(f"Warning: Consensus file not found for {sample_id}")

    print("\n--- Finalizing Phylogeny ---")
    all_consensus_genomes_file = "all_consensus_genomes.fa"
    genome_files = glob.glob("all_genomes/*.fa")

    if len(genome_files) < 2:
        print("\nWarning: Need at least 2 genomes to build a tree. Skipping tree building.")
    else:
        with open(all_consensus_genomes_file, 'w') as outfile:
            for fname in genome_files:
                with open(fname) as infile:
                    outfile.write(infile.read())

        command = f"iqtree -s {all_consensus_genomes_file} -m GTR+G -nt AUTO"
        run_command(command)
        print(f"\nYour final evolutionary tree is in the file: {all_consensus_genomes_file}.treefile")

    print("\n🚀 Pipeline Finished Successfully!")

if __name__ == "__main__":
    main()

