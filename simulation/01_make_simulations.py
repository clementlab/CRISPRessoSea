import os
import shutil
import subprocess

import numpy as np


RANDOM_SEED = 7
OUTPUT_FOLDER = "01_make_simulations_output"
TARGET_EDITING_RATES = [x / 100 for x in [100, 50, 25, 10, 5, 1, 0.1, 0.01, 0.001, 0.0001]]
BACKGROUND_MUTATION_RATES = [0, 0.01, 0.05, 0.1]
TOTAL_REPLICATES = 5
READS_PER_TARGET = 10000
TARGET_LENGTH = 100
TARGET_EDITING_POSITIONS = [50, 51]


def make_sample(
    reference_seq,
    mutation_inds,
    rng,
    mutation_rate=0,
    background_mutation_rate=0,
    n_reads=10000,
):
    reads = []

    real_mutations = rng.random(n_reads) < mutation_rate
    background_mutations = rng.random(n_reads) < background_mutation_rate
    mutated_reads = real_mutations | background_mutations
    mutate_inds = rng.choice(mutation_inds, size=n_reads)
    mutate_half_lens = rng.choice([1, 2, 3], size=n_reads)

    for i in range(n_reads):
        if mutated_reads[i]:
            mutate_ind = int(mutate_inds[i])
            mutate_half_len = int(mutate_half_lens[i])
            this_read_seq = (
                reference_seq[: mutate_ind - mutate_half_len]
                + reference_seq[mutate_ind + mutate_half_len :]
            )
        else:
            this_read_seq = reference_seq

        reads.append(
            f"@read_{i}_real{int(real_mutations[i])}_background{int(background_mutations[i])}\n"
            f"{this_read_seq}\n+\n{'F' * len(this_read_seq)}\n"
        )
    return reads


def write_reference_files(output_folder, target_references, target_guides):
    reference_fasta = os.path.join(output_folder, "reference.fasta")
    with open(reference_fasta, "w") as f:
        for idx, reference_seq in enumerate(target_references):
            f.write(f">chr{idx}\n{reference_seq}\n")

    bowtie2_build = shutil.which("bowtie2-build")
    if bowtie2_build is not None:
        subprocess.run(
            [bowtie2_build, reference_fasta, os.path.join(output_folder, "reference")],
            check=True,
        )
    else:
        print("Warning: bowtie2-build not found on PATH. Skipping Bowtie2 index generation.")

    guide_file = os.path.join(output_folder, "guides.txt")
    with open(guide_file, "w") as f:
        f.write("Guide\tTarget\tSequence\tPAM\t#MM\tLocus\n")
        for idx, guide_seq in enumerate(target_guides):
            pam = target_references[idx][
                TARGET_EDITING_POSITIONS[0] + 3 : TARGET_EDITING_POSITIONS[0] + 6
            ]
            f.write(
                f"guide{idx}\tGuide{idx}_EditingRate{TARGET_EDITING_RATES[idx]}"
                f"\t{guide_seq}\t{pam}\t0\tchr{idx}:{TARGET_EDITING_POSITIONS[0]}\n"
            )


def write_sample_files(output_folder, sample_records):
    sample_file = os.path.join(output_folder, "samples.txt")
    with open(sample_file, "w") as f:
        f.write("Name\tr1\tGroup\tbackground_noise_rate\treplicate\n")
        for record in sample_records:
            f.write(
                f"{record['sample_name']}\t{record['fastq_path']}\t{record['group']}"
                f"\t{record['background_noise_rate']}\t{record['replicate']}\n"
            )


def main():
    rng = np.random.default_rng(RANDOM_SEED)

    os.makedirs(OUTPUT_FOLDER, exist_ok=True)

    target_references = []
    target_guides = []
    for _ in TARGET_EDITING_RATES:
        this_reference = "".join(rng.choice(["A", "T", "C", "G"], size=TARGET_LENGTH))
        target_references.append(this_reference)
        target_guides.append(
            this_reference[
                TARGET_EDITING_POSITIONS[0] - 17 : TARGET_EDITING_POSITIONS[1] + 3
            ]
        )

    write_reference_files(OUTPUT_FOLDER, target_references, target_guides)

    sample_records = []
    sample_count = 0
    for replicate in range(TOTAL_REPLICATES):
        for background_mutation_rate in BACKGROUND_MUTATION_RATES:
            for group_name, mutation_rate in [("edited", "target"), ("unedited", 0)]:
                sample_name = (
                    f"sample_{sample_count}_{group_name}_{background_mutation_rate}"
                    f"_noise_rep{replicate}"
                )
                fastq_path = os.path.join(OUTPUT_FOLDER, f"{sample_name}.fq")

                with open(fastq_path, "w") as f:
                    for idx, editing_rate in enumerate(TARGET_EDITING_RATES):
                        this_mutation_rate = editing_rate if mutation_rate == "target" else 0
                        reads = make_sample(
                            target_references[idx],
                            mutation_inds=TARGET_EDITING_POSITIONS,
                            rng=rng,
                            mutation_rate=this_mutation_rate,
                            background_mutation_rate=background_mutation_rate,
                            n_reads=READS_PER_TARGET,
                        )
                        f.writelines(reads)

                sample_records.append(
                    {
                        "sample_name": sample_name,
                        "fastq_path": fastq_path,
                        "group": group_name,
                        "background_noise_rate": background_mutation_rate,
                        "replicate": replicate,
                    }
                )
                sample_count += 1

    write_sample_files(OUTPUT_FOLDER, sample_records)


if __name__ == "__main__":
    main()
