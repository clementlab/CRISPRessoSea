import numpy as np
import os
import subprocess

def make_sample(reference_seq, mutation_inds, mutation_rate = 0, background_mutation_rate = 0, n_reads=10000):
    reads = []
    
    # Generate all boolean flags at once for efficiency 
    real_mutations = np.random.random(n_reads) < mutation_rate
    background_mutations = np.random.random(n_reads) < background_mutation_rate
    
    for i in range(n_reads):
        has_real_mutation = real_mutations[i]
        has_background_mutation = background_mutations[i]
        
        if has_real_mutation or has_background_mutation:
            mutate_ind = np.random.choice(mutation_inds)
            mutate_half_len = np.random.choice([1, 2, 3])
            this_read_seq = reference_seq[:mutate_ind - mutate_half_len] + reference_seq[mutate_ind + mutate_half_len:]
        else:
            this_read_seq = reference_seq[:]
        reads.append(f"@read_{i}_real{int(has_real_mutation)}_background{int(has_background_mutation)}\n{this_read_seq}\n+\n{'F'*len(this_read_seq)}\n")
    return reads


output_folder = "01_make_simulations_output"
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

target_editing_rates = [x/100 for x in [100, 50, 25, 10, 5, 1, 0.1, 0.01, 0.001, 0.0001]]
target_references = []
target_guides = []
target_length = 100
target_editing_positions = [50,51]
for i in range(len(target_editing_rates)):
    this_reference = "".join(np.random.choice(["A", "T", "C", "G"], size=target_length))
    target_references.append(this_reference)
    target_guides.append(this_reference[target_editing_positions[0]-17:target_editing_positions[1]+3])

#generate genome reference file
reference_fasta = os.path.join(output_folder, "reference.fasta")
with open(reference_fasta, "w") as f:
    for idx, reference_seq in enumerate(target_references):
        f.write(f">chr{idx}\n{reference_seq}\n")

#index with bowtie2
subprocess.run(["bowtie2-build", reference_fasta, os.path.join(output_folder, "reference")], check=True)


#generate guide reference file
guide_file = os.path.join(output_folder, "guides.txt")
with open(guide_file, "w") as f:
    f.write("Guide\tTarget\tSequence\tPAM\t#MM\tLocus\n")
    for idx, guide_seq in enumerate(target_guides):
        f.write(f"guide{idx}\tGuide{idx}_EditingRate{target_editing_rates[idx]}\t{guide_seq}\t{target_references[idx][target_editing_positions[0]+3:target_editing_positions[0]+6]}\t0\tchr{idx}:{target_editing_positions[0]}\n")

background_mutation_rates = [0, 0.01, 0.05, 0.1]
sample_count = 0
sample_fastqs = []
sample_names = []
for replicate in range(3):
    for background_mutation_rate in background_mutation_rates:
        # generate edited sample
        edited_fastq_name = os.path.join(output_folder, f"sample_{sample_count}_edited_{background_mutation_rate}_noise_rep{replicate}.fq")
        sample_fastqs.append(edited_fastq_name)
        sample_names.append(f"sample_{sample_count}_edited_{background_mutation_rate}_noise_rep{replicate}")
        with open(edited_fastq_name, "w") as f:
            for idx, editing_rate in enumerate(target_editing_rates):
                reads = make_sample(target_references[idx], mutation_inds=target_editing_positions, mutation_rate=editing_rate, background_mutation_rate=background_mutation_rate, n_reads=10000)
                f.writelines(reads)
        sample_count += 1

        #generate unedited sample
        unedited_fastq_name = os.path.join(output_folder, f"sample_{sample_count}_unedited_{background_mutation_rate}_noise_rep{replicate}.fq")
        sample_fastqs.append(unedited_fastq_name)
        sample_names.append(f"sample_{sample_count}_unedited_{background_mutation_rate}_noise_rep{replicate}")
        with open(unedited_fastq_name, "w") as f:
            for idx, editing_rate in enumerate(target_editing_rates):
                reads = make_sample(target_references[idx], mutation_inds=target_editing_positions, mutation_rate=0, background_mutation_rate=background_mutation_rate, n_reads=10000)
                f.writelines(reads)
        sample_count += 1

#generate sample file
sample_file = os.path.join(output_folder, "samples.txt")
with open(sample_file, "w") as f:
    f.write("Name\tr1\n")
    for idx, sample_fastq in enumerate(sample_fastqs):
        f.write(f"{sample_names[idx]}\t{sample_fastq}\n")


        

        




            
            
