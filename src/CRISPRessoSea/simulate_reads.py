import argparse
import gzip
import os
import pandas as pd
import subprocess
import random

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N',
                  'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'}
    return "".join(complement.get(base, base) for base in reversed(seq))

def get_genome_seq(genome_fasta, chrom, start, end):
    # samtools faidx uses 1-based inclusive coordinates
    cmd = f"samtools faidx {genome_fasta} {chrom}:{start}-{end}"
    try:
        output = subprocess.check_output(cmd, shell=True).decode('utf-8')
        lines = output.strip().split('\n')
        if len(lines) > 1:
            return "".join(lines[1:])
        return ""
    except Exception as e:
        print(f"Error fetching sequence for {chrom}:{start}-{end}: {e}")
        return ""

def get_random_amplicons(genome_fasta, n, window):
    fai_path = genome_fasta + ".fai"
    if not os.path.exists(fai_path):
        print(f"Index {fai_path} not found. Cannot generate random reads.")
        return []
    
    chroms = []
    with open(fai_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            chroms.append((parts[0], int(parts[1])))
    
    amplicons = []
    while len(amplicons) < n:
        chrom, length = random.choice(chroms)
        if length <= window:
            continue
        start = random.randint(1, length - window)
        end = start + window - 1
        seq = get_genome_seq(genome_fasta, chrom, start, end)
        if seq and len(seq) == window:
            seq_upper = seq.upper()
            if seq_upper.count('N') / len(seq_upper) <= 0.05:
                amplicons.append(seq_upper)
    return amplicons

def main():
    parser = argparse.ArgumentParser(description="Generate simulated reads for amplicons.")
    parser.add_argument("-s", "--sample", required=True, help="Sample file (TSV with Name column)")
    parser.add_argument("-t", "--target", required=True, help="Target file (TSV with target and Locus columns)")
    parser.add_argument("-x", "--genome", required=True, help="Reference genome FASTA file")
    parser.add_argument("-a", "--alternate_alleles", default=None, help="Alternate alleles file")
    parser.add_argument("-n", "--n_sequences", type=int, default=100, help="Number of sequences per amplicon per sample")
    parser.add_argument("-w", "--window", type=int, default=200, help="Window size for genomic amplicons")
    parser.add_argument("--n_random_reads", type=int, default=100, help="Number of random off-target reads to add from the genome")
    parser.add_argument("--paired", action="store_true", help="Generate paired-end reads")
    parser.add_argument("--overlap", type=int, default=50, help="Overlap length for paired-end reads (default: 50)")
    parser.add_argument("-o", "--output_folder", default="simulated_reads", help="Output directory")

    args = parser.parse_args()

    os.makedirs(args.output_folder, exist_ok=True)

    # 1. Parse targets
    target_df = pd.read_csv(args.target, sep='\t')
    # find target and locus columns (case insensitive ideally, but let's assume 'target' and 'Locus')
    target_col = [c for c in target_df.columns if c.lower() == 'target'][0]
    locus_col = [c for c in target_df.columns if c.lower() == 'locus'][0]

    # 2. Parse alternate alleles
    alt_alleles = {}
    if args.alternate_alleles and os.path.exists(args.alternate_alleles):
        alt_df = pd.read_csv(args.alternate_alleles, sep='\t')
        if 'target_name' in alt_df.columns and 'reference_seqs' in alt_df.columns:
            for _, row in alt_df.iterrows():
                alt_alleles[row['target_name']] = row['reference_seqs']
        elif 'region_name' in alt_df.columns and 'reference_seqs' in alt_df.columns:
             for _, row in alt_df.iterrows():
                alt_alleles[row['region_name']] = row['reference_seqs']

    # 3. Generate amplicon sequences
    amplicons = {}
    for _, row in target_df.iterrows():
        t_name = row[target_col]
        locus = row[locus_col]
        if pd.isna(locus) and t_name not in alt_alleles:
            continue
        
        if t_name in alt_alleles:
            seq = alt_alleles[t_name]
        else:
            # Parse locus: chr10:-14934460 or chr10:14934460
            locus = str(locus).replace('+', '')
            parts = locus.split(':')
            chrom = parts[0]
            pos_str = parts[1]
            # handle negative sign which might indicate strand
            pos = int(pos_str.replace('-', ''))
            start = pos - args.window // 2
            end = pos + args.window // 2 - 1
            seq = get_genome_seq(args.genome, chrom, start, end)
        
        if seq:
            amplicons[t_name] = seq.upper()
    
    print(f"Generated {len(amplicons)} amplicon sequences.")
    
    random_amplicons = []
    if args.n_random_reads > 0:
        print(f"Generating {args.n_random_reads} random background reads...")
        random_amplicons = get_random_amplicons(args.genome, args.n_random_reads, args.window)

    # 4. Parse samples and write fastqs
    sample_df = pd.read_csv(args.sample, sep='\t')
    sample_col = [c for c in sample_df.columns if c.lower() == 'name'][0]
    
    new_sample_rows = []

    for _, row in sample_df.iterrows():
        sample_name = row[sample_col]
        
        r1_path = os.path.join(args.output_folder, f"{sample_name}_R1.fastq.gz")
        r2_path = os.path.join(args.output_folder, f"{sample_name}_R2.fastq.gz")
        
        new_row = {'Name': sample_name, 'fastq_r1': r1_path}
        
        with gzip.open(r1_path, 'wt') as f1:
            if args.paired:
                new_row['fastq_r2'] = r2_path
                with gzip.open(r2_path, 'wt') as f2:
                    read_idx = 1
                    for t_name, seq in amplicons.items():
                        L = len(seq)
                        R = min((L + args.overlap) // 2, L)
                        seq_r1 = seq[:R]
                        seq_r2 = reverse_complement(seq[-R:])
                        qual_r1 = "I" * len(seq_r1)
                        qual_r2 = "I" * len(seq_r2)
                        
                        for _ in range(args.n_sequences):
                            header1 = f"@{sample_name}_{read_idx} 1:N:0:1"
                            header2 = f"@{sample_name}_{read_idx} 2:N:0:1"
                            
                            f1.write(f"{header1}\n{seq_r1}\n+\n{qual_r1}\n")
                            f2.write(f"{header2}\n{seq_r2}\n+\n{qual_r2}\n")
                            read_idx += 1
                            
                    for seq in random_amplicons:
                        L = len(seq)
                        R = min((L + args.overlap) // 2, L)
                        seq_r1 = seq[:R]
                        seq_r2 = reverse_complement(seq[-R:])
                        qual_r1 = "I" * len(seq_r1)
                        qual_r2 = "I" * len(seq_r2)
                        
                        header1 = f"@{sample_name}_{read_idx} 1:N:0:1"
                        header2 = f"@{sample_name}_{read_idx} 2:N:0:1"
                        
                        f1.write(f"{header1}\n{seq_r1}\n+\n{qual_r1}\n")
                        f2.write(f"{header2}\n{seq_r2}\n+\n{qual_r2}\n")
                        read_idx += 1
            else:
                read_idx = 1
                for t_name, seq in amplicons.items():
                    qual = "I" * len(seq)
                    for _ in range(args.n_sequences):
                        header = f"@{sample_name}_{read_idx} 1:N:0:1"
                        f1.write(f"{header}\n{seq}\n+\n{qual}\n")
                        read_idx += 1
                        
                for seq in random_amplicons:
                    qual = "I" * len(seq)
                    header = f"@{sample_name}_{read_idx} 1:N:0:1"
                    f1.write(f"{header}\n{seq}\n+\n{qual}\n")
                    read_idx += 1
        
        new_sample_rows.append(new_row)
        print(f"Generated {read_idx - 1} reads for sample {sample_name}.")
        
    # Write a new sample file for convenience
    new_sample_df = pd.DataFrame(new_sample_rows)
    new_sample_file = os.path.join(args.output_folder, "simulated_samples.txt")
    new_sample_df.to_csv(new_sample_file, sep='\t', index=False)
    print(f"Done. Wrote new sample file to {new_sample_file}")
    
    cmd = f"CRISPRessoSea Process -t {args.target} -s {new_sample_file} -x {args.genome} -o Output/"
    if args.alternate_alleles:
        cmd += f" --alternate_alleles {args.alternate_alleles}"
    
    print("\n" + "="*80)
    print("You can run CRISPRessoSea with the simulated reads using this command:")
    print(cmd)
    print("="*80 + "\n")

if __name__ == '__main__':
    main()
