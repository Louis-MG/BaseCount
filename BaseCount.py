from collections import Counter
from Bio import SeqIO
import gzip
from multiprocessing import Pool, Value
import os
import argparse
from typing import List, Union


def create_output(output_path: Union[str, bytes, os.PathLike]):
    if os.path.isdir(output_path):
        print('Warning: output folder already exists')
    else:
        os.makedirs(output_path)


def write_output(fastq_file_path: Union[str, bytes, os.PathLike], output_path: Union[str, bytes, os.PathLike],
                 counts: List):
    annotated_reads_file = open(output_path + "/" + fastq_file_path.split("/")[-1].rstrip(".fastq") + "_base_counts.tsv",
                                "w")
    annotated_reads_file.write("count_A" + "\t" + "count_T" + "\t" + "count_C" + "\t" + "count_G" + "\n")
    for i in counts:
        annotated_reads_file.write(i["A"] + "\t" + i["T"] + "\t" + i["C"] + "\t" + i["G"] + "\n")
    print(f"Finished writing results for {fastq_file_path}", flush=True)
    annotated_reads_file.close()


def count(fasta_file_path: Union[str, bytes, os.PathLike], output_path: Union[str, bytes, os.PathLike]):
    results = list()
    with gzip.open(fasta_file_path, "rt") as handle:
        for record in SeqIO.parse(handle, 'fastq'):
            freq = dict()
            counts = Counter(record.seq)
            freq["A"] = counts['A']/len(record.seq)
            freq["C"] = counts['C']/len(record.seq)
            freq["G"] = counts['G']/len(record.seq)
            freq["T"] = counts['T']/len(record.seq)
            results.append(freq)
    write_output(fasta_file_path, output_path, counts=results)


def main():
    parser = argparse.ArgumentParser(description="Counts base proportions of each sequences of multiple fasta files in parallel.",
                                     epilog="Ex: python3 BaseCount.py --fasta $(find ./test/*.fastq) --threads 10 --output ./output")
    parser.add_argument("-f", "--fasta", help="Fasta file(s)", type=str, required=True, action="append")
    parser.add_argument("-t", "--threads", help="Number of threads", type=int, default=1, required=False)
    parser.add_argument("-o", "--output", help="Output directory", type=str, required=True, action='store')

    args = parser.parse_args()

    create_output(args.output)
    with Pool(processes=args.threads) as pool:
        pool.apply_async(count, args=(args.fasta, args.output))
        pool.close()
        pool.join()
    print(f"Done", flush=True)



if "__main__" == __name__:
    main()