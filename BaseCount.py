from collections import Counter
import Bio.SeqIO
import gzip
from multiprocessing import Pool
import os
import argparse
from typing import List, Union, Dict


def check_input(fastq_list: List[Union[str, bytes, os.PathLike]]):
    for i in fastq_list:
        try:
            os.path.isfile(i)
        except Exception as e:
            print(f"Exception occurred for file {i}: {e}")
            exit(1)

def create_output(output_path: Union[str, bytes, os.PathLike]):
    if os.path.isdir(output_path):
        print(f"Warning: output folder '{output_path}'already exists.")
    else:
        try :
            os.makedirs(output_path)
        except Exception as e:
            print(f"Exception creating output folder: {e}")
            exit(1)

def write_output(fastq_file_path: Union[str, bytes, os.PathLike], output_path: Union[str, bytes, os.PathLike],
                 counts: List):
    output_count_file = open(output_path + "/" + fastq_file_path.split("/")[-1].rstrip(".gz").rstrip(".fastq") +
                             "_base_counts.tsv", "w")
    output_count_file.write("count_A" + "\t" + "count_T" + "\t" + "count_C" + "\t" + "count_G" + "\n")
    for i in counts:
        output_count_file.write(f"{i['A']}\t{i['T']}\t{i['C']}\t{i['G']}\n")
    print(f"Finished writing results for {fastq_file_path}", flush=True)
    output_count_file.close()


def compute_frequencies(sequence: str) -> Dict[str, float]:
    freq = dict()
    counts = Counter(sequence)
    freq["A"] = counts['A'] / len(sequence)
    freq["C"] = counts['C'] / len(sequence)
    freq["G"] = counts['G'] / len(sequence)
    freq["T"] = counts['T'] / len(sequence)
    return freq


def count(fastq_file_path: Union[str, bytes, os.PathLike], output_path: Union[str, bytes, os.PathLike]):
    print(f"Counting reads in {fastq_file_path}", flush=True)
    results = list()
    if fastq_file_path.endswith(".gz"):
        with gzip.open(fastq_file_path, "rt") as handle:
            for record in Bio.SeqIO.parse(handle, 'fastq'):
                results.append(compute_frequencies(record))
    else:
        with open(fastq_file_path, "r") as handle:
            for record in Bio.SeqIO.parse(handle, 'fastq'):
                results.append(compute_frequencies(record))
    write_output(fastq_file_path, output_path, counts=results)


def main():
    parser = argparse.ArgumentParser(description="Counts base proportions of each sequences of multiple fasta files in parallel.",
                                     epilog="Ex: python3 BaseCount.py --fasta $(find ./test/*.fastq) --threads 10 --output ./output")
    parser.add_argument("-f", "--fastq", help="Fastq file(s), may be gziped, space-sep.", type=str, required=True, nargs="+")
    parser.add_argument("-t", "--threads", help="Number of threads", type=int, default=1, required=False)
    parser.add_argument("-o", "--output", help="Output directory", type=str, required=True, action='store')

    args = parser.parse_args()

    create_output(args.output)
    with Pool(processes=args.threads) as pool:
        for i in args.fastq:
            pool.apply_async(count, args=(i, args.output))
        pool.close()
        pool.join()
    print(f"Done", flush=True)



if "__main__" == __name__:
    main()