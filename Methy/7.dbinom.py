import csv
from scipy.stats import binom

def process_line(line):
    fields = line.strip().split("\t")
    n = fields[0]
    p = fields[1]
    d = fields[2]
    sc = int(fields[3])
    fc = int(fields[4])
    context = fields[5]
    re = fields[6]
    return n, p, d, sc, fc, context, re

def calculate_p_value(sc, fc):
    count = sc + fc
    p_value = binom.pmf(sc, count, 0.005)
    return p_value

def process_file(input_file_path, output_file_path):
    with open(input_file_path, 'r') as infile, open(output_file_path, 'w') as outfile:
        reader = csv.reader(infile, delimiter='\t')
        writer = csv.writer(outfile, delimiter='\t')
        
        mydata = []
        p_values = []

        for line in reader:
            n, p, d, sc, fc, context, re = process_line('\t'.join(line))
            if sc >= 3:
                p_value = calculate_p_value(sc, fc)
                if p_value <= 1e-4:
                    row = [n, p, d, sc, fc, context, re, p_value]
                    writer.writerow(row)

if __name__ == "__main__":
    input_file_path = "/public/home/phwan/Gossypium-barbadense/ViewBS1/Vpre/5H6-7-1/5H6-7-1_View1.tab"
    output_file_path = "/public/home/phwan/Gossypium-barbadense/Extractor2/5H6-7-1/5H6-7-1_dbionom2.CX.txt.CX_report.txt"

    process_file(input_file_path, output_file_path)

