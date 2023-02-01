import pysam
from pyfaidx import Fasta
import re
import fire
import os
import numpy as np

def extract_clean_sequence(ibam="input.bam",reference="hg38.fa",outpath="raw_reads",statistic="statistic.txt",metric="metrics.txt",leftS=0,leftE=113):
    ref_reader = Fasta(reference)
    bam_file = pysam.AlignmentFile(ibam, 'rb')
    library_count = {}
    for record in ref_reader:
        raw_count = 0
        clean_count = 0
        indel_count = 0
        mismatch_count = 0
        out_file = os.path.join(outpath, record.name+".out")
        with open(out_file,"w+") as f:
            ref_sequence = str(record)
            for read in bam_file.fetch(record.name):
                raw_count += 1
                left_flank = ref_sequence[leftS:leftE].upper()
                pattern = re.compile(r'{}[AGCT]+GTTT'.format(left_flank))
                search_result = pattern.findall(read.query_sequence)
                if len(search_result) > 0:
                    clean_count += 1
                    out_sequence = search_result[0][leftE:-4]
                    control_sequence = ref_sequence[leftE:]
                    if len(out_sequence)!=len(control_sequence):
                        indel_count += 1
                    else:
                        if out_sequence != control_sequence:
                            mismatch_count += 1
                    f.write(out_sequence+"\n")
        library_count[record.name] = [raw_count,clean_count,mismatch_count,indel_count,clean_count/raw_count*100]
    ## write statistic file
    with open(statistic,"w+") as s:
        s.write("library\ttotal\tclean\tmismatch\tindel\tdepth\n")
        for (key,value) in library_count.items():
            s.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(key,value[0],value[1],value[2],value[3],value[4]))
    ## write metrics file
    with open(metric,"w+") as m:
        hit_numbers = len([key for (key,value) in library_count.items() if value[1] > 0])
        hit_ratio = hit_numbers/len(library_count)*100
        average_depth = np.mean([value[4] for (_,value) in library_count.items()])
        m.write("hit_numbers\t{}\n".format(hit_numbers))
        m.write("hit_ratio\t{}%\n".format(hit_ratio))
        m.write("average_depth\t{}%\n".format(average_depth))

if __name__ == '__main__':
    fire.Fire(extract_clean_sequence)



