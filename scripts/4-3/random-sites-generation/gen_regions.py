# this program generates m random regions across the chromosomes
# first it picks a value between 0 and the sum of the lenght of the chromosomes, in bp.
# then it locates this position in the genome
# then it extends this position to create a region
# the extension length is taken from a random sample of lenghts
#   taken from the R script linked
# the distribution of lenghts is modelled after the distribution of the lenghts of
# the MACS2 regions

import random

# libraries for interfacing with R
import pandas as pd
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri

# setting up the interface with R
r = robjects.r
r['source']('sample_len.r')
get_sample_r = robjects.globalenv['get_sample']

# lengths of the chromosomes
chr_len = { '1' : 249250621,
            '2' : 243199373,
            '3' : 198022430,
            '4' : 191154276,
            '5' : 180915260,
            '6' : 171115067,
            '7' : 159138663,
            '8' : 146364022,
            '9' : 141213431,
            '10' : 135534747,
            '11' : 135006516,
            '12' : 133851895,
            '13' : 115169878,
            '14' : 107349540,
            '15' : 102531392,
            '16' : 90354753,
            '17' : 81195210,
            '18' : 78077248,
            '19' : 59128983,
            '20' : 63025520,
            '21' : 48129895,
            '22' : 51304566,
            'X' : 155270560,
            'Y' : 59373566 }
            
# link between chromosome name and index
index_to_chr_name = { 0 : '1',
                      1 : '2',
                      2 : '3',
                      3 : '4',
                      4 : '5',
                      5 : '6',
                      6 : '7',
                      7 : '8',
                      8 : '9',
                      9 : '10',
                      10 : '11',
                      11 : '12',
                      12 : '13',
                      13 : '14',
                      14 : '15',
                      15 : '16',
                      16 : '17',
                      17 : '18',
                      18 : '19',
                      19 : '20',
                      20 : '21',
                      21 : '22',
                      22 : 'X',
                      23 : 'Y' }
            
            
def main(m):
    # length of the genome, in bp
    total_genome_len = sum(chr_len.values())
    
    # random sample of m position across the genome, uniformly distributed
    random_spots = random.sample(range(0, total_genome_len), m)

    # converts the m random positions in 'random_spots' in m pairs 
    # (chr index, pos in chr)
    pairs = [absolute_to_chr(x, list(chr_len.values())) for x in random_spots]
    
    # converts from chromosome index to chromosome name in the pairs
    pairs = [(index_to_chr_name.get(x[0]), x[1]) for x in pairs]
    
    # get a sample of m lengths, these will be the lenghts of the produced
    # random regions: the random spot will be extended upstream and downstream
    # to match this length
    sample_lengths = get_sample_r(m)
    
    # extending each spot (x in pairs) to have the length in sample_lenghts
    extended = [extend(x, round(y/2)) for x,y in zip(pairs, sample_lengths)]
    
    filename = "output_random_regions.bed"
    write_bed(filename, extended)
    print(filename + ": " + str(m) + " regions written")
    
# given an absolute position, it returns the
# chrosomes in which it is located and the position
# on that chromosome
def absolute_to_chr(abspos, chrs):
    cumulative_chrs = build_cumulative_chrs(chrs)
    chr_index = find_pos(abspos, cumulative_chrs)
    rel_pos = get_relative_pos(cumulative_chrs, chr_index, abspos)
    return (chr_index, rel_pos)
    
# given the list with the length of the
# chromosomes, build the cumulative version:
# the start of the nth is the end of the (n-1)th+1
def build_cumulative_chrs(chrs):
    cumulative = chrs.copy()
    for i in range(1, len(cumulative)):
        cumulative[i] += cumulative[i-1]
    return cumulative
                  
 
# given pos absolute position and cumulative_chrs
# list with the ending position of the chromosomes
# returns the index of the chromosome that contains
# the position pos
def find_pos(pos, cumulative_chrs):
    i = 0
    while (pos > cumulative_chrs[i]):
        i+=1
    return i

# given the array with the end position of the chromosomes,
# the absolute position of a spot and the index of the chromosome
# in which it is located, returns the relative position of the 
# spot on the chromosome it is located
def get_relative_pos(cumulative_chrs, chr_index, abspos):
    if chr_index == 0:
        return abspos
    # now index > 0
    start_pos = cumulative_chrs[chr_index-1]+1
    return(abspos - start_pos)
    
# given a tuple (chr, pos), extendes pos by ext in both directions
# returns a tuple (chr, (startpos, endpos)) = (chr, (pos-ext, pos+ext))
def extend(tup, ext):
    pos = tup[1]
    start = pos-ext
    if start < 0:
        start = 0
    ret = (tup[0], (start, pos+ext))
    return ret
    
    
# given the list of tuples, writes the corresponding bed file
def write_bed(filename, tuples):
    file = open(filename, "w")
    for i in range(0, len(tuples)):
        current_chr = (tuples[i])[0]
        current_start = ((tuples[i])[1])[0]
        current_end = ((tuples[i])[1])[1]
        line = "chr" + current_chr + "\t" + str(current_start) + "\t" + str(current_end) + "\t" + "random_peak_" + str(i) + "\n"
        file.write(line)
    file.close()

    
 
if __name__ == "__main__" :
    m = 11028
    main(m)
