from bitarray import bitarray
import gzip
import codecs
import argparse
import sys
from Bio import SeqIO
from scipy.interpolate import interp1d

def read_processed_data(hap_file_name, pos_file_name):
    """
    Reads data.
    """
    data = list()
    with gzip.open(hap_file_name, 'rt') as hap_file:
        with gzip.open(pos_file_name, 'rt') as pos_file:        
            for info, variants in zip(pos_file, hap_file):
                chrom, lab, gen_pos, phys_pos = info.split()
                ba = bitarray(variants.strip())
                data.append((int(phys_pos), int(phys_pos), ba))

    return data


def read_cosi_data(hap_file_name, pos_file_name):
    """
    Reads data.
    """
    positions = list()
    snps = list()
    with open(hap_file_name, 'r') as hap_file:
        for line in hap_file:
            _id, pop, *alleles = line.split()
            snps.append(list(alleles))
    with open(pos_file_name, 'r') as pos_file:        
        header = next(pos_file)
        assert header.startswith('SNP')
        for line in pos_file:
            snp, chrom, phys_pos, allele1, freq1, allele2, freq2 = line.split()
            positions.append(phys_pos)
    data = list()
    for pos, alleles in zip(positions, map(list, zip(*snps))):
        ba = bitarray(''.join(alleles).replace('2', '0'))
        data.append((int(pos), int(pos), ba))

    return data


def read_rec_map(file_name):

    with open(file_name, 'r') as f:

        phys_positions = list()
        map_positions = list()

        for line in f:
            chrom, _, map_pos, phys_pos = line.split()
            map_pos, phys_pos = float(map_pos), int(phys_pos)
            phys_positions.append(phys_pos)
            map_positions.append(map_pos)

    return phys_positions, map_positions


def read_phased_vcf_data(vcf_file_name, ancestral_seq_file_name, recmap_file_name=None, haploid=False):

    ancestral = str(SeqIO.read(ancestral_seq_file_name, "fasta").seq)

    if recmap_file_name:
        phys_pos, map_pos = read_rec_map(recmap_file_name)
        phys2map = interp1d(phys_pos, map_pos)

    if vcf_file_name == 'stdin':
        vcf_file = sys.stdin
    else:
        vcf_file = open(vcf_file_name)

    data = list()
    for line in sys.stdin:
        if line.startswith('#'):
            continue
        chrom, pos, snpid, ref, alt, qual, filt, info, form, *genotypes = line.split()
        pos = int(pos)
        hap_list = list()

        for genotype in genotypes:
            if haploid != (len(genotype) == 1):
                # not in accordance with haploid bool
                continue
            if haploid:
                hap_list.append(genotype)
            else:
                a, _, b = genotype
                hap_list.append(a)
                hap_list.append(b)
        if not hap_list:
            continue

        ba = bitarray(''.join(hap_list))

        ancestral_base = ancestral[pos-1]
        if ancestral_base not in 'ATGC':
            # only uppercase are called reliably
            continue
        if alt == ancestral_base:
            ba.invert()

        if recmap_file_name:
            gen_pos = phys2map(pos)
            data.append((pos, gen_pos, ba))
        else:
            data.append((pos, pos, ba))

    return data


def sq_freq_pairs(lst, pos, end, size, max_ancest_freq, min_tot_freq, gen_dist=False):
    """
    Computes pairs of squared freqencies of derived variants to the
    left or right of the focal SNP. Called twice for each focal SNP.
    """
    focal_pos, focal_gen_pos, focal_deriv = lst[pos]
    focal_ancest = ~focal_deriv
    # focal_deriv &= mask
    # focal_ancest &= mask
    focal_deriv_count = focal_deriv.count()
    focal_ancest_count = focal_ancest.count()
    
    f_vals = list()        
    positions = pos < end and range(pos+1, end) or range(pos-1, end-1, -1)
    for i in positions:
        deriv_pos, deriv_gen_pos, deriv = lst[i]
        # deriv &= mask

        if gen_dist and abs(focal_gen_pos - deriv_gen_pos) > size:
            break
        if not gen_dist and abs(focal_pos - deriv_pos) > size:
            break
                
        f_d = (focal_deriv & deriv).count() / focal_deriv_count
        f_a = (focal_ancest & deriv).count() / focal_ancest_count
#        f_tot = deriv.count() / deriv.length()
        f_tot = deriv.count() / len(deriv)

        f_d2 = 0
        f_a2 = 0

        if f_d > 0.0000001 and f_d < 1:
            f_d2 = 1

        if f_a > 0.0000001 and f_a < 1:
            f_a2 = 1
            
        f_vals.append((f_d2,f_a2))

            #print(f_tot,f_d,f_a)
    return f_vals


epilog = """
Example command line with simulated cosi input:
    python hapdaf.py --pos tmp.pos --hap tmp.hap
Example command line with filtered VCF input piped to stdin: 
    vcftools --gzvcf /home/kmt/simons/faststorage/data/1000Genomes/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
--remove-indels --remove-filtered-all --max-alleles 2 --recode --recode-INFO-all --stdout | python hapdaf.py --vcf stdin \
--ancestral /home/kmt/simons/faststorage/data/1000Genomes/ancestral_alignmens/human_ancestor_GRCh37_e59/human_ancestor_22.fa \
--recmap data/plink_maps/plink.GRCh37.map/plink.chr22.GRCh37.map --window 0.05 --outfile stats_chr22.txt
"""

parser = argparse.ArgumentParser(description='Computes hapDAF statistic by David Enard.', 
                                 epilog=epilog, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-o", "--outfile", dest="out_file",
                    type=str,
                    help="Name of output file.")
parser.add_argument("-w", "--windowsize",
                    dest="window_size",
                    type=float,
                    default=50000,
                    help="Size of window. In bases or cM/Mb as appropriate.")
parser.add_argument("-a", "--max-ancest-freq",
                    dest="max_ancest_freq",
                    type=int,
                    default=1,
                    help="Maximum frequency of derived flanking variants on background.")
parser.add_argument("-t", "--min-tot-freq",
                    dest="min_tot_freq",
                    type=int,
                    default=0,
                    help="Minimum frequency of derived flanking variants.")
parser.add_argument("-m", "--min-focal-freq",
                    dest="min_focal_freq",
                    type=int,
                    default=0.25,
                    help="Minimum frequency of derived focal variant.")
parser.add_argument("-M", "--max-focal-freq",
                    dest="max_focal_freq",
                    type=int,
                    default=0.95,
                    help="Maximum frequency of derived focal variant.")
group1 = parser.add_argument_group('cosi formatted input', 'Give files formatted as returned from cosi.')                    
group1.add_argument("--hap", dest="hap_file",
                    type=str,
                    help="File name for haplotype file.")
group1.add_argument("--pos", dest="pos_file",
                    type=str,
                    help="File name for map file.")
group2 = parser.add_argument_group('Phased VCF input', 'Give input data as phased VCF, ancestral sequence and recombination map')                    
group2.add_argument("--vcf", dest="vcf_file",
                    type=str,
                    help="File name for VCF file.")
group2.add_argument("--ancestral", dest="ancestral_file",
                    type=str,
                    help="File name for ancestral sequence file. (E.g. from http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/supporting/ancestral_alignments/)")                    
group2.add_argument("--recmap", dest="recmap_file",
                    type=str,
                    help="File name for recombination map file. Physical positions will be converted to genetic positions using this map. (E.g. from http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/)")                    

args = parser.parse_args()

if bool(args.vcf_file) == bool(args.hap_file):
    print("Specify either VCF or iHS input parameters")
    sys.exit()
if args.hap_file:
    data = read_processed_data(args.hap_file, args.pos_file)
else:
    data = read_phased_vcf_data(args.vcf_file, args.ancestral_file, args.recmap_file)

results = list()
for i, (phys_pos, gen_pos, ba) in enumerate(data):
    # gen_pos is same as phys_pos unless a recombination map is provided

#    freq = ba.count() / ba.length()
    freq = ba.count() / len(ba)
    if freq < args.min_focal_freq or freq > args.max_focal_freq:
        continue

    #print(freq)

    gen_dist = bool(args.recmap_file)
    sq_freqs = sq_freq_pairs(data, i, 0, args.window_size/2, args.max_ancest_freq, args.min_tot_freq, gen_dist) + \
        sq_freq_pairs(data, i, len(data), args.window_size/2, args.max_ancest_freq, args.min_tot_freq, gen_dist)

    if sq_freqs:
#        num = sum(x - x + y for x,y in sq_freqs)
#        den = sum(y - y + x for x,y in sq_freqs) + 0.001
        num = sum(x - x + y + 1 for x,y in sq_freqs)
        den = sum(y - y + x + 1 for x,y in sq_freqs) # redefine to add one to get rid of blowup issue introduced by adding 0.001 to denominator
        hapdaf = num / den
        results.append((phys_pos, hapdaf, freq))
    
if args.out_file:
    out = open(args.out_file, 'w')
else:
    out = sys.stdout
print('#', 'phys_pos', 'gen_pos', 'hapdaf', 'deriv_freq', file=out)
for tup in results:
    print(*tup, file=out)
