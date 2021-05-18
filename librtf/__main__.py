import argparse
from librtf.rtutils import *
from librtf.sequence  import *
from librtf.tree  import *
from librtf.rtcounts  import *
from librtf.rtdistribution  import *




parser = argparse.ArgumentParser(description='Return time frquency based phylogeny -beta-',prog='librtf')
#positional arguments
parser.add_argument('filename', metavar="filename", help="path to input file")
parser.add_argument('ksize',metavar="ksize", help="set k-mer size (or kmer range ex: 1:10)")
parser.add_argument('type', metavar="type", help = 'specify the sequence type. Protein:P Nucleotide:N')
#optional argumnets
parser.add_argument('-p','--positions', metavar=' ', type = str, required = False, help = 'binary pattern which defines the on and')
parser.add_argument('-o', '--out', metavar=' ', type = str, required = False, help = 'path to output')
parser.add_argument('-index', metavar=' ', type = str, required = False, help = 'path to a file containing a pre-specified set of k-mers')
parser.add_argument("-d", metavar=' ',type = str, required = False, help='use an alternative distance function')
parser.add_argument("-ad", metavar=' ', type = str, required = False, help="calculate average distance using n distance functions")
parser.add_argument("-nj", help="use nj instead of rapidnj", action="store_true", required = False)
parser.add_argument("-w", help="write k-mer index to a file", action ="store_true", required = False)
parser.add_argument("-rtd",  help="use return time distribution method", action="store_true", required = False)
parser.add_argument("-write", help="select alternative distance function", action ="store_true", required = False)
args = parser.parse_args()

print(f'\nreading.... {args.filename}')

sizes = args.ksize.split(':')

if len(sizes) > 1:
    start,end = int(sizes[0]),int(sizes[1])+1

else:
    start,end = int(sizes[0]), int(sizes[0])+1

#check k-mer length and positions vector length
if args.positions:
    if end - start != 1 :
        raise SystemExit('gapped method can only be used for a specific k-mer size')

    if len(args.positions) != start:
        raise SystemExit('binary vector length does not match kmer size ')


#return times based phylogeny
if args.rtd:
    for k in range(start,end):

        X_mean,X_sd,y_,sorted_global = ReturnTimeDistribution.cal_rtd_vectors(Sequences(args.filename, format_='fasta', type_ = args.type), k = k,  positions = args.positions if args.positions else None)

        if len(sorted_global) != 0:

            if args.write:
                WriteDicToFile(sorted_global)
            
            if args.out:
                GenerateTree.ReturnTimeDistributionTree(X_mean, X_sd, y_,  args.out, nj=args.nj if args.nj else False)
            else:
                Fname = f'RTD-k{k}-{Datetime()}'
                GenerateTree.ReturnTimeDistributionTree(X_mean, X_sd, y_,  Fname, nj=args.nj if args.nj else False)
        else:
            print(f'No common k-mers were found for the given ksize {k}')       
#return time frequency based phylogeny   
else:   
    for k in range(start,end):
        X, y_, sorted_global = ReturnTimeCounts.cal_rt_vectors(Sequences(args.filename, format_='fasta', type_=args.type), k = k, positions = args.positions if args.positions else None)
        if len(sorted_global) != 0:    
            if args.write:
                WriteDicToFile(sorted_global)
            #X, y_ = dna2vec.cal_dna2vec_vectors(Sequences(args.filename,'fasta'))
            #GenerateTree.ReturnTimeCountsTree(X, y_, args.out)
            #m0 = process.memory_info().rss/(1024*1024) # in bytes 
            #X, y_, sorted_global = ReturnTimeCounts.cal_rt_vectors(Sequences('Hep.txt','fasta'), k = int(1))
            #m1 = process.memory_info().rss/(1024*1024)  # in bytes 
            if args.out:
                GenerateTree.ReturnTimeCountsTree(X, y_, args.out, nj=args.nj if args.nj else False)
            else:
                Fname = f'RTF-k{k}-{Datetime()}'
                GenerateTree.ReturnTimeCountsTree(X, y_, Fname, nj=args.nj if args.nj else False)
        else:
            print(f'No common k-mers were found for the given ksize {k}') 
    #m2 = process.memory_info().rss/(1024*1024)
    #print(process, m2 - m0)

