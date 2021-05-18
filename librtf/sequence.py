from Bio import SeqIO 
import re

class Sequences:
    
    '''This is class extends SeqIO class available on biopython '''
    
    def __init__(self, filename, format_='fasta', type_=None):
        
        self.filename  = filename
        self.format_ = format_
        if type_ == 'N' or type_=='P':
            self.type = type_
        else:
            raise ValueError(f'{type} is not a valid sequence type')

        '''if type == 'N':
            #Nucleotide ambiguous codes
            self.ambi = ['Y', 'R', 'W', 'S', 'K', 'M', 'D', 'V', 'H', 'B', 'X', 'N']
        elif type == 'P':
            #Protein ambiguous codes
            self.ambi = ['X','B','Z','J']
        else:
            raise ValueError(f'{type} is not a valid sequence type')'''

    def __iter__(self):
        
        fasta_sequences = SeqIO.parse(open(self.filename),self.format_)
        for seq in fasta_sequences:

            sequence = str(seq.__dict__['_seq']).upper()
            id_ = seq.__dict__['id']

            #search for ambiguous characters in sequences
            if self.type == 'N':
                matchobj = re.search(r'[^ATGC]', sequence)
            else:
                matchobj = re.search(r'[^ACDEFGHIKLMNPQRSTVWY]', sequence)

            if matchobj:
                #prints a warning message when a sequence has an ambiguous character
                print('WARNING ! ' + id_ + ' has invalid character/s ', matchobj.group(),'therefore, will not be included in the analysis')
                pass

            else:
                yield sequence, id_
