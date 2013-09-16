# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
import argparse
import random

ALL_PROPERLY_ALIGNED = 0x2
IM1=9
IP1=10
IC1=11


def mapping_parser(m):
    '''
    Parse a read in SAM format, return a dictionary with filled in fields of interest.
    '''
    if isinstance(m, str):
        m = m.strip().split('\t')
        d = {}
        d['flag'] = int(m[1])   # flags
        d['chr'] = m[2]         # chr
        d['pos'] = int(m[3])    # pos
        d['mapq'] = int(m[4])   # mapping quality
        d['cigar'] = m[5]       # cigar string
        d['seq'] = m[9]         # sequence
        d['qual'] = m[10]       # sequencing quality

    return d

def parse_cigar_string(s):
    '''
    Parse given CIGAR string to a list of operators.
    '''
    res = []
    crt_len = ''
    i = 0
    while i < len(s):
        if str.isdigit(s[i]):
            crt_len += s[i]
        else:
            res.append([s[i], int(crt_len)])
            crt_len = ''
        i += 1  
    return res

def upper_bound(arr, ch, pos, q, p):

    if (p-q==1):
        if arr[q]['chr']==ch and arr[q]['pos']==pos:
            return q
        else:
            return p
    
    mid=(q+p)/2;

    if arr[mid]['chr']>ch:
        return upper_bound(arr, ch, pos, q, mid);
    elif arr[mid]['chr']<ch:
        return upper_bound(arr, ch, pos, mid, p);

    if arr[mid]['pos']>pos:
        return upper_bound(arr, ch, pos, q, mid);
    else:
        return upper_bound(arr, ch, pos, mid, p);

 
def main():
    parser = argparse.ArgumentParser(description='Filter SNP positions by call quality and min. coverage. Awaits filenames for M, P .vcf files, and M, P .sam files.')
    parser.add_argument('filenames', type=str, nargs='+', help='paths to .vcf files with M, P SNPs and to corresponding .sam files')
    args = parser.parse_args()



    reads_file=open(args.filenames[0], "r")
    snips_file=open(args.filenames[1], "r")
    goalHaplotype=args.filenames[2]


    snips=[]

    nsnips=0

    # copy snips from file to RAM
    for line in snips_file:
        if line[0]=='#':
            continue
        if isinstance(line, str):
            line = line.strip().split('\t')
            haplos=[line[3], line[4]]
            phasing=line[IM1].strip().split(':')[0].strip().split('|');
            if (phasing[0] != phasing[1]):
                snips.append({"pos":int(line[1]),
                              "chr":('chr'+line[0]),
                              "HA": (haplos[int(phasing[0])]),
                              "HB": (haplos[int(phasing[1])])
                             })
#            line = line.strip().split('\t')
#            snips.append({"pos":int(line[2].strip().split(':')[1]),
#                          "chr":(line[2].strip().split(':')[0]),
#                          "HA": (line[5][0]),
#                          "HB": (line[6][0])
#                         })
#    print upper_bound(snips, 'chr20', 10181440, 0, len(snips));
#    print upper_bound(snips, 'chr20', 10281440, 0, len(snips));
#    exit()
    #for each read:
    for read in reads_file:
        parsed_read= mapping_parser(read)
        if parsed_read['flag'] & ALL_PROPERLY_ALIGNED == 0: continue 
        if parsed_read['mapq'] < 20: continue 

        pos_qr = 0
        pos_db = parsed_read['pos']
        op = parse_cigar_string(parsed_read['cigar'])

        snip = upper_bound(snips, parsed_read['chr'], parsed_read['pos'], 0, len(snips));
        countA=0
        countB=0
        countMis=0

        for o in op:
            while pos_db>snips[snip]['pos'] and snip<len(snips) and snips[snip]['chr']==parsed_read['chr']:
                snip += 1
            if snip==len(snips) or snips[snip]['chr']!=parsed_read['chr']:
                break

            if o[0] == 'H': continue
            elif o[0] in 'SI': pos_qr += o[1]
            elif o[0] in 'ND': pos_db += o[1]
            elif o[0] in 'M=X':
                for i in range(o[1]):
                    #if the read position is of sufficient quality, record this info
                    if pos_db==snips[snip]['pos'] and ord(parsed_read['qual'][pos_qr]) >= 33+20:
                        if snips[snip]['HA']==parsed_read['seq'][pos_qr]:
                            countA += 1
                        elif snips[snip]['HB']==parsed_read['seq'][pos_qr]:
                            countB += 1
                        else:
                            countMis += 1
                    pos_db += 1
                    pos_qr += 1

        if countA>countB and goalHaplotype=='A':
            print read,
        if countA<countB and goalHaplotype=='B':
            print read,
#        if countA==countB and random.randrange(0,2)==1:
#            print read,


if __name__ == '__main__':
    main()
