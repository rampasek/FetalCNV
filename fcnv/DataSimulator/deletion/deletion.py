# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
import argparse
import random

ALL_PROPERLY_ALIGNED = 0x2
sourceIndex={'IM1':9,'IP1':10,'IC1':11}

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
    fetusRate=float(args.filenames[2])
    goalHaplotype=args.filenames[3]

    snips_f=[]
    snips_m=[]

    nsnips=0

    # Copy snips from file to RAM
    for line in snips_file:
        if line[0]=='#':
            continue
        if isinstance(line, str):
            line = line.strip().split('\t')
            haplos=[line[3], line[4]]
            phasing_m=line[9].strip().split(':')[0].strip().split('|')
            phasing_f=line[11].strip().split(':')[0].strip().split('|')
            snips_m.append({"pos":int(line[1]),
                          "chr":('chr'+line[0]),
                          "HA": (haplos[int(phasing_m[0])]),
                          "HB": (haplos[int(phasing_m[1])])
                         })
            snips_f.append({"pos":int(line[1]),
                          "chr":('chr'+line[0]),
                          "HA": (haplos[int(phasing_f[0])]),
                          "HB": (haplos[int(phasing_f[1])])
                         })
    # Use this if the previous format for the phasing file is being used
#            line = line.strip().split('\t')
#            snips.append({"pos":int(line[2].strip().split(':')[1]),
#                          "chr":(line[2].strip().split(':')[0]),
#                          "HA": (line[5][0]),
#                          "HB": (line[6][0])
#                         })

    # Number of Snips in the region:
#   print upper_bound(snips, 'chr20', 10181440, 0, len(snips));
#   print upper_bound(snips, 'chr20', 10281440, 0, len(snips));
#   exit()

    # For each read count the snips in the read and check if it belongs to the target hapoltype:
    for read in reads_file:
        parsed_read= mapping_parser(read)

        # If the read is not aligned properly, ignore it
        if parsed_read['flag'] & ALL_PROPERLY_ALIGNED == 0: continue 
        if parsed_read['mapq'] < 20: continue 

        pos_qr = 0
        pos_db = parsed_read['pos']
        op = parse_cigar_string(parsed_read['cigar'])

        # Search for the start of the related snips
        snip = upper_bound(snips_f, parsed_read['chr'], parsed_read['pos'], 0, len(snips_f));

        count=0

        # Count the snips for each haplotype using the cigar string
        for o in op:
            while pos_db>snips_f[snip]['pos'] and snip<len(snips_f) and snips_f[snip]['chr']==parsed_read['chr']:
                snip += 1
            if snip==len(snips_f) or snips_f[snip]['chr']!=parsed_read['chr']:
                break

            correctHaplo=True
            if o[0] == 'H': continue
            elif o[0] in 'SI': pos_qr += o[1]
            elif o[0] in 'ND': pos_db += o[1]
            elif o[0] in 'M=X':
                for i in range(o[1]):
                    if pos_db==snips_f[snip]['pos'] and ord(parsed_read['qual'][pos_qr]) >= 33+20:
                        plasma=parsed_read['seq'][pos_qr]

                        if snips_f[snip]["H"+goalHaplotype]==plasma:
                            correctHaplo=False

                        if snips_f[snip]['HA']==plasma:
                            count += fetusRate/2.0

                        if snips_f[snip]['HB']==plasma:
                            count += fetusRate/2.0

                        if snips_m[snip]['HA']==plasma:
                            count += (1-fetusRate)/2.0

                        if snips_m[snip]['HB']==plasma:
                            count += (1-fetusRate)/2.0

                    pos_db += 1
                    pos_qr += 1

        if count==0:
            count=1
        rate=fetusRate/count
        # Find the haplotype for the read 
        if correctHaplo==False or random.random()>rate:
            print read,

if __name__ == '__main__':
    main()
