

# dict of amino acid name codes
aa = {'Ala':'A',
        'Arg':'R',
        'Asn':'N',
        'Asp':'D',
        'Cys':'C',
        'Glu':'E',
        'Gln':'Q',
        'Gly':'G',
        'His':'H',
        'Ile':'I',
        'Leu':'L',
        'Lys':'K',
        'Met':'M',
        'Phe':'F',
        'Pro':'P',
        'Ser':'S',
        'Thr':'T',
        'Trp':'W',
        'Tyr':'Y',
        'Val':'V'}

# parameters for chou-fasman from the pdf
parameters = """
Glu 1.53 Met 1.67
Ala 1.45 Val 1.65
Leu 1.34 Ile 1.60
His 1.24 Cys 1.30
Met 1.20 Tyr 1.29
Gln 1.17 Phe 1.28
Trp 1.14 Gln 1.23
Val 1.14 Leu 1.22
Phe 1.12 Thr 1.20
Lys 1.07 Trp 1.19
Ile 1.00 Ala 0.97
Asp 0.98 Arg 0.90
Thr 0.82 Gly 0.81
Ser 0.79 Asp 0.80
Arg 0.79 Lys 0.74
Cys 0.77 Ser 0.72
Asn 0.73 His 0.71
Tyr 0.61 Asn 0.65
Pro 0.59 Pro 0.62
Gly 0.53 Glu 0.26
""".split()

# dict for storing Chou-Fasman parameters
pa = {}
pb = {}
for i in range(0,len(parameters),4):
    pa[aa[parameters[i]]] = float(parameters[i+1])
    pb[aa[parameters[i+2]]] = float(parameters[i+3])
pa = dict(sorted( pa.items(),key= lambda x:-x[1]))
pb = dict(sorted( pb.items(),key= lambda x:-x[1]))


def chou_fasman_alpha(protein_sequence):
    """
    This function takes a protein sequence as input and returns the sites which are likely to be alpha helices.
    returns a boolean array which is true at index with likelihood of alpha helix
    """
    # boolean array
    h = [False]*len(protein_sequence)
    window_size,residues_needed = 6,4
    # main loop for going over all windows starting from index l
    for l in range(0,len(protein_sequence)-window_size+1):
        c = 0
        # count the number of residues with propensity greater than 1
        for i in range(l,l+window_size):
            if pa[protein_sequence[i]] > 1:
                c += 1
        if c>=residues_needed:
            # extending the window to the right
            r = l+window_size
            while r < len(protein_sequence):
                sum = 0
                # calculating the sum of propensity of residues in the extension
                for a in protein_sequence[r-3:r+1]:
                    sum += pa[a]
                if sum >= 4:
                    r += 1
                else:
                    break
            # extending the window to the left
            while l > 0:
                sum = 0
                for a in protein_sequence[l-1:l+3]:
                    sum += pa[a]
                if sum >= 4:
                    l -= 1
                else:
                    break
            # marking the residues in the window as alpha helices
            for i in range(l,r):
                h[i] = True
    return h

def chou_fasman_beta(protein_sequence):
    """
    This function takes a protein sequence as input and returns the sites which are likely to be beta strands.
    returns a boolean array which is true at index with likelihood of beta strand
    """
    # boolean array
    s = [False]*len(protein_sequence)
    window_size,residues_needed = 5,3
    # main loop for going over all windows starting from index l
    for l in range(0,len(protein_sequence)-window_size):
        c = 0
        # count the number of residues with propensity greater than 1
        for i in range(l,l+window_size):
            if pb[protein_sequence[i]] > 1:
                c += 1
        if c>=residues_needed:
            # extending the window to the right
            r = l+window_size
            while r < len(protein_sequence):
                sum = 0
                for a in protein_sequence[r-3:r+1]:
                    sum += pb[a]
                if sum >= 4:
                    r += 1
                else:
                    break
            # extending the window to the left
            while l > 0:
                sum = 0
                for a in protein_sequence[l-1:l+3]:
                    sum += pb[a]
                if sum >= 4:
                    l -= 1
                else:
                    break
            # marking the residues in the window as beta strands
            for i in range(l,r):
                s[i] = True
    return s


def chou_fasman_algorithm(protein_sequence):
    """
    This function takes a protein sequence as input and returns the Chou-Fasman
    prediction of the secondary structure of the protein.
    """
    # getting the alpha helix and beta strand prediction
    h = chou_fasman_alpha(protein_sequence)
    s = chou_fasman_beta(protein_sequence)
    # answer string
    a = ''
    i = 0
    # main loop for going over all residues
    while i<len(h):
        # over-lapping alpha helix and beta strand. resolving conflicts
        if h[i] and s[i]:
            c=0
            ph,ps = 0,0 
            while i<len(h) and h[i] and s[i]:
                ph+=pa[protein_sequence[i]]
                ps+=pb[protein_sequence[i]]
                i+=1
                c+=1
            # set overlap with structure having greater sum of propensity
            if ph>ps:
                a=a+'H'*c
            else:
                a=a+'S'*c
        elif h[i] and not s[i]:
            i += 1
            a=a+'H'
        elif not h[i] and s[i]:
            i += 1
            a=a+'S'
        else:
            a=a+'_'
            i += 1
    # converting boolean array to string
    h=''.join('h' if i  else '_' for i in h)
    s=''.join('s' if i  else '_' for i in s)
    
    # printing the answer
    print('alpha and beta sites:\n')
    for i in range(0,len(protein_sequence),50):
        print(i+1,'\t',protein_sequence[i:i+50])
        print('\t',h[i:i+50])
        print('\t',s[i:i+50])
        print()
    print('\nfinal output after resolving conflicts:\n')
    for i in range(0,len(protein_sequence),50):
        print(i+1,'\t',protein_sequence[i:i+50])
        print('\t',a[i:i+50])
        # print('\t',h[i:i+50])
        # print('\t',s[i:i+50])
        print()
    return a

# function to match chou-fasman algorithm with the stride result
def match(c,s):
    m = 0
    a = ''
    for i,j in zip(c,s):
        if i==j:
            m+=1
            a = a+'.'
        else:
            a=a+'*'
    
    for i in range(0,len(c),50):
        print('chou-f','\t',i+1,'\t',c[i:i+50])
        print('stride','\t','\t',s[i:i+50])
        print('mismatch','\t',a[i:i+50])
        print()
        
    print('matches','\t',m,'/',len(c))
    
if __name__=='__main__':
    # print('pa:',pa)
    # print('pb:',pb)
    print()
    protein_sequence = 'SGFRKMAFPSGKVEGCMVQVTCGTTTLNGLWLDDTVYCPRHVICTAEDMLNPNYEDLLIRKSNHSFLVQAGNVQLRVIGHSMQNCLLRLKVDTSNPKTPKYKFVRIQPGQTFSVLACYNGSPSGVYQCAMRPNHTIKGSFLNGSCGSVGF'
    print('Answer 1:','\n')
    c = chou_fasman_algorithm(protein_sequence)
    
    print('Answer 2:\n\ncomparing with stride result\n')
    
    stride = 'TTTT     HHHHHH EEEEEETTEEEEEEEETTEEEEEGGGG  HHHHH   HHHHHHH  GGG EEEETTEEE EEEEEEETTEEEEEE   TTTT        TTTEEEEEEEEETTEEEEEEEEEETTTT B    TTTTTTTEE '
    stride = stride.replace(' ','_')
    stride = stride.replace('E','S')
    match(c,stride)
    