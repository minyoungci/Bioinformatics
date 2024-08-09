def validate_dna(dna_seq):
    seqm = dna_seq.upper()
    valid = seqm.count("A") + seqm.count("T") + seqm.count("G") + seqm.count("C")
    if valid == len(seqm) : return True
    else: return False

def frequency(seq):
    dic = {}
    for s in seq.upper():
        if s in dic: dic[s] += 1
        else: dic[s] = 1
    return dic

def gc_content(dna_seq):
    gc_count = 0
    for s in dna_seq:
        if s in "GCgc" : gc_count += 1  
    return gc_count / len(dna_seq)

def gc_content_subseq(dna_seq, k=10):
    res = []
    for i in range(0, len(dna_seq)-k+1, k):
        subseq = dna_seq[i:i+k]
        gc = gc_content(subseq)
        res.append(gc)
    return res  

def transcription(dna_seq):
    assert validate_dna(dna_seq), "not a valid DNA sequence"    
    return dna_seq.upper().replace("T", "U")    

def reverse_complement(dna_seq):# DNA 서열의 역상보 서열을 반환
    assert validate_dna(dna_seq), "not a valid DNA sequence"    
    comp = ''
    for c in dna_seq.upper():
        if c == 'A': comp = 'T' + comp
        elif c == 'T': comp = 'A' + comp
        elif c == 'G': comp = 'C' + comp
        elif c == 'C': comp = 'G' + comp
    return comp

def translate_codon (cod):
    """Translates a codon into an aminoacid using an internal dictionary with the standard genetic code."""
    tc = {"GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A", 
      "TGT":"C", "TGC":"C",
      "GAT":"D", "GAC":"D",
      "GAA":"E", "GAG":"E",
      "TTT":"F", "TTC":"F",
      "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",
      "CAT":"H", "CAC":"H",
      "ATA":"I", "ATT":"I", "ATC":"I",
      "AAA":"K", "AAG":"K",
      "TTA":"L", "TTG":"L", "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
      "ATG":"M", "AAT":"N", "AAC":"N",
      "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
      "CAA":"Q", "CAG":"Q",
      "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R", "AGA":"R", "AGG":"R",
      "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S", "AGT":"S", "AGC":"S",
      "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
      "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
      "TGG":"W",
      "TAT":"Y", "TAC":"Y",
      "TAA":"_", "TAG":"_", "TGA":"_"}
    if cod in tc: return tc[cod]
    else: return None

def translate_seq (dna_seq, ini_pos = 0):
    assert validate_dna(dna_seq), "Invalid DNA sequence"
    seqm = dna_seq.upper()
    seq_aa = ""
    for pos in range(ini_pos,len(seqm)-2,3):
        cod = seqm[pos:pos+3]
        seq_aa += translate_codon(cod)
    return seq_aa

def codon_usage(dna_seq, aa):
    assert validate_dna(dna_seq), "Invalid DNA sequence"
    seqm = dna_seq.upper()
    dic = {}
    total = 0
    for i in range(0, len(seqm)-2, 3):
        cod = seqm[i:i+3]
        if translate_codon(cod) == aa:
            if cod in dic: 
                dic[cod] += 1
            else: dic[cod] = 1
            total += 1
    if total >0:
        for k in dic:
            dic[k] /= total
    return dic