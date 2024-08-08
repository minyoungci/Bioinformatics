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