def Motif (seq, motif):
    counter =[]
    for i in range(len(seq)-len(motif)+1):
        #A = any (seq [i:i+len(motif)]==motif)
        if seq [i:i+len(motif)]==motif:
            counter.append(i+1)
    return counter

ans = Motif ("GATTACAT","AT")
print (ans)


print (10<=5)
print (20>=5)
