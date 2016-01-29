from sys import argv

code=[ [['G']*4,['E']*2+['D']*2,['A']*4,['V']*4] ,
       [['R']*2+['S']*2,['K']*2+['N']*2,['T']*4,['M']+['I']*3] , 
       [['R']*4,['Q']*2+['H']*2,['P']*4,['L']*4] , 
       [['W','opal']+['C']*2,['amber','ochre']+['Y']*2,['S']*4,['L']*2+['F']*2] ]
order={'A':1 , 'C':2 , 'G':0 , 'U':3 }



def detectstart(ARNm):
    startposition=-1
    for i in xrange(len(ARNm)):
        a=ARNm[i:i+3]
        if a=="AUG" and startposition==-1:
            startposition=i
    return startposition

def detectstop(ARNm):
    stopposition=-1
    for i in xrange(startposition, len(ARNm)-2,3):
        a=ARNm[i:i+3]
        if stopposition==-1 and (a=="UAA" or a=="UAG" or a=="UGA"):
            stopposition=i
    return stopposition

def translatecodon(c):
    return code[ order[c[0]] ][ order[c[1]] ][ order[c[2]] ]

if len(argv)<2:
    print "Syntaxe: python Ribosome.py AGAUGGGCGCUUUGAUUGCUGUAGUAUAAU"
    exit("Pas de sequence a traduire")

    
ARNm=argv[1]

startposition=detectstart(ARNm)
stopposition=detectstop(ARNm)

if startposition==-1:
    print "Erreur: Pas de codon start."
    exit("Pas de codon start")
if stopposition==-1:
    print "Erreur: Pas de codon stop."
    exit("Pas de codon stop")

print "Codon start detecte:",ARNm[startposition : startposition + 3],"en position",startposition+1
print "Codon stop detecte:",ARNm[stopposition : stopposition +3],'(',translatecodon(ARNm[stopposition : stopposition +3]),"), la traduction s\'arretera en",stopposition

TR=ARNm[startposition : stopposition]

Prot=""
for i in range(0,len(TR)-2,3):
    codon=TR[i:i+3]
    AA=translatecodon(codon)
    Prot=Prot + AA

print "Proteine traduite:",Prot
