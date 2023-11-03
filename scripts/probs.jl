# Beta

pac = 0.14
pag = 0.05
pat = 0.05
#@assert pac+pag+pat == 0.25

pca = 0.13
pcg = 0.02
pct = 0.04
@assert pca + pcg + pct == 0.19

pga = 0.04
pgc = 0.08
pgt = 0.12
@assert pga + pgc + pgt == 0.24

pta = 0.08
ptc = 0.15
ptg = 0.09
#@assert pta + ptc + ptg == 0.33

@assert pca + pga + pta == 0.25
#@assert pac + ptc + pgc == 0.38
@assert pag + pcg + ptg == 0.16
@assert pat + pct + pgt == 0.21

pc = pg = 0.205
pa = pt = (1-0.41)/2
@assert pc + pg + pa + pt == 1.0

# purines: A, G
# Pyrimidines: C, T 

# purine to purine, A->G and G->A 
ppupu = pa*pag + pg*pga 

# purine to pyrimidine
# G->C, G->T, A->C, A->T  
ppupy = pg*pgc + pg*pgt + pa*pac + pa*pat  

# pyrimidine to purine 
# C->G, C->A, T->G, T->A 
ppypu = (pc*pg*pcg + pc*pa*pca + pt*pg*ptg + pt*pa*pta) / totalprob  

# pyrimidine to pyrimidine
# C->T, T->C 
ppypy = (pc*pt*pct + pt*pc*ptc) / totalprob  

@assert ppupu + ppupy + ppypu + ppypy ≈ 1.0



# Helicobacter

pac = 0.25
pag = 0.05
pat = 0.10
#@assert pac+pag+pat == 0.41

pca = 0.07
pcg = 0.02
pct = 0.06
@assert pca + pcg + pct == 0.15

pga = 0.02
pgc = 0.04
pgt = 0.18
@assert pga + pgc + pgt == 0.24

pta = 0.04
ptc = 0.10
ptg = 0.06
@assert pta + ptc + ptg == 0.20

@assert pca + pga + pta == 0.14
@assert pac + ptc + pgc == 0.39
@assert pag + pcg + ptg == 0.14
@assert pat + pct + pgt == 0.34

pc = pg = 0.22
pa = pt = 0.28

