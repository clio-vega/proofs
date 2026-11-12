from mosaic import *
def area(G,poly):
    z=[G.phi(P) for P in poly]; return 0.5*sum((z[i].conjugate()*z[(i+1)%len(z)]).imag for i in range(len(z)))
n=4
left,right,bottom=side_edges(n)
Ps=puzzles(n,{})
print('puzzles n=4 (all boundaries):',len(Ps))
for sign in (1,-1):
  for rule in (0,1):
    bad=0; okc=0
    for lab in Ps:
        if any(sum(lab[E(*e)] for e in side)!=sum(lab[E(*e)] for e in left) for side in (right,bottom)): continue
        M=build_mosaic(n,lab,sign,rule); G=Geo(sign)
        sq=[t for t in M['tiles'] if t[0]=='sq']
        for _,q in sq:
            z=[G.phi(P) for P in q]; d=((z[1]-z[0]).conjugate()*(z[3]-z[0])).real
            if abs(d)>1e-9: bad+=1
        hexv=[]
        for c in M['corners']: hexv+= [c['Xp'],c['X']]
        z=[G.phi(P) for P in hexv]
        cr=[((z[(i+1)%6]-z[i]).conjugate()*(z[(i+2)%6]-z[(i+1)%6])).imag for i in range(6)]
        conv= all(x<-1e-9 for x in cr) or all(x>1e-9 for x in cr)
        A=sum(abs(area(G,t[1])) for t in M['tiles'])+sum(abs(area(G,r[1])) for nest in M['nests'] for r in nest)
        okc+= conv and abs(A-abs(area(G,hexv)))<1e-6
    print('sign',sign,'rule',rule,'non-square rhombus pieces',bad,'convex & area-exact mosaics',okc)
