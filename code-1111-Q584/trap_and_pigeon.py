import sys, itertools
sys.path.insert(0,__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)),'..','code-1023c2-Q569'))
from puzzles import TILES, tile_edges, D, edge_key
from star import partition, tri_edges
# trap placed with wall edges on a=0: trap steps (1,0),(2,1),(4,0),(4,1),(0,1); find anchor with both dir-4 edges at a=0
T=TILES['trap']
for a in range(0,4):
  for b in range(-3,4):
    E=tile_edges(T,(a,b))
    wall=[(e,l) for e,l in E if e[1]==1 and e[0][0]==0]
    if len(wall)==2:
        bulk=[(e,l) for e,l in E if not (e[1]==1 and e[0][0]==0)]
        print('trap anchor',(a,b),' wall legs',wall,' bulk legs (class,label):',sorted((e[1],l) for e,l in bulk))
        tris=[('U',0,b-1+0),('D',0,b-1),('U',0,b)]
        # identify the three unit triangles covered (v=(0,b+?)) by matching edges
        cand=[('U',0,y) for y in range(b-3,b+3)]+[('D',0,y) for y in range(b-3,b+3)]
        Eset={e for e,_ in E}
        TR=None
        for c in itertools.combinations(cand,3):
            es=[x for t in c for x in tri_edges(*t)]
            bd={x for x in es if es.count(x)==1}
            if bd==Eset: TR=c
        print('trap = triangles',TR)
        B=dict(E)
        print('KZJ-triangle fillings of trap (no K / K-up / K-dn / both):',
              [partition(list(TR),B,u,d) for u,d in ((0,0),(1,0),(0,1),(1,1))])
        break
# pigeonhole: both sides of the type-B braid relation, all K orientations
L=[('A_lo','B_lo'),('B_hi','A_lo')]; R=[('A_hi','B_lo'),('A_hi','B_hi')]
for oa,ob in itertools.product(((0,2),(2,0)),repeat=2):
    cl={'A_lo':oa[0],'A_hi':oa[1],'B_lo':ob[0],'B_hi':ob[1]}
    print('A',oa,'B',ob,' LHS same-class crossings:',[p for p in L if cl[p[0]]==cl[p[1]]],
          ' RHS:',[p for p in R if cl[p[0]]==cl[p[1]]])
