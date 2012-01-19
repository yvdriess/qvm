from qutip import Bloch, Qobj
from sexpr import str2sexpr
import numpy

def getqmem(str):
    return str2sexpr(str)[0]
def getqids(exp):
    return exp[0]
def getnodes(exp):
    return exp[1]

qids = []
q = {}

f = open("out")
s = f.read()
f.close()

parsed = getqmem(s)

for qid in getqids(parsed):
    qids.append( int(qid) )
    
for node in getnodes(parsed):
    q[int(node[0])] = complex(node[1].replace('i','j'))

a = numpy.zeros((2**len(qids),1),complex)

for base,amp in q.iteritems():
    print base
    print amp
    a[base] = [amp]

b = Bloch()
b.add_states(Qobj(a))
b.show()
