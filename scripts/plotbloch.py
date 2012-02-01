#!/usr/bin/env python

import qutip
import sexpr
import numpy
import sys

def getqmem(str):
    return sexpr.str2sexpr(str)[0]
def getqids(exp):
    return exp[0]
def getnodes(exp):
    return exp[1]

def parse_file(file):
    qids = []
    q = {}
    f = open(file)
    s = f.read()
    f.close()
    parsed = getqmem(s)

    for qid in getqids(parsed):
        qids.append( int(qid) )
    
    for node in getnodes(parsed):
        q[int(node[0])] = complex(node[1].replace('i','j'))

    a = numpy.zeros((2**len(qids),1),complex)

    for base,amp in q.iteritems():
        a[base] = [amp]
    return a

def main():
    b = qutip.Bloch()
    for file in sys.argv[1:]:
        a = parse_file(file)
        print "plotting: "
        print a
        b.add_states(qutip.Qobj(a))
    b.show()

main()
