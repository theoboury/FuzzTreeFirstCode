import varnaapi, pickle

def traduction(label_char):
    resu = []
    if label_char[0] == 'C':
        resu.append('cis')
    else:
        resu.append('trans')
    for i in range(1, 3): 
        if label_char[i] == 'W':
            resu.append('wc')
        elif label_char[i] == 'H':
            resu.append('h')
        else:
            resu.append('s')
    return resu[1], resu[2], resu[0]

def draw(GT, nodes_target, mapping):
    #nodes_target should be ordered in lexicographic order
    secondary_structures = [(nodes_target.index(i), nodes_target.index(j)) for (i, j, t) in GT.edges.data() if t['label'] == "CWW"]
    noncanonicaledges = [(nodes_target.index(i), nodes_target.index(j), t['label']) for (i, j, t) in GT.edges.data() if (nodes_target.index(i), nodes_target.index(j)) not in secondary_structures]
    sequence  = []
    if nodes_target != []:
        strain = nodes_target[0][0]
    for node in  nodes_target:
        (i, t) = [(i, t) for (i, t) in GT.nodes.data() if i == node][0]
        if node[0] != strain:
            sequence += "&"
            strain = node[0]
        sequence += t['nc']
    sequence = "".join(sequence)
    v  = varnaapi.Structure(sequence, secondary_structures)
    for (i,j,t) in noncanonicaledges:
        e5,e3,ster = traduction(t)
        v.add_aux_BP( i, j, edge5=e5, edge3=e3, stericity=ster, color='blue', thickness=1) 
    v.savefig("varnatest.png")


with open("1Y27.nxpickle",'rb') as fT:
    GT = pickle.load(fT)

nodes_target = list(GT.nodes())
nodes_target.sort()
draw(GT, nodes_target, [])


