'''
goal of this script is to convert a networkx formated network to an
aracne-readable adjacency matrix
'''
import pickle, networkx,re

filepath='/Users/sgosline/Desktop/tf2gene.pkl'

net = pickle.load(open(filepath,'r'))

outfile=open('tf2gene.adj','w')

for node in net.nodes():
    sucs = net.successors(node)
    if(len(sucs)==0):
        continue

    outfile.write(node)
    for suc in sucs:
        edge_weight=net.get_edge_data(node,suc)['weight']
        outfile.write('\t'+re.sub('mrna','',suc)+'\t'+str(edge_weight))
    outfile.write('\n')

outfile.close()
