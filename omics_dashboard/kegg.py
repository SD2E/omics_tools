from itertools import chain
import pickle
from urllib.request import urlopen
from urllib.error import HTTPError
import os


def symbol_to_entrez(taxid=511145):
    curdir = os.path.dirname(os.path.abspath(__file__))
    with open(curdir + '/../data/Bacteria.gene_info') as of:
        data = of.readlines()
    symbol_to_entrez = {}
    #symbol_to_geneID = {}
    for line in data[1:]:
        line = line.split()
        if line[0] == str(taxid):
            #ncbi = line[1]
            symbol = line[2]
            locustag = line[3]
            # symbol_to_geneID[symbol] = ncbi
            symbol_to_entrez[symbol] = locustag
    return symbol_to_entrez


# Using Request and geneIDconv from 'sharepathway' package
# https://pypi.org/project/sharepathway/
# Author: Guipeng, Li guipeng.lee@gmail.com
def Request(*args, **kwargs):
    """return and save the blob of data that is returned
    from kegg without caring to the format"""
    downloaded = kwargs.get('force', False)
    save = kwargs.get('force', True)

    # so you can copy paste from kegg
    '''URL form
        http://rest.kegg.jp/<operation>/<argument>/<argument2 or option>
        <operation> = info | list | find | get | conv | link
        <argument> = <database> | <dbentries>
    '''
    args = list(chain.from_iterable(a.split('/') for a in args))
    args = [a for a in args if a]
    request = 'http://rest.kegg.jp/' + "/".join(args)
    print(("KEGG API: " + request))
    filename = "KEGG_" + "_".join(args)
    try:
        if downloaded:
            raise IOError()
        print(("loading the cached file " + filename))
        with open(filename, 'rb') as f:
            data = pickle.load(f)
    except IOError:
        print("downloading the library,it may take some time")
        try:
            req = urlopen(request)
            data = req.read().decode()
            data = [tuple(d.split('\t')) for d in data.split('\n')][:-1]
            if save:
                with open(filename, 'wb') as f:
                    print(("saving the file to " + filename))
                    pickle.dump(data, f)
        # clean the error stacktrace
        except HTTPError as e:
            raise e
    return data


def geneIDconv(*args, **kwargs):
    species = kwargs.get('species', 'eco')
    ngid = kwargs.get('NCBI-GeneID', 'ncbi-geneid')
    # KEGG ID conv from NCBI-GeneID to KEGG ID, though most of them are the same
    keggd = Request('conv', species, ngid)
    Kgid = {}
    for line in keggd:
        Kgid[line[0].split(':')[1]] = line[1]
    return Kgid
