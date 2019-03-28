from goatools.associations import read_ncbi_gene2go
from goatools.obo_parser import GODag
from goatools.go_enrichment import GOEnrichmentStudy
from goatools.godag_plot import plot_results
#from collections import defaultdict, namedtuple
import os


def build_network(taxid=511145, alpha=0.05):
    curdir = os.path.dirname(os.path.abspath(__file__))
    gene2go = curdir + '/../data/gene2go'
    taxid2asscs = {}#defaultdict(lambda: defaultdict(lambda: defaultdict(set)))
    # Associate ncbi ID to GO terms for the taxid
    read_ncbi_gene2go(gene2go, taxids=[taxid], taxid2asscs=taxid2asscs)

    # Create the goatools datastructure that holds name relationships
    GeneID2nt = make_geneID2nt(taxid)
    # Get GeneID2GOs association for current species
    geneid2gos = taxid2asscs[taxid]['GeneID2GOs']
    # Load in the ontology graph
    obodag = GODag(curdir+'/../data/go-basic.obo')
    goeaobj = make_goa_object(GeneID2nt, geneid2gos, obodag, alpha)
    return goeaobj, GeneID2nt


def id_to_genelist(GeneID2nt, gene_df):
    gene_names = list(gene_df.index)
    go_names = [x[5] for x in GeneID2nt.values()]
    set_names = {}
    for gene in go_names:
        for degene in gene_names:
            if gene == degene:
                set_names[gene] = [(gene, degene)]

    de_data = {}
    for e, row in gene_df.iterrows():
        dat = row.values
        de_gene = e.lower()
        for gene in set_names:
            if gene.lower() == de_gene:
                # gene, logFC, Pval, FDR
                de_data[gene] = [dat[0], dat[3], dat[4]]
    return de_data


def get_ID2nt(GeneID2nt, genelist, sign):
    geneid2symbol = {}
    geneids_study = []
    for geneid in GeneID2nt:
        genesymbol = GeneID2nt[geneid][5]
        if genesymbol in genelist[sign]:
            if geneid not in geneid2symbol:
                geneid2symbol[geneid] = [genesymbol]
            geneid2symbol[geneid].append(genesymbol)
            geneids_study.append(geneid)
    return geneids_study


def make_geneID2nt(taxid):
    curdir = os.path.dirname(os.path.abspath(__file__))
    with open(curdir + '/../data/Bacteria.gene_info') as of:
        data = of.readlines()

    species = get_species_from_taxid(taxid)
    #NtData = namedtuple('NtData',
    #                    ['tax_id', 'Org_name', 'GeneID', 'CurrentID', 'Status', 'Symbol', 'Aliases', 'description',
    #                    'other_designations', 'map_location', 'chromosome', 'genomic_nucleotide_accession_version',
    #                    'start_position_on_the_genomic_accession', 'end_position_on_the_genomic_accession',
    #                    'orientation', 'exon_count', 'OMIM', 'no_hdr0'])

    GeneID2nt = {}
    for genedata in data[1:]:
        k = genedata.strip().split()
        geneid = k[1]
        file_taxid = k[0]
        if int(file_taxid) == int(taxid):
            """GeneID2nt[int(geneid)] = NtData(tax_id=k[0],
                                            Org_name=species,
                                            GeneID=k[1],
                                            CurrentID=0,
                                            Status=k[12],
                                            Symbol=k[2],
                                            Aliases=k[4].replace('|', ','),
                                            description=k[8],
                                            other_designations=k[13],
                                            map_location=k[7],
                                            chromosome=k[6],
                                            genomic_nucleotide_accession_version=None,
                                            start_position_on_the_genomic_accession=None,
                                            end_position_on_the_genomic_accession=None,
                                            orientation='',
                                            exon_count=None,
                                            OMIM='',
                                            no_hdr0='')
            """
            GeneID2nt[int(geneid)] = [None, None, None, None, None, k[2]]
    return GeneID2nt


def make_goa_object(GeneID2nt, geneid2gos, obodag, alpha=0.05):
    goeaobj = GOEnrichmentStudy(
        GeneID2nt.keys(),
        geneid2gos,
        obodag,
        propagate_counts=False,
        alpha=alpha,
        methods=['fdr_bh'])
    return goeaobj


def get_species_from_taxid(taxid):
    curdir = os.path.dirname(os.path.abspath(__file__))
    with open(curdir + '/../data/scientific_taxids.txt') as taxes:
        data = taxes.readlines()
    for line in data:
        taxid_tmp = line.split()[0]
        species = ' '.join(line.split()[1:]).strip()
        if int(taxid) == int(taxid_tmp):
            return species
    if not species:
        return 'Unknown'
