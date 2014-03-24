"""Prints ranks for a set of genes in the given sample.

Usage:
  ranks.py --sample=<path> --atlas=<path> --gene_list=<path> --tissue=<tissue>

Options:
  -h --help     Show this screen.
  --sample	path to sample tabbed-csv, containing at least 'hgnc_symbol' and 'rpkm' columns
  --atlas	path to atlas tabbed-csv, containing at least 'hgnc_symbol' and <tissue> columns
  --gene_list	path to file containing list of genes to measure ranks for, one per line
  --tissue	tissue name as in RNA-seq atlas (e.g. "liver")
"""


import csv
import docopt
import operator


args = docopt.docopt(__doc__)
sample_path = args.get('--sample')
atlas_path = args.get('--atlas')
gene_list_path = args.get('--gene_list')
tissue_name = args.get('--tissue')


def sample_ranks():
  print 'Sample counts'
  reader = csv.reader(open(sample_path, 'r'), delimiter='\t')
  headers = reader.next()
  hgnc_pos = headers.index('hgnc_symbol')
  rpkm_pos = headers.index('rpkm')
  return genes_sorted_by_rpkm(reader, hgnc_pos, rpkm_pos)


def atlas_ranks():
  print 'Atlas counts'
  reader = csv.reader(open(atlas_path, 'r'), delimiter='\t')
  headers = reader.next()
  hgnc_pos = headers.index('hgnc_symbol')
  rpkm_pos = headers.index(tissue_name)
  return genes_sorted_by_rpkm(reader, hgnc_pos, rpkm_pos)


# returns a list of genes sorted by rpkm, descending. for repeat
# genes, uses the first value seen in the file.
def genes_sorted_by_rpkm(reader, hgnc_pos, rpkm_pos):
  vals = {}
  empty_count = 0
  repeat_count = 0
  total = 0
  for row in reader:
    total += 1
    gene_name = row[hgnc_pos]
    if gene_name == '' or gene_name == 'NA':
      empty_count += 1
      continue
    # data has gene repeats, fine only if the values are the same
    rpkm = eval(row[rpkm_pos])
    if gene_name in vals:
      if abs(rpkm - vals[gene_name]) > 0.00001:
        repeat_count += 1
    else:
      vals[gene_name] = rpkm
  print '%d unnamed, %d repeated, %d total rows' % (empty_count, repeat_count, total)
  return [x[0] for x in sorted(
      vals.iteritems(), key=operator.itemgetter(1), reverse=True)]


def read_gene_list():
  return open(gene_list_path).read().split('\n')


def print_ranks(genes, atlas_ranked_genes, sample_ranked_genes):
  print '\nABSOLUTE RANK'
  for gene in genes:
    try:
      print '%s: sample %d, atlas %d' % (
          gene, sample_ranked_genes.index(gene), atlas_ranked_genes.index(gene))
    except ValueError:
      print '%s: not found' % gene

  print '\nRANK FRACTION'
  for gene in genes:
    try:
      sample_fraction = float(sample_ranked_genes.index(gene))/len(sample_ranked_genes)
      atlas_fraction = float(atlas_ranked_genes.index(gene))/len(atlas_ranked_genes)
      print '%s: sample %f, atlas %f' % (gene, sample_fraction, atlas_fraction)
    except ValueError:
      print '%s: not found' % gene

  print '\nTotal distinct genes: %d sample, %d atlas' % (
      len(sample_ranked_genes), len(atlas_ranked_genes))


if __name__ == "__main__":
  sample = sample_ranks()
  atlas = atlas_ranks()
  genes = read_gene_list()
  print_ranks(genes, atlas, sample)
