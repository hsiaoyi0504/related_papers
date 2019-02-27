from Bio import Entrez, Medline
import csv

Entrez.email = 'hsiaoyi0504@gmail.com'


def collect_example():
    for date in ['2015', '2016', '2017', '2018']:
        f = open('./example_{}.csv'.format(date), 'w')
        writer = csv.writer(f)
        writer.writerow(['PMID', 'title', 'abstract', 'label'])
        search_handler = Entrez.esearch(
            db='pubmed', term='CLOCK Proteins/metabolism[MESH] OR Circadian Rhythm/*physiology[MESH] OR Circadian Clocks/physiology*[MESH]', mindate='{}/01/01'.format(date), maxdate='{}/12/31'.format(date), retmax=100000, usehistory='y')
        record = Entrez.read(search_handler)
        start = 0
        fetch_handler = Entrez.efetch(
                    db='pubmed', rettype='medline', retmode='text',
                    retstart=start, retmax=10000,
                    webenv=record['WebEnv'], query_key=record['QueryKey'])
        records = Medline.parse(fetch_handler)
        for record in records:
            pmid = record.get('PMID')
            title = record.get('TI')
            abstract = record.get('AB')
            if pmid is not None and title is not None and abstract is not None:
                writer.writerow([pmid, title, abstract, ''])


if __name__ == '__main__':
    collect_example()
