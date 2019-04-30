from datetime import datetime
from Bio import Entrez, Medline
from utils import get_recent_articles

Entrez.email = 'hsiaoyi0504@gmail.com'


def measure_annotation_time(start_date, end_date, count, webenv, query_key, mesh_terms):
    filtered_ids = []
    batch_size = 10000
    # debug code (modify count)
    # count = 10001
    durations = {}  # durations by year
    for i in range(int(start_date), int(end_date) + 1):
        durations[str(i)] = []
    for start in range(0, count, batch_size):
        # maximum retmax of efetch is 10000
        fetch_handler = Entrez.efetch(
                db='pubmed', rettype='medline', retmode='text',
                retstart=start, retmax=10000,
                webenv=webenv, query_key=query_key)
        records = Medline.parse(fetch_handler)
        for record in records:
            if record.get('MH') is not None:
                # although every record have MeSH date, many of them are not labelled
                crdt = datetime.strptime(record.get('CRDT')[0], '%Y/%m/%d %H:%M')
                edat = datetime.strptime(record.get('EDAT'), '%Y/%m/%d %H:%M')
                mhda = datetime.strptime(record.get('MHDA'), '%Y/%m/%d %H:%M')
                durations[crdt.strftime('%Y')].append((mhda - edat).days)
    for i in range(int(start_date), int(end_date) + 1):
        if len(durations[str(i)]) != 0:
            print(str(i), sum(durations[str(i)])/len(durations[str(i)]))
    return filtered_ids


if __name__ == '__main__':
    start_date = '2015'
    end_date = '2018'
    MeSH_TERMS = [
        'Algorithms', 'Animals', 'CLOCK Proteins/metabolism',
        'Circadian Rhythm/*physiology', 'Circadian Clocks/physiology*',
        'Databases, Genetic', 'Gene Expression Profiling/*methods', 'Humans',
        'Liver/metabolism/physiology', 'Liver Neoplasms/metabolism',
        'Lung/metabolism/physiology', 'Machine Learning', 'Mice',
        'Statistics as Topic/*methods', 'Transcription, Genetic/genetics'
    ]
    count, webenv, query_key = get_recent_articles(start_date, end_date)
    pubmed_ids = measure_annotation_time(start_date, end_date, count, webenv, query_key, MeSH_TERMS)
