import http
from urllib.error import HTTPError
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
    num_total_records = {} # number of records by year
    num_mesh_records = {} # number of records with MeSH term by year
    for i in range(int(start_date), int(end_date) + 1):
        durations[str(i)] = []
        num_total_records[str(i)] = 0
        num_mesh_records[str(i)] = 0
    total_miss_crdt = 0
    for start in range(0, count, batch_size):
        attempt = 1
        while attempt < 3:
            try:
                # maximum retmax of efetch is 10000
                fetch_handler = Entrez.efetch(
                        db='pubmed', rettype='medline', retmode='text',
                        retstart=start, retmax=10000,
                        webenv=webenv, query_key=query_key)
                break
            except HTTPError as err:
                print("Received error from server %s" % err)
                print("Attempt %i of 3" % attempt)
                if attempt == 3:
                    raise
                else:
                    attempt += 1
        total = 0
        while True:
            try:
                record = Medline.read(fetch_handler)
                if record.get('CRDT') is not None:
                    crdt = datetime.strptime(record.get('CRDT')[0], '%Y/%m/%d %H:%M')
                    num_total_records[crdt.strftime('%Y')] += 1
                else:
                    total_miss_crdt += 1
                if record.get('MH') is not None:
                    # although every record have MeSH date, many of them are not labelled
                    edat = datetime.strptime(record.get('EDAT'), '%Y/%m/%d %H:%M')
                    mhda = datetime.strptime(record.get('MHDA'), '%Y/%m/%d %H:%M')
                    durations[crdt.strftime('%Y')].append((mhda - edat).days)
                    num_mesh_records[crdt.strftime('%Y')] += 1
                total += 1
            except http.client.IncompleteRead:
                print('Error: IncompleteRead, following is the detail of the error')
                print(record, start + total)
                break
            except StopIteration:
                print('Finished {}/{} of records'.format(start + total, count))
                break
    print('---------- summary --------------------------------')
    print('Average annotation time by year (day)')
    for i in range(int(start_date), int(end_date) + 1):
        if len(durations[str(i)]) != 0:
            print(str(i), sum(durations[str(i)])/len(durations[str(i)]))
    print('Total number of records by year:')
    for i in range(int(start_date), int(end_date) + 1):
        print(str(i), num_total_records[str(i)])
    print('Total number of records without creating date:', total_miss_crdt)
    print('Total number of records with Mesh terms by year')
    for i in range(int(start_date), int(end_date) + 1):
        print(str(i), num_mesh_records[str(i)])
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
