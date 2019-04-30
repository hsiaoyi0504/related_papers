import pickle
from collections import Counter
from Bio import Entrez, Medline
import numpy as np
from utils import abstract2words, get_recent_articles


Entrez.email = 'hsiaoyi0504@gmail.com'


def filter_articles(count, webenv, query_key, mesh_terms, threshold):
    filtered_ids = []
    clf, word_dict = pickle.load(open('model.pkl', 'rb'))
    batch_size = 10000
    # TODO: batch prediction
    for start in range(0, count, batch_size):
        # maximum retmax of efetch is 10000
        fetch_handler = Entrez.efetch(
                db='pubmed', rettype='medline', retmode='text',
                retstart=start, retmax=10000,
                webenv=webenv, query_key=query_key)
        records = Medline.parse(fetch_handler)
        for record in records:
            pubmed_id = record['PMID']
            record_mesh_terms = record.get('MH')
            if record_mesh_terms is not None:
                count = 0
                for mesh in mesh_terms:
                    if mesh in record_mesh_terms:
                        count += 1
                if count >= threshold:
                    filtered_ids.append(pubmed_id)
                    continue
            abstract = record.get('AB')
            if abstract is not None:
                words = abstract2words(abstract)
                c = Counter(words)
                values = []
                for i in range(len(word_dict)):
                    values.append(c[word_dict[i]])
                total_words = sum(values)
                if total_words == 0:
                    continue
                values = [v / total_words for v in values]
                values = np.array(values).reshape(1, -1)
                # TODO: set a cutoff to ensure the predicted value is useful
                if clf.predict(values) == 1:
                    filtered_ids.append(pubmed_id)
            # Don't need to consider article without abstract
    return filtered_ids


if __name__ == '__main__':
    MeSH_TERMS = [
        'Algorithms', 'Animals', 'CLOCK Proteins/metabolism',
        'Circadian Rhythm/*physiology', 'Circadian Clocks/physiology*',
        'Databases, Genetic', 'Gene Expression Profiling/*methods', 'Humans',
        'Liver/metabolism/physiology', 'Liver Neoplasms/metabolism',
        'Lung/metabolism/physiology', 'Machine Learning', 'Mice',
        'Statistics as Topic/*methods', 'Transcription, Genetic/genetics'
    ]
    THRESHOLD = 3
    start_date = '2015'
    end_date = '2018'
    count, webenv, query_key = get_recent_articles(start_date, end_date)
    pubmed_ids = filter_articles(count, webenv, query_key, MeSH_TERMS, THRESHOLD)
    print(pubmed_ids)
