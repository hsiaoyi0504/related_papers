import pickle
from collections import Counter
from Bio import Entrez, Medline
import numpy as np
from utils import abstract2words


Entrez.email = 'hsiaoyi0504@gmail.com'


def get_recent_articles():
    search_handler = Entrez.esearch(
        db='pubmed', term='Circadian Rhythm[ALL] OR Circadian Clocks[ALL]', reldate=60, retmax=100000)
    record = Entrez.read(search_handler)
    return record['IdList']


def filter_articles(pubmed_ids, mesh_terms, threshold):
    filtered_ids = []
    clf, word_dict = pickle.load(open('model.pkl', 'rb'))
    for pubmed_id in pubmed_ids:
        fetch_handler = Entrez.efetch(
            db='pubmed', rettype='medline', retmode='text', id=pubmed_id)
        record = Medline.read(fetch_handler)
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
            values = [v / total_words for v in values]
            values = np.array(values).reshape(1, -1)
            # TODO: set a cutoff to ensure the predicted value is useful
            if clf.predict(values) == 1:
                filtered_ids.append(pubmed_id)
        else:
            # TODO: handle articles without abstract
            pass
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
    pubmed_ids = get_recent_articles()
    pubmed_ids = filter_articles(pubmed_ids, MeSH_TERMS, THRESHOLD)
    print(pubmed_ids)
