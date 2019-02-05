from collections import Counter
import pickle
from Bio import Entrez, Medline
from sklearn.neighbors import KNeighborsClassifier
from utils import abstract2words

Entrez.email = 'hsiaoyi0504@gmail.com'


def get_abstracts(file_name):
    pubmed_ids = []
    with open(file_name) as f:
        for line in f:
            pubmed_ids.append(int(line.rstrip('\n')))
    abstracts = []
    for pubmed_id in pubmed_ids:
        fetch_handler = Entrez.efetch(
            db='pubmed', rettype='medline', retmode='text', id=str(pubmed_id))
        record = Medline.read(fetch_handler)
        abstracts.append(record['AB'])
    return abstracts


if __name__ == '__main__':
    pos_abstracts = get_abstracts('positive_examples')
    neg_abstracts = get_abstracts('negative_examples')
    all_abstracts = pos_abstracts + neg_abstracts
    all_words = set()
    for i, a in enumerate(all_abstracts):
        words = abstract2words(a)
        all_words.update(words)
        all_abstracts[i] = words
    word_dict = dict(zip(range(len(all_words)), all_words))
    X = []
    for a in all_abstracts:
        c = Counter(a)
        values = []
        for i in range(len(all_words)):
            values.append(c[word_dict[i]])
        total_words = sum(values)
        values = [v / total_words for v in values]
        X.append(values)
    y = [1] * len(pos_abstracts) + [0] * len(neg_abstracts)
    clf = KNeighborsClassifier(n_neighbors=3)
    clf.fit(X, y)
    pickle.dump([clf, word_dict], open('model.pkl', 'wb'))
