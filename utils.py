from math import ceil
from nltk import word_tokenize
from nltk.corpus import stopwords
from nltk.stem.snowball import SnowballStemmer
from Bio import Entrez


def get_recent_articles(start_date, end_date):
    # TODO: handle case with over 100000 result
    search_handler = Entrez.esearch(
        db='pubmed', term='{}:{}[CRDT]'.format(start_date, end_date), retmax=100000, usehistory='y')
    record = Entrez.read(search_handler)
    id_list = record['IdList']
    web_env = record['WebEnv']
    query_key = record['QueryKey']
    if int(record['Count']) > 100000:
        times = ceil(int(record['Count']) / 100000) - 1
        start = 100000
        for i in range(times):
            print(i)
            start += 100000
            search_handler = Entrez.esearch(
                db='pubmed', term='{}:{}[CRDT]'.format(start_date, end_date), retmax=100000, usehistory='y', webenv=web_env, query_key=query_key)
            id_list += record['IdList']
            web_env = record['WebEnv']
            query_key = record['QueryKey']
    print(record['Count'])
    print(id_list)
    return (id_list, web_env, query_key)


def abstract2words(abstract):
    all_stopwords = stopwords.words('english')
    stemmer = SnowballStemmer('english')
    words = word_tokenize(abstract)
    # Remove single-character tokens (mostly punctuation)
    words = [word for word in words if len(word) > 1]
    # Remove numbers
    words = [word for word in words if not word.isnumeric()]
    # Stemming
    words = [stemmer.stem(word) for word in words]
    # Remove stopwords
    words = [word for word in words if word not in all_stopwords]
    return words
