from nltk import word_tokenize
from nltk.corpus import stopwords
from nltk.stem.snowball import SnowballStemmer


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
