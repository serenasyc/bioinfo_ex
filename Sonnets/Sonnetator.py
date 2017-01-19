#!/usr/bin/python
import nltk
from collections import defaultdict
import math
import random

def get_rhymes(tokens):
  '''
  tokens  LIST
  Returns 
  '''
  d = nltk.corpus.cmudict.dict()
  rhymes = defaultdict(list)
  allwords = list(set(tokens))
  
  for word in allwords:
    if word in d:
      rhymes[d[word][0][-1]].append(word)
    # if word not in cmudict, subword might be
    else:
      for i in range(len(word)):
        subword = word[i:]
        if subword in d:
          x = d[subword][0][-1]
          rhymes[x].append(word)
          d[word] = [[x]]
          break
  return rhymes, d


def random_ngram(fdist,last_word):
    l = [x for x in fdist if x[-1] == last_word]
    return random.choice(l)
  

def next_ngram(ngram, fdist):
  '''
  ngram   tuple
  fdist   nltk.probability.FreqDist of ngrams
  Returns an negram starting with the last (n-1) tokens of ngram
  '''
  n1gram = ngram[1:]
  l = list(filter(lambda x : x[0][:-1] == n1gram, fdist.items()))

  if ngram[0] == "s/":
      return random.choice(l)[0]
  
  maxfreq = max(map(lambda x: x[1], l))
  maxfreq_grams = list(map(lambda y: y[0], 
                             filter(lambda x: x[1] > (maxfreq-5), l)))
  next_g = random.choice(maxfreq_grams)
  # avoid the black hole of infinite loops!
  if next_g[1] == n1gram[0]:
      return random.choice(l)[0]
  else:
      return next_g


def prev_ngram(ngram, fdist):
  '''
  ngram   TUPLE
  fdist   nltk.probability.FreqDist of ngrams
  Returns an negram ending with the first (n-1) tokens of ngram
  '''
  n1gram = ngram[:-1]
  l = list(filter(lambda x : x[0][1:] == n1gram, fdist.items()))
  
  if ngram[-1] == "/s":
      return random.choice(l)[0]
  
  maxfreq = max(map(lambda x: x[1], l))
  maxfreq_grams = list(map(lambda y: y[0], 
                             filter(lambda x: x[1] > (maxfreq-5), l)))
  
  prev_g = random.choice(maxfreq_grams)
  # avoid the black hole of infinite loops!
  if prev_g[0] == ngram[1]:
      return random.choice(l)[0]
  else:
      return prev_g


# ----- MAIN -----

# read Sonnets, each line starts with "s/ s/" and ends with "/s"
f = open('ShakespearesSonnets3.txt')
raw = f.read().lower()
f.close()

tokens = nltk.word_tokenize(raw)
trigrams = nltk.trigrams(tokens)
fdist3 = nltk.FreqDist(trigrams)

rhymes,d = get_rhymes(tokens)

sonnet, sentence = [], []
i = 0
while i < 7: 
    ngram = ("s/","s/","s/")    

    while "/s" not in ngram:
        ngram = next_ngram(ngram, fdist3)
        sentence.append(ngram[-1])
    if len(sentence) > 5:
        sentence[0] = sentence[0][0].upper() + sentence[0][1:]
        sonnet.append((" ".join(sentence[:-1])))
        i += 1
        
        last_index = -1
        while sentence[last_index] in [",",".","!","?",";",":","'","/s"]:
            last_index -= 1
        last_word = sentence[last_index]
        
        syllable = d[last_word][0][-1]
        new_last_word = random.choice(rhymes[syllable])
        ngram = random_ngram(fdist3, new_last_word)
        sentence = [new_last_word]
        
        while True:
            while "s/" not in ngram:
                ngram = prev_ngram(ngram, fdist3)
                sentence.append(ngram[0])
            if len(sentence) > 5:
                sentence = sentence[::-1]
                sentence[1] = sentence[1][0].upper() + sentence[1][1:]
                sonnet.append(" ".join(sentence[1:]))
                sentence = []
                break
            else:
                sentence.pop()
                ngram = random_ngram(fdist3,"/s")
    else:
        sentence.pop()

order = [0, 2, 1, 3, 4, 6, 5, 7, 8, 10, 9, 11, 12, 13]

temp = []
print()
for i in order:
    temp.append(sonnet[i])
    
for i in range(4):
    print(temp[i])
print()
for i in range(4,8):
    print(temp[i])
print()
for i in range(8,12):
    print(temp[i])
print()
for i in range(12,14):
    print(temp[i])
print()
