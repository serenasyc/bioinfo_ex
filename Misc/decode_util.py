import heapq
import math
from collections import namedtuple
import numpy as np

scores = namedtuple("scores", "lm, reorder, phr_fe, lex_fe, phr_ef, lex_ef")
phrase = namedtuple("phrase", "english, phr_fe, lex_fe, phr_ef, lex_ef")
phraseloc = namedtuple("phraseloc", "s, t, phrase")
state = namedtuple('state', 'bigram, b, r, scores, futurecost')
parameters = namedtuple("parameters", "b, d, e, nbest, s, v")
candidate = namedtuple('candidate','sentence, features, bleu')

## CLASS

class Stack:
  '''
  actually a max priority queue
  '''
  def __init__(self):
      self.queue = []

  def push(self, state):
      score = state.futurecost + state.scores.lm
      heapq.heappush(self.queue, (-score, state))

  def pop(self):
      return heapq.heappop(self.queue)[1]

  def __len__(self):
      return len(self.queue)

## FUNCTIONS

def get_phrases(q, f, tm, d=6):
  '''
  ref: ph(q) in Collins
  q   TUPLE   state/partial hypothesis=( bigram, bitstring, endindex, score)
  f   TUPLE   french sentence split by words
  tm  DICT    translation model (phrase table)
  d   INT     distortion limit
  Given a bigram (in q), find valid phrases that extend the partial hypothesis
  valid if:
    1. phrase does not overlap with translations in parital hypothesis
    2. distortion limit is not violated
  '''
  phrases = []

  n = len(f)
  bitstring = bitmap2str(q.b, n)
  pbit = prefix1bits(q.b)
  maxindex = min(n+1, pbit+d+1)
  for i in range(pbit, maxindex):
    for j in range(i+1, maxindex):
      if bitstring[i:j] != "."*(j-i):
        break
      else:
        if f[i:j] in tm:
          for phrase_lexicon in tm[f[i:j]]:
            # phrase_lexicon = (english, logprobs)
            p = phraseloc(i, j, phrase_lexicon)
            phrases.append(p)
  return phrases


def next_hypothesis(p, q, lm, cost, eta):
  '''
  ref: next(q,p) from Collins
  p     TUPLE   phrase=(startindex, endindex, english, logprob)
  q     TUPLE   state/partial hypothesis=( bigram, bitstring, endindex, score)
  lm    CLASS   language model; given bigram and word, return language model score
  cost  DICT    dictionary of future cost estimates
  eta   FLOAT   factor to penalize distortion
  Given a partial hypothesis and a phrase to extend it with,
  return a new partial hypothesis
  '''
  lm_score = 0.0
  new_bigram = q.bigram
  for e in p.phrase.english.split():
    (new_bigram, partial_score) = lm.score(new_bigram, e)
    lm_score += partial_score
  new_bitstring = q.b + bitmap(range(p.s, p.t))
  new_endindex = p.t-1
  qscores = np.array(q.scores)
  pscores = np.array([lm_score,
                      -eta*float(abs(q.r-p.s)),
                      p.phrase.phr_fe,
                      p.phrase.lex_fe,
                      p.phrase.phr_ef,
                      p.phrase.lex_ef])
  new_score = scores(*(qscores+pscores))
  new_futurecost = q.futurecost - cost[p.s, p.t]
  new_state = state(new_bigram, new_bitstring, new_endindex,
                    new_score, new_futurecost)
  return new_state

def eq(q1,q2):
  '''
  q1  TUPLE   state/partial hypothesis=( bigram, bitstring, endindex, score)
  q2  TUPLE   another state/partial hypothesis
  test if the two partial hypotheses are equal
  '''
  if (q1.bigram == q2.bigram
      and q1.b == q2.b
      and q1.r == q2.r):
    return True
  else:
    return False

def combine_scores(scores, futurecost, weights):
  ''' 
  scores  NAMEDTUPLE  lm, reorder, phr_fe, lex_fe, phr_ef, lex_ef
  futurecost  FLOAT   future cost estimate
  uniform weights to combine log probabilities with future cost estimate
  '''
  weighted_score = np.array(scores)*weights
  return sum(weighted_score, futurecost)

def beam(stack, beta, stack_size, weights):
  '''
  stack LIST    of states with number of words translated
  beta  FLOAT
  Only keep hypotheses that are a max of beta difference from the maxscore
  '''
  Q = Stack()
  best_state = stack.pop()
  maxscore = combine_scores(best_state.scores, best_state.futurecost, weights)
  Q.push(best_state)

  while len(stack) != 0 and len(Q) < stack_size:
    state = stack.pop()
    if combine_scores(state.scores, state.futurecost, weights) >= (maxscore - beta):
      Q.push(state)
    else:
      break
  return Q

def Add(stack, q1, q0, p, backpointers):
  '''
  stack LIST  of states with number of words translated
  q1    TUPLE state
  '''
  in_stack = False
  for q2 in stack.queue:
    q2 = q2[1]
    if eq(q2, q1):
      in_stack = True
      if q2.scores.lm < q1.scores.lm:
        backpointers[q1] = (q0, p.phrase.english)
        stack.queue.append((-(q1.scores.lm + q1.futurecost), q1))
        stack.queue.remove((-(q2.scores.lm + q2.futurecost), q2))
      else:
        return stack
  if not in_stack:
    backpointers[q1] = (q0, p.phrase.english)
    stack.queue.append((-(q1.scores.lm + q1.futurecost), q1))
  heapq.heapify(stack.queue)
  return stack

def backtrack(lastq, backpointers):
  english = ''
  english_words = []
  scores = lastq.scores
  while lastq.bigram != ('<s>',):
    (lastq,english) = backpointers[lastq]
    english_words.append(english)
  sentence =  ' '.join(english_words[::-1])
  return sentence,scores

##borrowed from CMPT 413 @ SFU (Fall 2016) taught by Anoop Sarkar
def bitmap(sequence):
  """ Generate a coverage bitmap for a sequence of indexes """
  return reduce(lambda x,y: x|y, map(lambda i: long('1'+'0'*i,2), sequence), 0)

def bitmap2str(b, n, on='o', off='.'):
  """ Generate a length-n string representation of bitmap b """
  return '' if n==0 else (on if b&1==1 else off) + bitmap2str(b>>1, n-1, on, off)

def prefix1bits(b):
  """ Count number of bits encountered before first 0 """
  return 0 if b&1==0 else 1+prefix1bits(b>>1)
