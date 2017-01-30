#!/usr/bin/python
import unittest
from decode_util import *

class DecodeTest(unittest.TestCase):
    def test_stack(self):
        s = Stack()
        s.push(state("I am", -0.1))
        s.push(
        
