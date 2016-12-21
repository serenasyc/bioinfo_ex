#!/bin/usr/env python
import argparse
from collections import defaultdict, namedtuple
from argparse import RawTextHelpFormatter
from random import shuffle

clause = namedtuple("clause", "head, body_pos, body_neg")

def read_inputfile(filename):
    '''
    filename    STR     filename of file containing Definite Clauses
                        format: [[head] [body+] [body-]c]
    returns     DICT    clauses[head] = clause
                LIST    finitely failed atoms
    '''
    all_head_atoms = set()
    all_body_atoms = set()
    # a directory of clauses where keys are head atoms
    # value is a list because multiple clauses may have same head atom
    clause_dir = defaultdict(list)
    
    with open(filename, 'r') as f:
        lines = [line.strip() for line in f.readlines()]
        for line in lines:
            index = 0
            atom = ""
            cl = clause([],[],[])
            d = [cl.head, cl.body_pos, cl.body_neg]
            for i in line[1:-1]:
                if i == "[":
                    index += 1
                    atom = ""
                elif i in [" ", "]"]:
                    if len(atom) == 0:
                        continue
                    else:
                        d[index].append(atom)
                        atom = ""
                else:
                    atom += i
            
            print("{} <- {}".format(cl.head[0], " , ".join(cl.body_pos + [negate(x) for x in cl.body_neg])))
##            # tautologies are not used
##            if (len(cl.body_pos) != 0
##                and len(cl.body_neg) != 0
##                and len(set(cl.body_pos).intersection(set(cl.body_neg))) != 0):
##                continue
            
            clause_dir[cl.head[0]].append(cl)           
            
            all_head_atoms.add(cl.head[0])
            body_atoms = set(cl.body_pos + cl.body_neg)
            all_body_atoms = all_body_atoms.union(body_atoms)
    
    # atoms that finitely fails are those found in the body but not the head
    # of any clauses in KB
    finitely_fails = set(all_body_atoms.difference(all_head_atoms))
    print('atom(s) that finitely fail:')
    print("\t"+" , ".join([r for r in finitely_fails]))

    return clause_dir, set([negate(x) for x in finitely_fails])

def negate(atom):
    '''
    input:  STR    atom   
    output: STR    negated input atoms
    '''
    if atom[0] == "~":
        return atom[1:]
    else:
        return "~"+atom

def check_atoms(atoms, C):
    '''
    input:  LIST    atoms   of atoms
    input:  SET     C       of derived atoms
    output: 1st BOOL    True if atoms or its negation in C
            2nd BOOL    False if negated atoms in C or not in C
    '''
    intersect = set(atoms).intersection(C)

    # atoms in C
    if len(intersect) == len(atoms):
        return True,True
    neg_intersect = set(map(lambda x: negate(x), atoms)).intersection(C)
    # negation of atoms in C
    if len(neg_intersect) == len(atoms):
        return True,False
    # atoms have undetermined values
    return False, False

def split_clause(clause):
    '''
    input:  clause    TUPLE   Definite Clause
    output: head    STR     head atom
    output: body    LIST    of body atoms, with "~" if negative   
    '''
    head = clause.head[0]
    body_pos = clause.body_pos
    body_neg = list(map(negate, clause.body_neg))
    body = body_pos + body_neg

    return head, body

def get_clause(clause):
    AND = " ^ "
    return("{} <- {}{}{}".format(clause[0][0],
                                      AND.join(clause[1]),
                                      AND if (len(clause[2])!= 0 and len(clause[1]) != 0) else "",
                                      AND.join(map(negate,clause[2]))
                                      ))

def get_C(C):
    return("C = { " + " , ".join(list(filter(lambda x: x != None, C))) + " }")
    
def fittingoperator(clause_dir, C, verbose=True):
    '''
    input:
    clause_dir  DICT    clause_dir[head] = clause(head, body_pos, body_neg)
    C           LIST    of consequents (starts with finitely failed atoms)
    output:     SET     of derived atoms

    if verbose = True, prints out atoms as they're derived
    '''
    print(get_C(C))
    clauses = []
    for h in clause_dir:
        for c in clause_dir[h]:
            clauses.append(c)

    success = False
    while True:
        # to randomize selection
        shuffle(clauses)
        
        L = len(C)
        # if an atom is derived, I want to remove the rules with that atom as the head
        heads_to_rm = set([])
        
        for r in clauses:
            # if an atom is derived
            if success:
                print(get_C(C))

            # re-initialize loop
            success = False
            all_clauses = []
            # holds info to print out
            derived_atom = ""
            msg = ""
            
            head,body = split_clause(r)
            if (head not in C) and (negate(head) not in C):
                # if there's no body, the head atom is a fact
                if len(body) == 0:
                    C.add(head)
                    heads_to_rm.add(head)
                    
                    success = True
                    msg = "clause has no body (is a fact)"
                    derived_atom = head
                    
                # if each atom in the body is in C
                elif check_atoms(body,C) == (True,True):
                    C.add(head)
                    heads_to_rm.add(head)

                    success = True
                    msg = "all atoms in the body are derived"
                    derived_atom = head
                    
                # if at least one of the atoms in the body finitely fails
                else:
                    neg_atoms = set(map(negate, body))
                    # an atom should finitely fail ONLY IF you fail to prove
                    # all clauses with that atom as the head
                    if len(C.intersection(neg_atoms)) != 0:
                        complete = True
                        positive = False
                        # check all clauses with same atom as head
                        for clause in clause_dir[head]:
                            rhead, rbody = split_clause(clause)
                            rcomplete, rpositive = check_atoms(rbody, C)
                                                        
                            complete &= rcomplete
                            positive |= rpositive
                            all_clauses.append(clause)
                            if not complete:
                                break
       
                        if complete:
                            success = True
                            if positive:
                                C.add(head)
                                msg = "atom(s) in body of at least one clause is derived"
                                derived_atom = head
                            else:
                                C.add(negate(head))
                                msg = "all clauses finitely fails"
                                derived_atom = negate(head)
                if success:
                    heads_to_rm.add(head)
                    if verbose:
                        if len(all_clauses) != 0:
                            print("\nAll clauses with same head:")
                            print("\t> "+"\n\t> ".join(map(get_clause, all_clauses)))
                        else:
                            print("\nSelect: clause: " + get_clause(r))

                        print("> "+msg)
                        print("> derived: {}".format(derived_atom))                       
        if len(C) == L:
            break
        heads_to_rm = heads_to_rm.union(set(map(negate, heads_to_rm)))
        for r in clauses:
            head,body = split_clause(r)
            if head in heads_to_rm:
                clauses.remove(r)
    return C

def main():
    parser = argparse.ArgumentParser(description="Takes an input file of definite clauses (with exactly 1 head atom)\n"\
                                     "Outputs derived literals with negation as failure\n"\
                                     "\nDefinite Clause format:\n\t[HEAD_ATOM [BODY_POS_ATOMS] [BODY_NEG_ATOMS]]\n"\
                                     "Example clause:\n\tIf it is summer, sunny, and not raining, I will go to the beach.\n"\
                                     "\tbeach <- summer & sunny & ~rainy\n"\
                                     "\t[beach [summer sunny] [rainy]]\n", formatter_class=RawTextHelpFormatter) 
    parser.add_argument("inputfile", type=str, help="File containing a set of clauses")
    args = parser.parse_args()
    inputfile = args.inputfile

    print("Reading clauses ...")
    clause_dir, C = read_inputfile(inputfile)
    
    print("\nRunning Fitting Operator ...")
    C = fittingoperator(clause_dir, C)

if __name__ == '__main__':
    main()
