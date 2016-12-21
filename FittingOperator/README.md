Uses Python 3.5

# Fitting Operator
Script to carry out reasoning given an input file of definite clauses. Assumes a closed world via negation as failure.

	python fittingoperator.py -h
	usage: fittingoperator.py [-h] inputfile

	Takes an input file of definite clauses (with exactly 1 head atom)
	Outputs derived literals with negation as failure

	Definite Clause format:
		[HEAD_ATOM [BODY_POS_ATOMS] [BODY_NEG_ATOMS]]
	Example clause:
		If it is summer, sunny, and not raining, I will go to the beach.
		beach <- summer & sunny & ~rainy
		[beach [summer sunny] [rainy]]

	positional arguments:
	  inputfile   File containing a set of clauses

	optional arguments:
	  -h, --help  show this help message and exit


# Example

Consider the following clauses:

* If someone eats meat, then they are a carnivore.
* If someone eats vegetables and does not eat meat, then they are a herbivore.
* There is an individual named Sam.
* Sam eats ham.
* If someone eats ham, then they eat meat.

The input file would appear as below:

	[Carnivore [EatsMeat] []]
	[Herbivore [EatsVegs] [EatsMeat]]
	[Sam [] []]
	[EatsHam [Sam] []]
	[EatsMeat [EatsHam] []]
	
An example execution:

	Reading clauses ...
		Sam <- 
		Carnivore <- EatsMeat
		Herbivore <- EatsVegs , ~EatsMeat
		EatsHam <- Sam
		EatsMeat <- EatsHam

	atom(s) that finitely fail:
		EatsVegs

	Running Fitting Operator ...
	C = { ~EatsVegs }

	Select clause: Sam <- 
	> clause has no body (is a fact)
	> derived: Sam
	C = { Sam , ~EatsVegs }

	Select clause: EatsHam <- Sam
	> all atoms in the body are derived
	> derived: EatsHam
	C = { EatsHam , Sam , ~EatsVegs }

	Select clause: EatsMeat <- EatsHam
	> all atoms in the body are derived
	> derived: EatsMeat
	C = { EatsMeat , EatsHam , Sam , ~EatsVegs }

	All clauses with same head:
		> Herbivore <- EatsVegs ^ ~EatsMeat
	> all clauses finitely fails
	> derived: ~Herbivore
	C = { EatsMeat , ~Herbivore , EatsHam , Sam , ~EatsVegs }

	Select clause: Carnivore <- EatsMeat
	> all atoms in the body are derived
	> derived: Carnivore
	C = { EatsMeat , EatsHam , Sam , ~EatsVegs , Carnivore , ~Herbivore }
	
It is derived that Sam is a carnivore.