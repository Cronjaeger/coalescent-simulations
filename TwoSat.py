import networkx as nx
import unittest

# We start out by some unit tests for the Two-sat Solver
class test_two_sat_solver(unittest.TestCase):

    def test_no_solution(self):
        '''case: (a or a) and (not a or not a) [no solution]'''
        clauses = [( bool_var('A',True), bool_var('A',True) ), (bool_var('A',False), bool_var('A',False))]
        sol = TwoSatSolver(clauses)
        self.assertFalse( sol.solvable , msg = 'problem not marked as unsolvable')
        self.assertTrue(sol.solution is None, msg = 'unsolvable problem assigned solution')

    def test_no_solution_two_vars(self):
        '''Unsolvable, 2 variables
Try solving the following 2-SAT problem: (A or B) and (A or not B) and (not
A or B) and (not A or not B). This problem is unsatisfiable. Hence the solver
should return no solution'''
        clauses = [(bool_var('A',x), bool_var('B',y)) for x in [True, False] for y in [True, False]]
        sol = TwoSatSolver(clauses)
        self.assertFalse( sol.solvable , msg = 'problem not marked as unsolvable')
        self.assertTrue(sol.solution is None, msg = 'unsolvable problem assigned solution')

    def test_tautology(self):
        '''The problem "A or not A" should always admit a solution'''
        clauses = [(bool_var('A',True), bool_var('A',False))]
        sol = TwoSatSolver(clauses)
        self.assertTrue( sol.solvable )
        self.assertTrue( (bool_var('A',True) in sol.solution) or (bool_var('A',False) in sol.solution) )
        self.assertFalse( (bool_var('A',True) in sol.solution) and (bool_var('A',False) in sol.solution), msg = 'invalid solution')

    def test_A_or_A(self):
        '''The problem 'A or A' should have the solution A=True'''
        clauses = [(bool_var('A',True), bool_var('A',True))]
        sol = TwoSatSolver(clauses)
        target_solution = [bool_var('A',True)]
        self.assertTrue( sol.solvable )
        self.assertTrue( target_solution == sol.solution, msg = "solution should be %s; instead it is %s."%(target_solution,sol.solution) )

class bool_var(object):
    '''A class to represent boolean variables for use when solving 2-SAT.
bool_var(i,True) = "x_i"; bool_var(i,False) = "not x_i";
i should be a hashable type
    '''
    def __init__(self,name,value = True):
        #assert type(value) is bool
        self.name = name
        self.value = bool(value)

    def __hash__(self):
        return hash((self.name,self.value))

    def __str__(self):
        return 'x_%s'%str(self.name) if self.value else '-x_%s'%str(self.name)


    def __repr__(self):
        #return 'boolean variable(%s)'%str(self)
        return str(self)

    def __eq__(self,other):
        if not isinstance(other, bool_var):
            return False
        else:
            return self.name == other.name and self.value == other.value

    def __cmp__(self,other):
        if not isinstance(other,bool_var):
            return cmp(self.name,other)
        else:
            return cmp(self.name,other.name)

    def negated(self):
        return bool_var(self.name, not self.value)

    def getName(self):
        return self.name

    def index(self):
        '''alias for getName'''
        return self.getName()

    def as_bool(self):
        return self.value

def neg(bool_var):
    return bool_var.negated()

def name_val_pair_to_bool_var(t):

    assert len(t) == 2

    if type(t[0]) is bool and type(t[1]) is not bool:
        t = (t[1],t[0])

    return bool_var(t[0],t[1])

class TwoSatSolver(object):
    '''Implements 2-SAT algorithm from (Aspvall, Plass, and Tarjan 1979)
    '''

    def __init__(self,clauses):
        '''Solve 2-SAT given clauses in conjunctive normal form.
  The algorithm used is an implementation of the one proposed in [1].
  Clauses should be a list of pairs of bool_var-instances. Each pair (A, B) is
interpreted as the statement "A or B".
  A litteral (like 'A') can be either a bool_var, or a pair of the form
(name,sign) (where 'A' would be encoded '('A',True)', and 'not A' would be
('A',False)).
  Given input [(A1, B1), (A2,B2), ... ] we solve the associated 2-SAT problem
which in conjunctive normal form is expressed: (A1 or B1) and (A2 or B2) and (A3
or B3) ...

References:
[1] Bengt Aspvall, Michael F. Plass, Robert Endre Tarjan, A linear-time
algorithm for testing the truth of certain quantified boolean formulas,
Information Processing Letters, Volume 8, Issue 3, 1979, Pages 121-123, ISSN
0020-0190, http://dx.doi.org/10.1016/0020-0190(79)90002-4. '''
        if any(map(lambda clause: len(clause) != 2, clauses)):
            raise ValueError('clauses must be a list of pairs')

        if all(map(lambda clause: isinstance(clause[0],bool_var) and isinstance(clause[1],bool_var), clauses)):
            self.clauses = list(clauses)
        else:
            #raise TypeError('All clauses must be pairs of instances of bool_var')
            self.clauses = [map(name_val_pair_to_bool_var, clause) for clause in clauses]

        self.implications = flatten(map(clause_to_edge,self.clauses))

        # self.variable_names = set([])
        # for x,y in self.implications:
        #     self.variable_names.add(x.name)
        #     self.variable_names.add(y.name)

        self.variable_names = [e[i].name for e in self.implications for i in [0,1]]
        self.variable_names = list(set(self.variable_names)) #remove double-entries
        self.variable_names.sort() # sort the list of variable names

        '''D is an implication graph. each implication of the form A => B is
        represented as a directed edge from A to B'''
        self.D = nx.DiGraph()
        self.D.add_nodes_from([bool_var(i,val) for i in self.variable_names for val in [True, False]])
        self.D.add_edges_from(self.implications)
        try:
            assert len(self.D.nodes()) == 2*len(self.variable_names)
        except AssertionError:
            raise AssertionError('|D| != 2|var_names|\nD.nodes = %s\nvar_names = %s'%(str(self.D.nodes()),str(self.variable_names)))

        ## the strongly connected components of D
        #self.scc = map(list,nx.components.strongly_connected_components(self.D))

        self.solvable = True
        self.solution = None

        C = nx.condensation(self.D)
        D_node_to_C_node = dict(C.graph['mapping'])
        C_node_to_D_nodes = dict([(v,C.node[v]['members']) for v in C.nodes()])

        marked_C_components = []
        partial_solution = []
        for v in nx.topological_sort(C, reverse = True):
            '''
            Whenever an unmarked component is encountered, we set all
            literals in that component to true (which implies that all
            literals in the dual component must be false). The problem
            admits no solution if and only if some component is self-dual.
            if we encounder such a component, we mark the entire problem as
            unsolveable.

            If a component is already marked, we need to do noting furhter.
            The same applies if we already know the problem admits no solution.
            '''
            if (v not in marked_C_components) and self.solvable:

                affected_vars = list(C_node_to_D_nodes[v])
                v_dual = D_node_to_C_node[neg(affected_vars[0])]
                assert all([ v_dual == D_node_to_C_node[neg(x)] for x in affected_vars[1:]])

                if v == v_dual:
                    self.solvable = False
                else:
                    partial_solution.extend(affected_vars)
                    marked_C_components.extend([v,v_dual])

        if self.solvable:
            self.solution = partial_solution
            self.solution.sort()
            assert set(self.variable_names) == set([x.name for x in self.solution])
            assert len(self.variable_names) == len(self.solution)

        assert self.solvable is False or self.solution is not None

        # for c in self.scc:
        #     c_sorted = list(c)
        #     c_sorted.sort()
        #     '''Variables are now sorted by label; It is therefore sufficient to
        #     verify that c_sorted[i] != c_sorted[i+1] holds for all indices. This
        #     holding would imply that a variable is true iff its negation is
        #     true. No such cases should occur.
        #     '''
        #     for i in range(len(c) - 1):
        #         if c_sorted[i] == c_sorted[i+1]:
        #             self.sollvable = False
        #
        # if self.solvable is not False:
        #     self.solvable = True
        #
        # if self.solvable is True:
        #     solution_indices = []
        #     C = nx.condensation(self.D)
        #
        #     for v in nx.topological_sort(C):
        #         literals = C.node[v]['members']
        #         if set([x.index() for x in literals]).isdisjoint(set(solution_indices)):
        #             for x in literals:
        #                 self.solution.append(x)
        #                 #self.solution.append(x if x.value else neg(x))
        #                 solution_indices.append(x.index())
        #
        #     self.solution.sort()
        #     assert set(x.index() for x in self.solution) == set([x.index() for x in self.D.nodes()])

    def __str__(self):
        return '2-SAT problem\n clauses: %s\n solution: %s'%(str(self.clauses),'none exists' if self.solution is None else str(self.solution) )

    def asDict(self):
        if not self.solvable:
            return dict()
        return dict([(x.name,x.value) for x in self.solution])

    def solution_trues(self):
        if not self.solvable:
            return []
        return [x.name for x in self.solution if x.value is True]

    def solution_false(self):
        if not self.solvable:
            return []
        return [x.name for x in self.solution if x.value is False]

def flatten(LOL):
    '''Will lake LOL (a List Of Lists)
    = [[x_1, x_2, ... , x_k], [y_1, y_2, ..., y_r], ... ]
 and return
    [x_1,x_2, ... , x_k, y_1, y_2, ... , y_r, ...]'''
    return reduce(lambda x,y: x+y, LOL, [])

def clause_to_edge(clause):
    assert len(clause) == 2
    assert type(clause[0]) == type(clause[1]) == bool_var
    x,y = clause
    return [(neg(x),y) , (neg(y),x)]
