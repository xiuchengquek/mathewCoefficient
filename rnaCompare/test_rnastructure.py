import unittest

from RnaStructure import RnaStructure, StructureScore


class testRnaStructure(unittest.TestCase):
    def setUp(self):
        self.id = 'seqA'
        self.structure = '...(....)...'
        self.double_structure = '...(....)......(....)...'

    def test_rna_structure(self):
        rna_structure_a = RnaStructure(self.id, self.structure)

        self.assertEqual(rna_structure_a.id, self.id)
        self.assertEqual(rna_structure_a.structure, self.structure)

    def test_is_valid(self):
        rna_structure_a = RnaStructure(self.id, self.structure)
        self.assertTrue(rna_structure_a.is_valid())
        rna_structure_a.structure = '...(....(...)...'
        self.assertFalse(rna_structure_a.is_valid())

    def test_get_bp(self):
        rna_structure_a = RnaStructure(self.id, self.structure)
        bp_position = rna_structure_a.get_bp_position()
        self.assertCountEqual(bp_position, [(3,8)])

    def test_get_bp_double(self):
        rna_structure_a = RnaStructure(self.id, self.double_structure)
        bp_position = rna_structure_a.get_bp_position()
        self.assertCountEqual(bp_position, [(3,8), (15,20)])





class testTruePositive(unittest.TestCase):
    def setUp(self):
        self.id = 'seqA'
        self.structure = '...(....)...'
        self.test_structure  = RnaStructure(self.id, '..' + self.structure)
        self.reference_structure = RnaStructure(self.id, self.structure)

    def test_excpetion_raised(self):
        structure_score = StructureScore()
        self.assertRaises(TypeError, lambda : structure_score.add_structure_pair('lol', self.reference_structure)  )
        structure_score.add_structure_pair(self.reference_structure, self.reference_structure)

    def test_get_bp_position(self):
        structure_score = StructureScore()
        structure_score.add_structure_pair(self.reference_structure, self.reference_structure)
        bp = structure_score.get_bp_position()
        self.assertCountEqual(bp, [[(3,8)], [(3,8)]])


    def test_get_tp_1(self):
        structure_score = StructureScore()
        structure_score.add_structure_pair(self.reference_structure, self.reference_structure)
        bp = structure_score.get_bp_position()
        bp_score = structure_score.remove_and_get_tp(bp)
        self.assertEqual(1, bp_score)

    def test_get_tp_0(self):
        structure_score = StructureScore()
        structure_score.add_structure_pair(self.reference_structure, self.test_structure)
        bp = structure_score.get_bp_position()
        bp_score = structure_score.remove_and_get_tp(bp)
        self.assertEqual(0, bp_score)
        self.assertCountEqual([[(3,8)],[(5,10)]], bp)


class testFalseNegative(unittest.TestCase):

    def setUp(self):
        self.reference_structure = RnaStructure('ref' , '..((....))...') # (3,10) , (4,9)
         ## Missing a rna structure of the same bp
        self.one_false_negative = RnaStructure('one_false', '..(..(.).)...') # (3,10), (6,8)
        ## test missing a complate rna strucuture
        self.one_false_negative_b = RnaStructure('one_false', '..(......)...') # (3,10),
        ## Test no missing but test sequence has a additional structure
        self.no_false_negative  = RnaStructure('no_false','..((....))...((.....)).')
        ## Test missing both structure
        self.two_false_negative  = RnaStructure('no_false','((....))((.....)).')


    def test_one_false_negative(self):
        one_false = StructureScore()
        one_false.add_structure_pair(self.reference_structure,self.one_false_negative)
        bp = one_false.get_bp_position()
        fn = one_false.get_fn(bp)
        self.assertEqual(fn, 1)

    def test_one_false_negative_b(self):
        one_false = StructureScore()
        one_false.add_structure_pair(self.reference_structure,self.one_false_negative_b)
        bp = one_false.get_bp_position()
        fn = one_false.get_fn(bp)
        self.assertEqual(fn, 1)

    def test_no_false_negative(self):
        no_false = StructureScore()
        no_false.add_structure_pair(self.reference_structure,self.no_false_negative)
        bp = no_false.get_bp_position()
        fn = no_false.get_fn(bp)
        self.assertEqual(fn, 0)

    def test_two_false_negative(self):
        two_false = StructureScore()
        two_false.add_structure_pair(self.reference_structure, self.two_false_negative)
        bp = two_false.get_bp_position()
        fn = two_false.get_fn(bp)
        self.assertEqual(fn, 2)



class testFalsePositive(unittest.TestCase):


    def setUp(self):
        self.reference_structure = RnaStructure('ref', '..((.....).)..') # (3,12) , (4,10)
        self.two_false_positive = RnaStructure('two_false_positive', '..(.(...))....') # (3,10) , (5,9)

    def test_two_false_positive(self):
        two_false_positive = StructureScore()
        two_false_positive.add_structure_pair(self.reference_structure, self.two_false_positive)
        bp = two_false_positive.get_bp_position()
        fp = two_false_positive.get_fp(bp)
        self.assertEqual(fp, 2)




class testInconsistentPairs(unittest.TestCase):
    """


    Inconsistnet is when one the false negative is in the referennce:
        - for predicted  i-j, reference is j-k theres i-k , and/or h-j
        - where h != i , j != k
    contradicting is for  k-l in refernece, i-j in sequence, when k < i < l < j
    """

    def setUp(self):
        self.reference_structure = RnaStructure('ref', '..((.....).)..') # (3,12) , (4,10)

        # (3,10 and 3,12)
        self.one_inconsistent = RnaStructure('one_inconsistent', '..(......)....') # (3,10)
        #  (4,12 and 3,12)
        self.one_inconsistent_b = RnaStructure('one_inconsistent', '...(.......)...') # (4,12)
        # (2,12) and (5,10)
        self.two_inconsistent = RnaStructure('two_inconsistent' , '.(..(....).)..') # (2,12) , (5,10) -
        # No inconsistent, by 2 fp
        self.no_inconsistent_two_fp = RnaStructure('no_inconsistent_2fp' , '....(..)........(...)')  # (5,8) ,  (17,21)

    def test_one_inconsistent(self):
        ss = StructureScore()
        ss.add_structure_pair(self.reference_structure, self. one_inconsistent)
        bp = ss.get_bp_position()
        one_inconsistent = ss.get_inconsistent_pair(bp)
        self.assertCountEqual(one_inconsistent, [(2,9)])

    def test_one_inconsistent_b(self):
        ss = StructureScore()
        ss.add_structure_pair(self.reference_structure, self. one_inconsistent_b)
        bp = ss.get_bp_position()
        one_inconsistent = ss.get_inconsistent_pair(bp)
        self.assertCountEqual(one_inconsistent, [(3,11)])

    def test_two_inconsistent(self):
        ss = StructureScore()
        ss.add_structure_pair(self.reference_structure, self.two_inconsistent)
        bp = ss.get_bp_position()
        two_inconsistent = ss.get_inconsistent_pair(bp)
        self.assertCountEqual(two_inconsistent, [(1,11) , (4,9)])

    def test_no_inconsistent_two_fp(self):
        ss =  StructureScore()
        ss.add_structure_pair(self.reference_structure, self.no_inconsistent_two_fp)
        bp = ss.get_bp_position()
        fp = ss.get_fp(bp)
        self.assertEqual(fp, 2)
        no_inconsistent = ss.get_inconsistent_pair(bp)
        self.assertCountEqual(no_inconsistent, [])

class testContradictingPairs(unittest.TestCase):
    def setUp(self):
        self.reference_structure = RnaStructure('ref', '..((.....).)..') # (3,12) , (4,10)
        self.one_contradicting = RnaStructure('one_contradicting', '....(..........).') # (5,16)
        self.two_contradicting = RnaStructure('two_contradicting', '....(.(...)....).') # (5,16) (7,11)
        self.no_contradicting = RnaStructure('one_inconsistent', '..(......)....')


    def test_one_contradicting(self):
        ss = StructureScore()
        ss.add_structure_pair(self.reference_structure, self.one_contradicting)
        bp = ss.get_bp_position()
        one_contradicting = ss.get_contradicting_pairs(bp)
        self.assertListEqual(one_contradicting, [(4,15)])


    def test_two_contradicting(self):
        ss = StructureScore()
        ss.add_structure_pair(self.reference_structure, self.two_contradicting)
        bp = ss.get_bp_position()
        two_contradicting = ss.get_contradicting_pairs(bp)
        self.assertCountEqual(two_contradicting, [(4,15), (6,10)])

    def test_no_contradicting(self):
        ss = StructureScore()
        ss.add_structure_pair(self.reference_structure, self.no_contradicting)
        bp =  ss.get_bp_position()
        no_contradicting = ss.get_contradicting_pairs(bp)
        self.assertCountEqual(no_contradicting, [])

class testGetFPCorrected(unittest.TestCase):

    def setUp(self):
        self.reference_structure = RnaStructure('ref', '..((.....).)..') # (3,12) , (4,10)
        self.one_contradicting = RnaStructure('one_contradicting', '....(..........).') # (5,16)
        self.one_inconsistent = RnaStructure('one_inconsistent', '..(......)....') # (3,10)
        self.inconsistent_with_contradicting = RnaStructure('inconsistent_contradicting' , '..(.(....)....).' )

    def test_one_contradicting(self):
        ss = StructureScore()
        ss.add_structure_pair(self.reference_structure, self.one_contradicting)
        bp = ss.get_bp_position()
        one_contradicting = ss.get_fp_corrected(bp)
        self.assertEqual(one_contradicting, 1)

    def test_one_inconsistent(self):
        ss = StructureScore()
        ss.add_structure_pair(self.reference_structure, self. one_inconsistent)
        bp = ss.get_bp_position()
        one_inconsistent = ss.get_fp_corrected(bp)
        self.assertEqual(one_inconsistent, 1)

    def test_inconsistent_with_contradicting(self):
        ss = StructureScore()
        ss.add_structure_pair(self.reference_structure, self. inconsistent_with_contradicting)
        bp = ss.get_bp_position()
        inconsistent_with_contradicting = ss.get_fp_corrected(bp)
        self.assertEqual(inconsistent_with_contradicting, 2)

class testCalculateScore(unittest.TestCase):
    def setUp(self):
        self.reference_structure = RnaStructure('ref', '(...)..(.((.....).).)........') # (3,12) , (4,10)
        ## Thre will be 2tp, 2fp, 2fn,                  (..(.(....).))'
        self.predicted_structure = RnaStructure('inconsistent_contradicting' , '(...)..(..(.(....).))...(...)')

    def test_calculate_score(self):
        ss = StructureScore()
        ss.add_structure_pair(self.reference_structure, self.predicted_structure)
        score = ss.calculate_score()
        self.assertDictEqual(score,
            {
                'fp' : 2,
                'tp' : 2,
                'fn' : 2
            }
        )















































































































if __name__ == '__main__':
    unittest.main()


