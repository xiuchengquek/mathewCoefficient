from collections import deque, Counter

#from rnaCompare import rnaStructure


class RnaStructure:
    """
    Class to that encapsulate the a rna seconday structure and its id.
    Contains method to get all bp position from the rna structure

    """
    def __init__(self, id, structure):
        """
        Construct a RnaStructure Instance
        :param id: RNA strucuture id
        :param structure: seconday structure in string , ie '.....(.....)...
        :return:
        """
        self.id = id
        self.structure = structure

    def is_valid(self):
        """
        Method to check that the rna strcuture is correct.The number left and right parenthesis should balance
        example  : '..(....(...)...' returns false
        example : '..(....(...)...). ' returns true
        :return: true or false
        """
        bases = Counter(self.structure)
        return bases['('] == bases[')']

    def get_bp_position(self):
        """
        Find all exisitng base pairs in RNA seconday Structure.
        Return a list of tuples of base pair position
        eg.  '...(....)...' returns [(3,8)]

        :return: list containing tuples of base pair position
        """
        unmatched_left = []
        base_pairs = []
        for ix, bp in enumerate(self.structure):
            if bp == '(':
                unmatched_left.append(ix)
            elif bp == ')':
                current_left = unmatched_left.pop()
                base_pairs.append((current_left, ix))
        return base_pairs




class StructureScore:
    """
    Class that encapsulate 2 rna structure as a list. The first element is considered as the reference
    Contains methods that takes the first element as reference and calculate the following metrics

    TP : true positive
    FN : false negative
    FP : false positive
    """

    def add_structure_pair(self, rna_structure_a, rna_structure_b):
        if not (isinstance(rna_structure_a, RnaStructure) and isinstance(rna_structure_a, RnaStructure)):
            raise TypeError('Need Rnastrucutre class, got %s and %s' % (type(rna_structure_a), type(rna_structure_b)))
        self.rna_structure_pair = (rna_structure_a, rna_structure_b)



    def calculate_score(self):
        """

        find the following scores
        1) true positive
        2) false negative
        3) false positive
            - inconsistent + contradicting
        :return: score false position, true postiion, false negative and stimga
        """

        base_pairs = self.get_bp_position()
        tp = self.remove_and_get_tp(base_pairs)
        fn = self.get_fn(base_pairs)
        fp = self.get_fp_corrected(base_pairs)

        score = {
            'fp' : fp,
            'fn' : fn,
            'tp' : tp,
        }

        return score


    def get_bp_position(self):
        """
        Method to get get the bp position of each rna strcuture

        :return:
        """
        seq_a_bp = self.rna_structure_pair[0].get_bp_position()
        seq_b_bp = self.rna_structure_pair[1].get_bp_position()
        return [seq_a_bp, seq_b_bp]


    def get_fp_corrected(self, base_pairs):
        """
        method to get the number of false positive.
        false positive is define as predicted bp that is absent in the reference
        however not all fp are even : only contracting fp and inconsistent pairs are penalizin

        so fp with penalty is given as contradicting + inconsistent

        :param base_pairs: [bp in reference, bp in predicted ]
        :return:
        """
        bp_a = base_pairs[0]
        bp_b = base_pairs[1]

        contradicting_pairs = self.get_contradicting_pairs(base_pairs)
        inconsistent_pairs = self.get_inconsistent_pair(base_pairs)
        fp = set(contradicting_pairs) | set(inconsistent_pairs)
        return len(fp)


    ## Get true positive and remove them for better comparsion later
    def remove_and_get_tp(self, base_pairs):
        """
        if both base_pairs matches in references and predicted
        true positives in base_pairs are removed.
        ..-(...)...
        ...(...)...
        :param base_pairs: [bp in reference, bp in predicted ]
        :return: the number of true positive
        """
        bp_a  = base_pairs[0]
        bp_b = base_pairs[1]

        true_bp = set(bp_a) & set(bp_b)

        tp_score = len(true_bp)

        # remove tp pairs from eqaution


        return tp_score

    def get_fn(self, base_pairs):
        """
        False negative is defined as bp in reference but not detected in predicted structure
        :param base_pairs: [bp in reference, bp in predicted ]
        :return: the number of fn
        """
        bp_a = base_pairs[0]
        bp_b = base_pairs[1]
        missing_reference = set(bp_a)  - set(bp_b)
        return(len(missing_reference))

    def get_fp(self, base_pairs):
        """
        False positve is bp in predicted but not in reference
        :param base_pairs:  [bp in reference, bp in predicted ]
        :return: number of false positive
        """
        bp_a = base_pairs[0]
        bp_b = base_pairs[1]
        fp = set(bp_b) - set(bp_a)
        return len(fp)


    def get_inconsistent_pair(self, base_pairs):
        """
        Inconsistnet is when one the false negative is in the referennce:
        - for predicted  i-j, reference is j-k theres i-k , and/or h-j
        - where h != i , j != k
        :param base_pairs:  [bp in reference, bp in predicted ]
        :return: [inconsistent bp]
        """
        bp_a = base_pairs[0]
        bp_b = base_pairs[1]


        reference_left = [x[0] for x in bp_a]
        reference_right = [x[1] for x in bp_a]

        fp = set(bp_b) - set(bp_a)

        inconsistent_pair = []

        for x in fp :
            if x[0] in reference_left or x[1] in reference_right:
                inconsistent_pair.append(x)


        return inconsistent_pair

    def get_contradicting_pairs(self, base_pairs):
        """
        The function is used to get contradicting pair in the predicted alignment,
         Where :
            for  k-l in refernece,
            i-j in predicted
         Then
            k < i < l < j

        returns i-j



        :param base_pairs: [bp in reference, bp in predicted ]
        :return: [ contradicting bp ]
        """

        bp_a = base_pairs[0]
        bp_b = base_pairs[1]

        fp  = set(bp_b) - set(bp_a)
        contradicting = []

        for x in fp :
            for y in bp_a:
                if (y[0] < x[0] < y[1] < x[1] ):
                    contradicting.append(x)
        return list(set(contradicting))










