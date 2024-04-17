import copy
import random
from dataclasses import dataclass
from typing import List, Dict, Tuple
import numpy as np
from debugpy.common.compat import izip


@dataclass
class Swap:
    old_cell: int
    new_cell: int

@dataclass
class SingleCellSwap:
    events_amnt: int
    swap: Swap

@dataclass
class HotSpot:
    experiment_type: str
    experiment_name: str
    id: int
    hotspot_cells: List[int]
    best_swaps: List[Swap]

@dataclass
class SwapPermutation:
    swaps: List[Swap]


EXPERIMENTS_HOT_SPOTS = [

    HotSpot(
        experiment_type='mid_third_wild_type',
        experiment_name='sample_9',
        hotspot_cells=[9,13,18,20,21,24,32],
        id=0,
        best_swaps=[
            Swap(20, 5),
            Swap(24, 11),
            Swap(9, 0),
            Swap(32, 4),
            Swap(18, 12),
            Swap(13, 13),
            Swap(21, 1),
        ]
    ),
    HotSpot(
        experiment_type='mid_third_wild_type',
        experiment_name='sample_9',
        hotspot_cells=[2,4,5,30],
        id=1,
        best_swaps=[]
    ),

    HotSpot(
        experiment_type='mid_third_wild_type',
        experiment_name='sample_1',
        hotspot_cells=[11,12,17, 13,22,25,26],
        id=0,
        best_swaps=[]
    ),

    HotSpot(
        experiment_type='mid_third_wild_type',
        experiment_name='sample_7',
        hotspot_cells=[8,9,11,19,23,25,26,27,28,30,33,35,39,41],
        id=0,
        best_swaps=[]
    ),

    HotSpot(
        experiment_type='mid_third_wild_type',
        experiment_name='sample_2',
        hotspot_cells=[8,9,10,42],
        id=0,
        best_swaps=[
            Swap(9, 15),
            Swap(8, 36),
            Swap(42, 30),
            Swap(10, 27),
        ]
    ),

    HotSpot(
        experiment_type='mid_third_wild_type',
        experiment_name='sample_4',
        hotspot_cells=[19,23,24,28,40,41],
        id=0,
        best_swaps=[]
    ),

    HotSpot(
        experiment_type='mid_third_wild_type',
        experiment_name='sample_8',
        hotspot_cells=[29, 37, 41, 10, 34, 49],
        id=1,
        best_swaps=[
            Swap(29, 16),
            Swap(37, 28),
            Swap(41, 48),
            Swap(10, 2),
            Swap(34, 44),
            Swap(49, 6),
        ]
    ),

    HotSpot(
        experiment_type='mid_third_wild_type',
        experiment_name='sample_11',
        hotspot_cells=[7, 32, 36, 3, 6, 5],
        id=0,
        best_swaps=[
            Swap(3, 42),
            Swap(36, 16),
            Swap(32, 34),
            Swap(6, 55),
            Swap(7, 30),
        ]
    ),

    HotSpot(
        experiment_type='mid_third_wild_type',
        experiment_name='sample_10',
        hotspot_cells=[51, 99, 39, 63, 53],
        id=0,
        best_swaps=[
            Swap(51, 20),
            Swap(99, 22),
            Swap(39, 54),
        ]
    ),

    HotSpot(
        experiment_type='mid_third_wild_type',
        experiment_name='sample_10',
        hotspot_cells=[26, 29, 80, 90, 96],
        id=1,
        best_swaps=[]
    ),

    HotSpot(
        experiment_type='mid_third_wild_type',
        experiment_name='sample_10',
        hotspot_cells=[3,36,18,35,78,77,16],
        id=2,
        best_swaps=[
            Swap(3,82),
            Swap(36,32),
            Swap(18,57),
            Swap(35,7),
            Swap(78,102),
            Swap(77,76),
            Swap(16,27),
        ]
    ),

    HotSpot(
        experiment_type='mid_third_wild_type',
        experiment_name='sample_6',
        hotspot_cells=[30,33,43,44,48,84],
        id=1,
        best_swaps=[]
    ),

    HotSpot(
        experiment_type='mid_third_wild_type',
        experiment_name='sample_6',
        hotspot_cells=[0,1,4,5,  15,17,18, 22,23,25,26,  62,82,83,99],
        id=2,
        best_swaps=[]
    ),
    HotSpot(
        experiment_type='mid_third_wild_type',
        experiment_name='sample_6',
        hotspot_cells=[31,9,55,24,45,46,60,37,34,61,72,79,92,104],
        id=3,
        best_swaps=[]
    ),


    ####################
    ####### ZPG ########
    ####################
    HotSpot(
        experiment_type='mid_third_zpg_RNAi',
        experiment_name='sample_5',
        hotspot_cells=[1,2,5,8],
        id=0,
        best_swaps=[]
    ),

    HotSpot(
        experiment_type='mid_third_zpg_RNAi',
        experiment_name='sample_1',
        hotspot_cells=[0,4,5,6,11,15],
        id=0,
        best_swaps=[]
    ),
    HotSpot(
        experiment_type='mid_third_zpg_RNAi',
        experiment_name='sample_1',
        hotspot_cells=[8,9,13],
        id=1,
        best_swaps=[]
    ),


    HotSpot(
        experiment_type='mid_third_zpg_RNAi',
        experiment_name='sample_4',
        hotspot_cells=[1,12,17,19,23],
        id=0,
        best_swaps=[]
    ),


    HotSpot(
        experiment_type='mid_third_zpg_RNAi',
        experiment_name='sample_3',
        hotspot_cells=[2,8,13,17,19,20,24],
        id=0,
        best_swaps=[]
    ),



    HotSpot(
        experiment_type='mid_third_zpg_RNAi',
        experiment_name='sample_2',
        hotspot_cells=[14,36,37],
        id=0,
        best_swaps=[]
    ),

    HotSpot(
        experiment_type='mid_third_zpg_RNAi',
        experiment_name='sample_2',
        hotspot_cells=[9,26,30],
        id=1,
        best_swaps=[]
    ),


    HotSpot(
        experiment_type='mid_third_zpg_RNAi',
        experiment_name='sample_2',
        hotspot_cells=[0,5,11,24,27],
        id=2,
        best_swaps=[]
    ),




    ##########################
    ####### CBX 3.125 ########
    ##########################
    HotSpot(
        experiment_type='cbx_inhibitor_3.125',
        experiment_name='sample_2',
        hotspot_cells=[1,14,51],
        id=0,
        best_swaps=[]
    ),
    HotSpot(
        experiment_type='cbx_inhibitor_3.125',
        experiment_name='sample_2',
        hotspot_cells=[34,19,29,48,35,9,42,36,45,18,47],
        id=1,
        best_swaps=[]
    ),

    HotSpot(
        experiment_type='cbx_inhibitor_3.125',
        experiment_name='sample_1',
        hotspot_cells=[7,8,19],
        id=0,
        best_swaps=[]
    ),


    ############################
    ####### CBX washout ########
    ############################
    HotSpot(
        experiment_type='cbx_inhibitor_washout',
        experiment_name='sample_4',
        hotspot_cells=[2, 7, 11, 31],
        id=0,
        best_swaps=[]
    ),

    HotSpot(
        experiment_type='cbx_inhibitor_washout',
        experiment_name='sample_4',
        hotspot_cells=[10, 44, 14, 43],
        id=1,
        best_swaps=[]
    ),


    ############################
    ####### WT late 2nd ########
    ############################
    HotSpot(
        experiment_type='late_second_wild_type',
        experiment_name='sample_3',
        hotspot_cells=[2, 4, 6, 8],
        id=0,
        best_swaps=[]
    ),
    HotSpot(
        experiment_type='late_second_wild_type',
        experiment_name='sample_3',
        hotspot_cells=[3, 6, 7, 10],
        id=1,
        best_swaps=[]
    ),

    HotSpot(
        experiment_type='late_second_wild_type',
        experiment_name='sample_4',
        hotspot_cells=[0,1,2,3,4,6,8,9,11,14,16,18],
        id=0,
        best_swaps=[]
    ),

    HotSpot(
        experiment_type='late_second_wild_type',
        experiment_name='sample_1',
        hotspot_cells=[2,13,23,11,7,21,6],
        id=0,
        best_swaps=[]
    ),

    #############################
    ####### WT early 3rd ########
    #############################
    HotSpot(
        experiment_type='early_third_wild_type',
        experiment_name='sample_1',
        hotspot_cells=[67,52,22,11,2,48,20,3,8,6,21,7, 63,41],
        id=0,
        best_swaps=[]
    ),

    HotSpot(
        experiment_type='early_third_wild_type',
        experiment_name='sample_1',
        hotspot_cells=[16,17,24,9, 27, 64],
        id=1,
        best_swaps=[]
    ),

    HotSpot(
        experiment_type='early_third_wild_type',
        experiment_name='sample_3',
        hotspot_cells=[10,62,5,20,61],
        id=0,
        best_swaps=[]
    ),

    HotSpot(
        experiment_type='early_third_wild_type',
        experiment_name='sample_3',
        hotspot_cells=[0,9,12,60,81],
        id=1,
        best_swaps=[]
    ),

    HotSpot(
        experiment_type='early_third_wild_type',
        experiment_name='sample_3',
        hotspot_cells=[2,57,67,69],
        id=2,
        best_swaps=[]
    ),

    HotSpot(
        experiment_type='early_third_wild_type',
        experiment_name='sample_3',
        hotspot_cells=[3,24,31,40,45,46,47,53],
        id=3,
        best_swaps=[]
    ),

    #############################
    ####### ZPG late 2nd ########
    #############################


    ##############################
    ####### ZPG early 3rd ########
    ##############################

    HotSpot(
        experiment_type='early_third_zpg_RNAi',
        experiment_name='sample_1',
        hotspot_cells=[7,10,15],
        id=0,
        best_swaps=[]
    ),

]
@dataclass
class Candidate:
    cell_id: int
    delta_activations: float

def _find_cell_with_same_activations(original_cell_ix: int, activations: List[int], old_cells: List[int]) -> List[Candidate]:
    """
    Finding cells with SAME or GREATER activations as <original_cell_ix> (up to 120%).
    Candidate can be of the original hub.
    :param original_cell_ix:
    :param activations:
    :return:
    """
    res = []
    required_activation = activations[original_cell_ix]

    for ix, activation in enumerate(activations):
        if activation >= required_activation:
            res.append(Candidate(cell_id=ix, delta_activations=(activation-required_activation)/float(required_activation)))
    return res


def _rec_perm(l_of_l, n, current_perm=[]):
    if len(l_of_l) == 0:
        return [current_perm]

    else:
        l_of_l_copy = copy.deepcopy(l_of_l)
        curr_list = l_of_l_copy.pop(0)
        perms = []

        for i in curr_list:
            if i not in current_perm:
                tmp_perms = _rec_perm(l_of_l_copy, n, current_perm + [i])
                if len(tmp_perms) > 0:
                    perms.extend(tmp_perms)
                    n -= len(tmp_perms)
                if n <= 0:
                    return perms
        return perms

def _generate_permutations(possible_candidates, old_cells):
    """ Generate up to n permutations of possible candidates.
    Allow up to <same_hub_percentage> of the old hub to be in the new hub."""
    [random.shuffle(pc) for pc in possible_candidates]

    n = 10000
    same_hub_percentage = 0.5

    perms = _rec_perm(possible_candidates, n)

    filtered_perms = list(filter(lambda permutation: len(set(permutation)) == len(permutation)
                                and len(set(permutation).intersection(set(old_cells))) <= same_hub_percentage*len(old_cells), perms))

    if len(filtered_perms) > n:
        return random.sample(filtered_perms, n)
    else:
        return filtered_perms



def _check_permutation_mean_delta_in_acts(permutation, old_cells, candidates_per_cell):
    """ Calculate the mean delta in activations of a permutation. """
    perm_deltas = []
    for ix, cand_cell_id in enumerate(permutation):
        original_cell_id = old_cells[ix]
        delta = list(filter(lambda cand: cand.cell_id == cand_cell_id, candidates_per_cell[original_cell_id]))[
            0].delta_activations
        perm_deltas.append(delta)
    return np.mean(perm_deltas)


def _check_permutations_mean_delta_in_acts(permutations, old_cells, candidates_per_cell):
    """ Calculate the mean delta in activations of a list of permutations."""
    deltas = []
    for permutation in permutations:
        deltas.append(_check_permutation_mean_delta_in_acts(permutation, old_cells, candidates_per_cell))
    return np.mean(deltas)


def create_swap_candidates(old_cells: List[int], activations: List[int], n_samples: int) -> Tuple[List[SwapPermutation], List[float]]:
    """ Create a list of swap candidates for the hotspots. The overall delta in activations of each permutations is limited to 100%-120%."""

    def _get_all_swaps(activation_candidates, old_cells, candidates_per_cell,
                       existing_swaps, existing_permutations, existing_delta, n):
        candidates = [[cand.cell_id for cand in cands_list] for cands_list in activation_candidates.values()]
        sorted_old_cells = copy.deepcopy(old_cells)
        sorted_lists = sorted(izip(sorted_old_cells, candidates), reverse=False, key=lambda x: len(x[1]))
        sorted_old_cells, candidates = [[x[i] for x in sorted_lists] for i in range(2)]

        valid_permutations = _generate_permutations(candidates, sorted_old_cells)
        unsorted_to_sorted_map = {old_cells[ix]: sorted_old_cells.index(old_cells[ix]) for ix in range(len(old_cells))}
        unsorted_valid_permutation = [[permutation[unsorted_to_sorted_map[cell_id]] for cell_id in old_cells] for permutation in valid_permutations]
        valid_permutations = unsorted_valid_permutation
        valid_permutations = random.sample(valid_permutations, min(n, len(valid_permutations)))

        for permutation in valid_permutations:
            if permutation not in existing_permutations and n > 0:
                swap_permutation = SwapPermutation(
                    swaps=[Swap(old_cell=old_cells[ix], new_cell=permutation[ix])
                           for ix in range(len(permutation))])
                existing_swaps.append(swap_permutation)
                existing_permutations.append(permutation)
                existing_delta.append(_check_permutation_mean_delta_in_acts(permutation, old_cells, candidates_per_cell))
                n -= 1

        return existing_swaps, existing_permutations, existing_delta


    N = n_samples

    # Find candidates for each cell in the hub according to the activation level
    candidates_per_cell: Dict[int, List[Candidate]] = {
        cell_id: _find_cell_with_same_activations(cell_id, activations, old_cells) for cell_id in old_cells
    }

    # Divide candidates to groups based on the delta in activations
    same_activations_candidates = {}
    lte_activations_candidates = {}
    gte_activations_candidates = {}
    for cell_id, candidates in candidates_per_cell.items():
        same_activations_candidates[cell_id] = list(filter(lambda cand: cand.delta_activations == 0, candidates))
        lte_activations_candidates[cell_id] = list(filter(lambda cand: cand.delta_activations <= 0, candidates))
        gte_activations_candidates[cell_id] = list(filter(lambda cand: cand.delta_activations >= 0, candidates))

    existing_swaps, existing_permutations, existing_delta = _get_all_swaps(same_activations_candidates, old_cells, candidates_per_cell,
                                                                           [], [], [], n=N)
    N = N - len(existing_swaps)
    existing_swaps, existing_permutations, existing_delta = _get_all_swaps(gte_activations_candidates, old_cells, candidates_per_cell,
                                                                           existing_swaps, existing_permutations, existing_delta, n=N)
    return existing_swaps, existing_delta
