# -*- coding: utf-8 -*-
""" A class to manage conversion of SequenceVariants to repeat representation"""

from dataclasses import dataclass
import re
import hgvs
from hgvs.pretty.models import VariantCoords


@dataclass(eq=True, repr=True, frozen=True, order=True)
class RepeatUnit:
    repeat_count: int
    repeat_unit:str
    block_size:int
    block: str


class RepeatAnalyser:

    def __init__(self, fs: VariantCoords, reverse:bool=False) -> None:

        self.is_repeat = False
        self.repeat_units_ref = None
        self.repeat_units_alt = None
        self.ref_str = fs.ref
        self.alt_str = fs.alt
        self.reverse = reverse

        # protect from too large sequences getting analyzed.
        if len(self.ref_str) > hgvs.global_config.repeats.max_repeat_length:
            self.repeat_units_ref = []
            self.repeat_units_alt = []
            return

        self.repeat_units_ref = detect_repetitive_block_lengths(self.ref_str, reverse=self.reverse)
        self.repeat_units_alt = detect_repetitive_block_lengths(self.alt_str, reverse = self.reverse)

        if len(self.repeat_units_ref) == 0 and len(self.repeat_units_alt) ==0 :
            return
                
        # check longest repeat blocks:
        # we only look at ref to determine if there are repeats 
        # If ref has no repeat, we don't call this a repeat variant, even if alt would have a repetitive unit
        longest_r_unit = self._get_longest_repeat_unit(self.repeat_units_ref)
        if longest_r_unit is None:
            return
            
        # filter our too fragmented results
        expected_size = len(self.ref_str) / 3
        if longest_r_unit.block_size < expected_size:            
            return
            
        if longest_r_unit.repeat_unit not in self.alt_str:
            return

        self.repeat_units_alt = detect_repetitive_block_lengths(self.alt_str, longest_ref_unit = longest_r_unit, reverse = self.reverse)

        self.is_repeat = True

        ref_repeat = get_repeat_str(self.ref_str, self.alt_str, self.repeat_units_ref, self.repeat_units_alt, self.reverse)
        alt_repeat = get_repeat_str(self.alt_str, self.ref_str, self.repeat_units_alt, self.repeat_units_ref, self.reverse)
        self.ref_str = ref_repeat
        self.alt_str = alt_repeat

    def __repr__(self):
        return f"{self.ref_str}>{self.alt_str}"

    def _get_longest_repeat_unit(self, repeat_units:list[RepeatUnit])->RepeatUnit:
        lru = None
        for ru in repeat_units:
            if not lru:
                lru = ru
                continue

            if ru.block_size > len(lru.block):
                lru = ru

        return lru


def detect_repetitive_block_lengths(sequence: str, longest_ref_unit:RepeatUnit|None=None, reverse: bool = False) -> list[RepeatUnit]:
    """Detects the length of repetitive blocks in a string, with an option to search from left to right or reverse.

    In reverse mode, it creates the largest possible blocks of the smallest possible units.
    """
    result: list[RepeatUnit] = []
    seq_len = len(sequence)

    if longest_ref_unit is not None:
        # look for full containment of the longest ref repeat
        # this is so we can detect [2]>[1] (or really any repeat length variation to just 1)
        pattern = f'({re.escape(longest_ref_unit.repeat_unit)})+'
        match = re.search(pattern, sequence)
        
        if match:
            repeat_count = len(sequence) // len(longest_ref_unit.repeat_unit)
            block = repeat_count* longest_ref_unit.repeat_unit
            ru = RepeatUnit(repeat_count, longest_ref_unit.repeat_unit, len(block), block)                        
            result.append(ru)
            shuffleable_bases = sequence[repeat_count*len(ru.repeat_unit):]
            rus = detect_repetitive_block_lengths(shuffleable_bases, reverse = reverse)
            result.extend(rus)
            return result
            

        
    if reverse:
        i = seq_len  # Start from the end of the sequence
        while i > 0:
            matched = False
            # Iterate over block sizes from smallest to largest
            for block_size in range(1, seq_len // 2 + 1):
                if i - block_size < 0:
                    continue  # Not enough characters to form a repeat unit

                # Extract the potential repeat unit ending at position i
                repeat_unit = sequence[i - block_size:i]

                # Calculate the maximum possible number of repeats for this block size
                max_possible_repeats = i // block_size

                # Initialize repeat_count to the maximum possible and decrease
                for repeat_count in range(max_possible_repeats, 1, -1):
                    start_index = i - repeat_count * block_size
                    if start_index < 0:
                        continue  # Not enough characters to form the repetitive block

                    # Extract the substring that could be the repetitive block
                    substr = sequence[start_index:i]

                    # Build the regex pattern for the current repeat unit and count
                    pattern = rf'({re.escape(repeat_unit)})' + r'{' + f'{repeat_count}' + r'}'

                    # Check if the substring matches the pattern
                    if re.fullmatch(pattern, substr):
                        repetitive_block = substr
                        ru = RepeatUnit(
                            block=repetitive_block,
                            block_size=len(repetitive_block),
                            repeat_unit=repeat_unit,
                            repeat_count=repeat_count
                        )
                        result.append(ru)
                        # Move the index `i` backward by the length of the repetitive block
                        i -= len(repetitive_block)
                        matched = True
                        break  # Found the largest block for this unit size

                if matched:
                    break  # Proceed to the next position `i`

            if not matched:
                # No repeat found, remove one character
                ru = RepeatUnit(
                    block=sequence[i - 1],
                    block_size=1,
                    repeat_unit=sequence[i - 1],
                    repeat_count=1
                )
                result.append(ru)
                i -= 1  # Move back by one character if no match is found

    else:
        i = 0  # Start from the beginning of the sequence
        while i < seq_len:
            matched = False
            for block_size in range(1, seq_len // 2 + 1):
                # Build the regex pattern for the current block size
                pattern = rf'(.{{{block_size}}})\1+'
                match = re.match(pattern, sequence[i:])

                if match:
                    repetitive_block = match.group()  # The full repeating pattern
                    repeated_unit = match.group(1)    # The repeating unit
                    repetition_count = len(repetitive_block) // len(repeated_unit)

                    # Add the repetitive block and its details to the result
                    ru = RepeatUnit(
                        block=repetitive_block,
                        block_size=len(repetitive_block),
                        repeat_unit=repeated_unit,
                        repeat_count=repetition_count
                    )
                    result.append(ru)

                    # Move the index `i` forward by the length of the repetitive block
                    i += len(repetitive_block)
                    matched = True
                    break  # Found a match, break the loop for current block size

            if not matched:
                i += 1  # Move forward by one character if no match is found

    return result


def get_repeat_str(
    sequence: str,
    other_seq: str,
    primary_repeat_unit: list[RepeatUnit],
    other_repeat_unit: list[RepeatUnit],
    reverse: bool = False
) -> str:
    if len(primary_repeat_unit) == 0 and len(other_repeat_unit) == 0:
        return None

    if len(primary_repeat_unit) == 0 and len(other_repeat_unit) == 1 and sequence == other_repeat_unit[0].repeat_unit:
        return f"{sequence}[1]"
    elif len(primary_repeat_unit) == 0 and len(other_repeat_unit) == 1 and sequence != other_repeat_unit[0].repeat_unit:
        return None

    if len(primary_repeat_unit) > 0 and len(other_repeat_unit) > 0:
        return_str = assemble_repeat_string(sequence, primary_repeat_unit, reverse=reverse)
        return return_str

    if len(other_repeat_unit) == 0 and len(other_seq) > 0:
        return_str = assemble_repeat_string(sequence, primary_repeat_unit, reverse=reverse)
        if len(return_str) > 0:
            return return_str

    return None

def assemble_repeat_string(sequence: str, repeat_units: list[RepeatUnit], reverse: bool = False) -> str:
    return_str = ""
    primary_repeat_unit = repeat_units.copy()
    seq = sequence
    if reverse:
        while len(seq) > 0:
            found_unit = None
            for ru in primary_repeat_unit:
                if seq.endswith(ru.block):
                    return_str = f"{ru.repeat_unit}[{ru.repeat_count}]" + return_str
                    seq = seq[:-len(ru.block)]
                    found_unit = ru
                    break
            if not found_unit:
                # remove one character if no matching repeat unit is found
                count = 1
                seq_char = seq[0]
                seq = seq[:-1]

                # count consecutive repeating chars  from the end
                while seq and seq[-1] == seq_char:
                    count += 1
                    seq = seq[:-1]

                return_str = f"{seq_char}[{count}]" + return_str

    else: # forward direction
        while len(seq) > 0:
            found_unit = None
            for ru in primary_repeat_unit:
                if seq.startswith(ru.block):
                    return_str += f"{ru.repeat_unit}[{ru.repeat_count}]"
                    seq = seq[len(ru.block):]
                    found_unit = ru
                    break
            if not found_unit:
                # remove one character if no matching repeat unit is found
                count = 1
                seq_char = seq[0]
                seq = seq[1:]

                # count consecutive repeating chars 
                while seq and seq[0] == seq_char:
                    count += 1
                    seq = seq[1:]
                
                return_str += f"{seq_char}[{count}]"

    
    return return_str

