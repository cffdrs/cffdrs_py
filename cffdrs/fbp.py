"""
Fire Behavior Prediction System function

Calculates the outputs from the Canadian Forest Fire Behavior Prediction (FBP) System.

:param input: List of FBPInput or single FBPInput
:param output: "Primary", "Secondary", or "All"
:param m: Not used in Python version
:param cores: Not used in Python version

:returns: List of FBPPrimaryOutput, FBPSecondaryOutput, or FBPAllOutput
"""

from typing import List, Union
from cffdrs.models import FBPAllOutput, FBPInput, FBPPrimaryOutput, FBPSecondaryOutput
from cffdrs.fire_behaviour_prediction import fire_behaviour_prediction


def fbp(input: Union[FBPInput, List[FBPInput]] = None, output: str = "Primary") -> List[Union[FBPPrimaryOutput, FBPSecondaryOutput, FBPAllOutput]]:
    if input is None:
        input = FBPInput()
    if isinstance(input, FBPInput):
        input = [input]
    results = []
    for inp in input:
        result = fire_behaviour_prediction(inp, output)
        results.append(result)
    return results