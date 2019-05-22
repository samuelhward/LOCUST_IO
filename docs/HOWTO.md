# HOW TO

This is a collection of tutorials - enjoy!


Table of Contents
-----------------

* [Combining beam depositions](#Combine-input-particle-lists)


## Combining beam depositions

```python
import context
from classes.input_classes.beam_deposition import Beam_Deposition as bd 

bd_combined=bd(ID='total combined beam depo') #create an empty beam deposition

for beam_depo in beam_depos_to_combine: #if we have lots of other beam depositions stored in a list 'beam_depos_to_combine'
    bd_combined.combine(beam_depo)    
```  