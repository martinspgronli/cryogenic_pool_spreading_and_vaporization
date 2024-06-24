This repository contains Python code for simulating the spill of liquid ammonia or liquid hydrogen onto solid ground
with variable thermal properties. Both boiling correlations and perfect thermal contact have been implemented. 
The code uses the Python interface of [Clawpack](https://www.clawpack.org), called PyClaw, to solve the shallow water equations. 
Several parameters can be varied including initial spill velocity, ground topography, obstructions, and details 
regarding the thermal properties of the substrate. More details can be found in: 
[Spill, evaporation and substrate thermal transport model for liquid H2 and NH3 (2024)](ADD_LINK_TO_PUBLICATION.com).

After cloning, pip install requirements by typing 
```bash
pip install -r /path/to/requirements.txt
```
in terminal.

Choose fluid in `main.py`. For NH3, adjust parameters in `variables_NH3.py`. For H2, adjust parameters in `variables_H2.py`.

Run `main.py` to simulate the spill and generate plots.