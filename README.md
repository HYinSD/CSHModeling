# C-S-H Atomic Structure Simulation

## Introduction
This project simulates the atomic structure of C-S-H (Calcium Silicate Hydrate) materials. By adjusting atomic interactions and charge distributions, we can analyze the performance of C-S-H under different calcium-to-silicon ratios(C/S = 1.2–2.3).

## Table of Contents
- [Installation](#installation)
- [Usage](#usage)
- [Features](#features)
- [Contributing](#contributing)
- [License](#license)

## Installation
1. Ensure you have Python 3 and the NumPy library installed.
2. Clone this repository:
   ```bash
   git clone https://github.com/HYinSD/CSHModeling
   cd CSHModeling
3. Install the required dependencies:
   pip install numpy
   pip install tkinter
   pip install random

## Usage
1. Place the orthorhombic 11 Å tobermorite coordinate file inputcoord.txt in the project root directory.
2. Run the simulation script:

   python CSHModeling.py
   
3. The output data will be saved in the output.data file.

## Features
1. Load atomic data and adjust charge distributions.
2. Expand atomic cells to generate a three-dimensional model.
3. Calculate the number of target aggregates based on calcium-to-silicon ratios.

## Contributing
Contributions are welcome! Please raise issues or submit pull requests. Ensure you follow coding standards and guidelines.

## License
(Copyright (C) 2024  Shijie Wang (shijiewang54@gmail.com), Hang Yin (yinh@sdau.edu.cn)

**CSHModeling** is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program; if not, you can obtain it from https://www.gnu.org/licenses/old-licenses/gpl-2.0.html.)
