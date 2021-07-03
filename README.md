# Supplementary code for nonparametric inference for interventional effects with multiple mediators

[![MIT
license](http://img.shields.io/badge/license-MIT-brightgreen.svg)](http://opensource.org/licenses/MIT)

**Authors:** [David
Benkeser](https://davidbphd.com), Jialu Ran

-----

## Description

The `intermed_sims` repository contains the code needed to reproduce the simulation results of the manuscript, "Nonparametric inference for interventional effects with multiple mediators."

-----

## The `intermed` package

The simulations rely on the [`intermed`](https://github.com/benkeser/intermed) R package, which can be downloaded and installed from GitHub as follows.

```r
devtools::install_github("benkeser/intermed")
```

-----

## Additional simulation 3

Code for the additional simulation 3 in the appendix is included in the `add_sim3` directory, which includes its own `README` file detailing the contents. This simulation study compares the proposed approach with parametric approaches by [`Vansteelandt and Daniel (2017)`](https://pubmed.ncbi.nlm.nih.gov/27922534/).

-----

## Issues

If you encounter any bugs, please [file an issue](https://github.com/benkeser/intermed_sims/issues).

-----

## License

© 2021- [David Benkeser](https://davidbphd.com)

The contents of this repository are distributed under the MIT license.
See below for details:

    The MIT License (MIT)
    
    Copyright (c) 2021- David C. Benkeser
    
    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:
    
    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.
    
    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.

