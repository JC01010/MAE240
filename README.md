<p align="center"><h1 align="center">MAE240: Solar Radiation Pressure (SRP) Effect on L2 Orbital Motion: Application to James Webb Space Telescope (JWST)</h1></p>
<p align="center"><h2 align="center">Taewoo Kang, Jingtong Sheng, Jason Chang, UC San Diego</h2></p>
<p align="center">
    <em><code>❯ REPLACE-ME</code></em>
</p>
<p align="center">
    <img src="https://img.shields.io/github/license/JC01010/MAE240?style=default&logo=opensourceinitiative&logoColor=white&color=0080ff" alt="license">
    <img src="https://img.shields.io/github/last-commit/JC01010/MAE240?style=default&logo=git&logoColor=white&color=0080ff" alt="last-commit">
    <img src="https://img.shields.io/github/languages/top/JC01010/MAE240?style=default&color=0080ff" alt="repo-top-language">
    <img src="https://img.shields.io/github/languages/count/JC01010/MAE240?style=default&color=0080ff" alt="repo-language-count">
</p>
<p align="center"><!-- default option, no dependency badges. -->
</p>
<p align="center">
    <!-- default option, no dependency badges. -->
</p>
<br>

##  Table of Contents

- [ Project Structure](#-project-structure)
  - [ Project Index](#-project-index)
- [ Getting Started](#-getting-started)
  - [ Prerequisites](#-prerequisites)
  - [ Installation](#-installation)
- [ License](#-license)
- [ Acknowledgments](#-acknowledgments)

---


##  Project Structure

```sh
└── MAE240/
    ├── MAE240report.pdf
    ├── README.md
    ├── averageanalysis.m
    └── exactanalysis
        └── MAE 240 Project
```


###  Project Index
<details open>
    <summary><b><code>MAE240/</code></b></summary>
    <details open> <!-- __root__ Submodule -->
        <summary><b>__root__</b></summary>
        <blockquote>
            <table>
            <tr>
                <td><b><a href='https://github.com/JC01010/MAE240/blob/master/averageanalysis.m'>averageanalysis.m</a></b></td>
                <td><code>❯ REPLACE-ME</code></td>
            </tr>
            </table>
        </blockquote>
    </details>
    <details open> <!-- exactanalysis Submodule -->
        <summary><b>exactanalysis</b></summary>
        <blockquote>
            <details open>
                <summary><b>MAE 240 Project</b></summary>
                <blockquote>
                    <table>
                    <tr>
                        <td><b><a href='https://github.com/JC01010/MAE240/blob/master/exactanalysis/MAE 240 Project/inertial_to_syn.m'>inertial_to_syn.m</a></b></td>
                        <td><code>Converts coordinates from the inertial frame to the synodic (rotating) frame.</code></td>
                    </tr>
                    <tr>
                        <td><b><a href='https://github.com/JC01010/MAE240/blob/master/exactanalysis/MAE 240 Project/rgb.m'>rgb.m</a></b></td>
                        <td><code>Utility function for generating RGB color values for plotting.</code></td>
                    </tr>
                    <tr>
                        <td><b><a href='https://github.com/JC01010/MAE240/blob/master/exactanalysis/MAE 240 Project/synodic_eom.m'>synodic_eom.m</a></b></td>
                        <td><code>Defines the equations of motion in the synodic (rotating) frame.</code></td>
                    </tr>
                    <tr>
                        <td><b><a href='https://github.com/JC01010/MAE240/blob/master/exactanalysis/MAE 240 Project/syn_to_inertial.m'>syn_to_inertial.m</a></b></td>
                        <td><code>Converts coordinates from the synodic (rotating) frame to the inertial frame.</code></td>
                    </tr>
                    <tr>
                        <td><b><a href='https://github.com/JC01010/MAE240/blob/master/exactanalysis/MAE 240 Project/srp.m'>srp.m</a></b></td>
                        <td><code>Calculates the solar radiation pressure (SRP) force acting on the spacecraft.</code></td>
                    </tr>
                    <tr>
                        <td><b><a href='https://github.com/JC01010/MAE240/blob/master/exactanalysis/MAE 240 Project/inertial_eom.m'>inertial_eom.m</a></b></td>
                        <td><code>Defines the equations of motion in the inertial frame.</code></td>
                    </tr>
                    <tr>
                        <td><b><a href='https://github.com/JC01010/MAE240/blob/master/exactanalysis/MAE 240 Project/KepEqn.m'>KepEqn.m</a></b></td>
                        <td><code>Solves Kepler's Equation for orbital motion calculations.</code></td>
                    </tr>
                    <tr>
                        <td><b><a href='https://github.com/JC01010/MAE240/blob/master/exactanalysis/MAE 240 Project/tight_subplot.m'>tight_subplot.m</a></b></td>
                        <td><code>Creates subplots with minimal spacing for figure layouts.</code></td>
                    </tr>
                    <tr>
                        <td><b><a href='https://github.com/JC01010/MAE240/blob/master/exactanalysis/MAE 240 Project/zvc.m'>zvc.m</a></b></td>
                        <td><code>Computes zero-velocity curves for the restricted three-body problem.</code></td>
                    </tr>
                    </tr>
                    </table>
                </blockquote>
            </details>
        </blockquote>
    </details>
</details>

---
##  Getting Started

###  Prerequisites

Before getting started with MAE240, ensure your runtime environment meets the following requirements:

- **Programming Language:** MATLAB


###  Installation

Install MAE240 using one of the following methods:

**Build from source:**

1. Clone the MAE240 repository:
```sh
❯ git clone https://github.com/JC01010/MAE240
```

2. Navigate to the project directory:
```sh
❯ cd MAE240
```

##  Contributing

- ** The code for the exact analysis was written by Jingtong Sheng.
- ** The code for the averaged analysis was written by Jason Chang.
- ** Code was provided by Prof. Aaron Rosengren.
