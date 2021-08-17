NOTE: This CDM routine is working for single element models (i.e. it passes single element verification tests), but is highly unstable for multi-element models. If you spot a bug or any other inconsistency let me know (contact details below).

# 3D Continuum Damage Mechanics VUMAT for Composite Materials
Explicit material subroutine (VUMAT) implementing a continuum damage mechanics (CDM) model for composite materials in Abaqus (in fixed format Fortran 77).

## Summary
The continuum damage mechanics model implemented in this subroutine is based on the work of several authors. Failure stresses are calculated according to the three-dimensional failure criteria developed by Catalanotti et al. [1]. Damage evolution is adapted from the work of Maimi et al. [2][3] and Tan et al. [4][5].

The verification directory contains the input files required to run single-element simulations that verify the model correctly predicts the failure stress and energy dissipation for different types of loading. 

The examples directory contains the input files required to run a number of simple simulations illustrating the flexibility of the model.

## Usage
To run a simulation using subroutines your Abaqus installation must be linked with a Fortran compiler and compatible Visual Studio installation, see: 

The model requires the following material properties to be defined in the simulation input (.inp) file:  
**E<sub>11</sub>** = elastic modulus in the fiber direction  
**E<sub>22</sub>**= elastic modulus transverse direction (in-plane)  
**E<sub>33</sub>** = elastic modulus transverse direction (out-of-plane)  
**ν<sub>12</sub>** = Poisson's ratio 12 direction  
**ν<sub>13</sub>** = Poisson's ratio 13 direction  
**ν<sub>23</sub>** = Poisson's ratio 23 direction  
**G<sub>12</sub>** = shear modulus 12 direction  
**G<sub>13</sub>** = shear modulus 13 direction  
**G<sub>23</sub>** = shear modulus 23 direction  
**X<sub>T</sub>** = tensile strength fiber direction  
**X<sub>C</sub>** = compressive strength fiber direction   
**Y<sub>T</sub><sup>is</sup>** =  in-situ tensile strength transverse direction  
**Y<sub>C</sub><sup>is</sup>** = in-situ compressive strength transverse direction  
**S<sub>L</sub><sup>is</sup>** = in-situ longitudinal shear strength  
**S<sub>T</sub><sup>is</sup>** = in-situ transverse shear strength  
**η<sub>L</sub>** = shear friction coefficient longitudinal direction  
**α<sub>0</sub>** = failure plane angle pure transverse compression  
**G<sub>1+</sub>** =  tensile fracture toughness fiber direction  
**G<sub>1-</sub>** =  compressive fracture toughness fiber direction  
**G<sub>2+</sub>** = tensile fracture toughness transverse direction  
**G<sub>6</sub>** = shear fracture toughness  

## List of Fortran source code
- **composite_cdm.for** : Implementation of CDM model in fixed format Fortran 77 (compatible with all Abaqus installations)  

## List of Verification Models  
- **Tension_11** : Pure tension in fiber direction  
- **Tension_22** : Pure tension in transverse direction  
- **Compression_11** : Pure compression in fiber direction  
- **Shear_12** : Shear in the 12 direction  
- **Shear_23** : Shear in the 23 direction  

## List of Example Models  
- **placeholder** : TBD

***
Rutger Kok  
PhD Candidate  
email : rutger.kok@ed.ac.uk  

Institute for Infrastructure and Environment  
University of Edinburgh    
Thomas Bayes Road, King's Buildings, Edinburgh EH9 3FG   
United Kingdom

***
>[1] G. Catalanotti, P.P. Camanho, A.T. Marques  
>Three-dimensional failure criteria for fiber-reinforced laminates   
>Composite Structures 95 (2013) 63–79  
>http://dx.doi.org/10.1016/j.compstruct.2012.07.016  

>[2] P. Maimi, P.P. Camanho, J.A. Mayugo, C.G. Davila  
>A continuum damage model for composite laminates: Part I – Constitutive model  
>Mechanics of Materials 39 (2007) 897–908  
>http://dx.doi.org/10.1016/j.mechmat.2007.03.005  

>[3] P. Maimi, P.P. Camanho, J.A. Mayugo, C.G. Davila  
>A continuum damage model for composite laminates: Part II – Computational implementation and validation  
>Mechanics of Materials 39 (2007) 909–919  
>http://dx.doi.org/10.1016/j.mechmat.2007.03.006  

>[4] W. Tan, B. G. Falzon, L. N. S. Chiu, and M. Price  
>Predicting low velocity impact damage and Compression-After-Impact (CAI) behaviour of composite laminates  
>Composites Part A 71 (2015) 212–226.  
>http://doi.org/10.1016/j.compositesa.2015.01.025  

>[5] B. G. Falzon, H. Liu, and W. Tan  
>Comment on ‘A tensorial based progressive damage model for fibre reinforced polymers’  
>Composite Structures 176 (2017) 877–882.  
>http://doi.org/10.1016/j.compstruct.2017.06.011
