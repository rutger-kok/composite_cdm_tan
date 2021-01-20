# 3D Continuum Damage Mechanics VUMAT for Composite Materials
Explicit material subroutine (VUMAT) implementing a continuum damage mechanics (CDM) model for composite materials in Abaqus (in fixed format Fortran 77).

## Summary
The continuum damage mechanics model implemented in this subroutine is based on the work of several authors. Failure stresses are calculated according to the three-dimensional failure criteria developed by Catalanotti et al. [1]. Damage evolution is adapted from the work of Maimi et al. and Camanho et al. [2][3][4].

The verification directory contains the input files required to run single-element simulations that verify the model correctly predicts the failure stress and energy dissipation for different types of loading. 

The examples directory contains the input files required to run a number of simple simulations illustrating the flexibility of the model.

## Usage
To run a simulation using subroutines your Abaqus installation must be linked with a Fortran compiler and comptaible Visual Studio installation, see: 

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
**Y<sub>T</sub>** =  tensile strength transverse direction  
**Y<sub>C</sub>** = compressive strength transverse direction  
**S<sub>L</sub>** = shear strength  
**α<sub>0</sub>** = initial failure plane angle  
**G<sub>1+</sub>** =  tensile fracture toughness fiber direction  
**G<sub>1-</sub>** =  compressive fracture toughness fiber direction  
**G<sub>2+</sub>** = tensile fracture toughness transverse direction  
**G<sub>2-</sub>** =  compressive fracture toughness transverse direction  
**G<sub>6</sub>** = shear fracture toughness  

## List of Fortran source code
- **composite_cdm.for** : Implementation of CDM model in fixed format Fortran 77 (compatible with all Abaqus installations)
- **composite_cdm.f90** : Implementation of CDM model in free format Fortran 90 (requires modification of Abaqus environment file)

## List of Verification Models  
- **Tension_11** : Pure tension in fiber direction  
- **Tension_22** : Pure tension in transverse direction  
- **Tension_33** : Pure tension in transverse direction  
- **Compression_11** : Pure compression in fiber direction  
- **Compression_22** : Pure compression in fiber direction  
- **Compression_33** : Pure compression in fiber direction  
- **Shear_12** : Shear in the 12 direction  
- **Shear_13** : Shear in the 13 direction  
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

>[4] P.P. Camanho, M.A. Bessa, G. Catalanotti, M. Vogler, R. Rolfes  
>Modeling the inelastic deformation and fracture of polymer composites – Part II: Smeared crack model  
>Mechanics of Materials 59 (2013) 36–49  
>http://dx.doi.org/10.1016/j.mechmat.2012.12.001  
