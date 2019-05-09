# Description 

* This is a program to compute the Harborth constant g(G), that accompanies our article: "On the Harborth constant of `C_3 ⊕ C_3n`

* G is a finite abelian group.

* The Harborth constant g(G) is the smallest integer k such that each set over G whith size at least k, has a subset of size `e=exp(G)` that sums to 0.

* There is no interface but you only have to change this source code to get results for the Harborth constant for the different finite abelian groups.

* This program is valid for any finite abelian group with rank at most four. With the hardware at our disposal it is possible to compute the Harborth constant for finite abelian groups of order up to about 45. The programme does not use in a substantial way that the group is of rank at most four. The restriction is imposed for convenience as we could not treat many groups of rank five anyway.

* The main limiting factor is memory.

* In order to increase the size of accessible groups,currently, we are working on a another version, more efficient based on data compression.

* The group intervenes only in this step `[Initialization]`.

* The subsets of G are represented by a bitmap. 


# Program requirements

Please note that the program has been tested using the following elements:

* The used system is `ubuntu 16.04.10` with `kernel 4.15.0-39`
* The used SHELL is `/bin/bash` version `4.3.48`
* The `gcc` compiler is used with `5.4.0-6 version`

# Sample to run the program

* Compile the program first `make -f makefile`
* Run the program with the compiled file `./sommenulle`

# The user manual

For now, there is no interface but you only have to change this source code to get results for the Harborth constant for the different finite abelian groups.

You have to change the parameters of: `MODULO1, MODULO2, MODULO3 and MODULO4`, according to the finite abelian group of which you are looking for its Harborth constant. 

## Example: 
For the group `C_3 ⊕ C_12`, you have to enter the following parameters:
 
## Initialization: 
```
define RANG 4

define MODULO1 1
define MODULO2 1
	if (RANG>=3)
		define MODULO3 3
	endif
	if (RANG==4)
		define MODULO4 12
	endif
	if (RANG==2)
		define CARDINAL (MODULO1*MODULO2)
		define EXPOSANT MODULO2
	endif
	if (RANG==3)
		define CARDINAL (MODULO1*MODULO2*MODULO3)
		define EXPOSANT MODULO3
	endif
	if(RANG==4)
		define CARDINAL (MODULO1*MODULO2*MODULO3*MODULO4)
		define EXPOSANT MODULO4
	endif 
```
